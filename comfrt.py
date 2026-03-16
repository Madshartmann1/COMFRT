#!/usr/bin/env python3
"""
Competitive Mapping Filter
Recovers reads from a BAM aligned to a combined reference (target + non-target
scaffolds).  Single-pass over the input BAM using pysam:

  Unique reads   MAPQ > threshold on a target scaffold
                 → written directly to a new BAM

  Ambiguous reads  MAPQ = 0 and touches a target scaffold (primary or XA tag)
                 → classified by comparing best NM across target vs non-target
                 → recovered reads (target equally good or better) written as FASTQ
                 → discarded reads (non-target strictly better) dropped

No temp files.  Only dependency beyond pysam is samtools for indexing the
output BAM (samtools index).

Originally by Kirstine Tersbøl Melsen, s215096
Generalized and rewritten with pysam
"""

import argparse
import contextlib
import gzip
import os
import subprocess
import sys

import pysam


# ============================================================================
# CLASSIFICATION HELPERS
# ============================================================================

def parse_xa(xa_value):
    """
    Parse XA tag payload (the value returned by pysam get_tag('XA'), without
    the 'XA:Z:' prefix) into a list of {'scaffold': ..., 'mismatches': ...}.
    """
    entries = [e for e in xa_value.rstrip(';').split(';') if e]
    out = []
    for e in entries:
        parts = e.split(',')
        if len(parts) != 4:
            continue
        try:
            out.append({'scaffold': parts[0], 'mismatches': int(parts[3])})
        except ValueError:
            continue
    return out


def touches_target(read, target_scaffolds):
    """True if read maps to a target scaffold in its primary or any XA alignment."""
    if read.reference_name in target_scaffolds:
        return True
    if read.has_tag('XA'):
        for entry in parse_xa(read.get_tag('XA')):
            if entry['scaffold'] in target_scaffolds:
                return True
    return False


def classify(read, target_scaffolds):
    """
    Classify a MAPQ=0 read that touches a target scaffold.

    Uses NM (edit distance: mismatches + indels) to compare alignments.
    AS (alignment score) would be more accurate — it reflects BWA's actual
    gap-open/extend penalties — but the XA tag only reports NM for secondary
    alignments, so NM is used for consistency across all alignment comparisons.

    Returns:
        'keep'    — target alignment is the best or equal best
        'discard' — non-target alignment is strictly better
        None      — no NM tag, skip
    """
    if not read.has_tag('NM'):
        return None

    primary_nm = read.get_tag('NM')
    sec_list = parse_xa(read.get_tag('XA')) if read.has_tag('XA') else []

    all_alignments = [{'scaffold': read.reference_name, 'mismatches': primary_nm}]
    all_alignments.extend(sec_list)

    target_nms    = [a['mismatches'] for a in all_alignments if a['scaffold'] in target_scaffolds]
    nontarget_nms = [a['mismatches'] for a in all_alignments if a['scaffold'] not in target_scaffolds]

    if not target_nms:
        return 'discard'
    if not nontarget_nms:
        return 'keep'
    min_t  = min(target_nms)
    min_nt = min(nontarget_nms)
    if min_t < min_nt:
        return 'keep'
    elif min_t == min_nt:
        return 'tie'
    else:
        return 'discard'


def read_to_fastq_str(read):
    """Convert a pysam AlignedSegment to a FASTQ-formatted string."""
    seq = read.query_sequence
    if seq is None:
        return None
    quals = read.query_qualities
    qual_str = ''.join(chr(q + 33) for q in quals) if quals is not None else 'I' * len(seq)
    return f"@{read.query_name}\n{seq}\n+\n{qual_str}\n"


# ============================================================================
# REFERENCE LOADING
# ============================================================================

def _load_scaffold_set(path, label='Scaffolds file'):
    if not os.path.exists(path):
        print(f"ERROR: {label} not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        scaffolds = set(line.strip() for line in f if line.strip())
    if not scaffolds:
        print(f"ERROR: No scaffold names found in {path}", file=sys.stderr)
        sys.exit(1)
    return scaffolds


def load_single_ref(scaffolds_file, name):
    """Return a one-entry refs list from a plain scaffold names file."""
    return [{'name': name, 'scaffolds': _load_scaffold_set(scaffolds_file)}]


def load_refs_file(refs_file):
    """
    Parse a two-column references file:
        /path/to/scaffolds.txt    nickname

    Lines starting with # are ignored.
    Returns list of {'name': str, 'scaffolds': set}.
    """
    refs = []
    with open(refs_file) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                print(
                    f"ERROR: {refs_file}:{lineno}: expected 2 columns "
                    f"(scaffold_file  nickname), got: {line!r}",
                    file=sys.stderr,
                )
                sys.exit(1)
            scaffolds_path, nickname = parts[0], parts[1]
            scaffolds = _load_scaffold_set(scaffolds_path, label=f"scaffold file '{nickname}'")
            refs.append({'name': nickname, 'scaffolds': scaffolds})
    if not refs:
        print(f"ERROR: No entries found in {refs_file}", file=sys.stderr)
        sys.exit(1)
    return refs


# ============================================================================
# PIPELINE
# ============================================================================

def run_pipeline(bam, refs, outdir, stats_only, mapq, threads, samtools_path):
    """
    Single pass over `bam`, classifying reads against all scaffold sets simultaneously.

    refs: list of {'name': str, 'scaffolds': set}
    """
    n = len(refs)
    os.makedirs(outdir, exist_ok=True)

    # Build per-ref output paths and stats containers
    for ref in refs:
        # Full output → subdir per ref when multi-ref.  Stats-only → flat.
        if n > 1 and not stats_only:
            ref_dir = os.path.join(outdir, ref['name'])
            os.makedirs(ref_dir, exist_ok=True)
        else:
            ref_dir = outdir

        ref['summary_file'] = os.path.join(ref_dir, f"{ref['name']}_summary.txt")

        if not stats_only:
            ref['unique_bam']   = os.path.join(ref_dir, f"{ref['name']}_unique.bam")
            ref['r1_out']       = os.path.join(ref_dir, f"{ref['name']}_recovered_R1.fq.gz")
            ref['r2_out']       = os.path.join(ref_dir, f"{ref['name']}_recovered_R2.fq.gz")
            ref['merged_out']   = os.path.join(ref_dir, f"{ref['name']}_recovered_merged.fq.gz")

        ref['stats'] = {
            'unique':      0,
            'not_target':  0,
            'target_only': 0,
            'recovered':   0,
            'ties':        0,
            'discarded':   0,
            'skipped_nm':  0,
        }
        ref['r1_count'] = ref['r2_count'] = ref['merged_count'] = 0

    total_seen = 0
    print(f"\n[1/2] Scanning BAM ({n} reference{'s' if n > 1 else ''})...")

    with contextlib.ExitStack() as stack:
        inbam = stack.enter_context(pysam.AlignmentFile(bam, 'rb', threads=threads))

        if not stats_only:
            for ref in refs:
                ref['bam_out'] = stack.enter_context(
                    pysam.AlignmentFile(ref['unique_bam'], 'wb', template=inbam, threads=threads)
                )
                ref['f1']   = stack.enter_context(gzip.open(ref['r1_out'],     'wt'))
                ref['f2']   = stack.enter_context(gzip.open(ref['r2_out'],     'wt'))
                ref['fm']   = stack.enter_context(gzip.open(ref['merged_out'], 'wt'))

        for read in inbam:
            if read.is_unmapped:
                continue
            total_seen += 1

            for ref in refs:
                target_scaffolds = ref['scaffolds']
                st = ref['stats']

                # Unique: MAPQ > threshold AND primary alignment on a target scaffold
                if read.mapping_quality > mapq and read.reference_name in target_scaffolds:
                    if not stats_only:
                        ref['bam_out'].write(read)
                    st['unique'] += 1
                    continue  # skip ambiguous classification for this ref

                # Ambiguous: MAPQ = 0
                if read.mapping_quality == 0:
                    if not touches_target(read, target_scaffolds):
                        st['not_target'] += 1
                        continue

                    result = classify(read, target_scaffolds)

                    if result is None:
                        st['skipped_nm'] += 1
                        continue

                    if result in ('keep', 'tie'):
                        if not stats_only:
                            rec = read_to_fastq_str(read)
                            if rec is not None:
                                if read.is_paired:
                                    if read.is_read1:
                                        ref['f1'].write(rec); ref['r1_count'] += 1
                                    elif read.is_read2:
                                        ref['f2'].write(rec); ref['r2_count'] += 1
                                    else:
                                        ref['fm'].write(rec); ref['merged_count'] += 1
                                else:
                                    ref['fm'].write(rec); ref['merged_count'] += 1

                        if result == 'tie':
                            st['ties'] += 1
                        else:
                            xa = read.get_tag('XA') if read.has_tag('XA') else None
                            sec_list = parse_xa(xa) if xa else []
                            if any(e['scaffold'] not in target_scaffolds for e in sec_list):
                                st['recovered'] += 1
                            else:
                                st['target_only'] += 1
                    else:
                        st['discarded'] += 1

    # Index unique BAMs
    if not stats_only:
        print(f"[2/2] Indexing {n} unique BAM{'s' if n > 1 else ''}...")
        for ref in refs:
            subprocess.run([samtools_path, 'index', ref['unique_bam']], check=True)

    # Print and write summaries
    for ref in refs:
        st = ref['stats']
        r1 = ref['r1_count']
        r2 = ref['r2_count']
        rm = ref['merged_count']
        touching = (st['unique'] + st['target_only'] + st['recovered']
                    + st['ties'] + st['discarded'] + st['skipped_nm'])
        pct = 100.0 * touching / total_seen if total_seen else 0.0

        print(f"\n{'=' * 45}")
        print(f"Reference: {ref['name']}")
        print(f"  Total mapped reads:            {total_seen}")
        print(f"  Touch target scaffold:         {touching}  ({pct:.1f}%)")
        print(f"  Unique target reads:           {st['unique']}")
        print(f"  Ambiguous \u2192 target-only:       {st['target_only']}")
        print(f"  Ambiguous \u2192 recovered:         {st['recovered']}")
        print(f"  Ambiguous \u2192 ties (kept):       {st['ties']}")
        print(f"  Ambiguous \u2192 discarded:         {st['discarded']}")
        print(f"  Ambiguous \u2192 not touching:      {st['not_target']}")
        if st['skipped_nm']:
            print(f"  Skipped (no NM tag):           {st['skipped_nm']}")
        if not stats_only:
            print(f"  Recovered FASTQ: R1={r1}  R2={r2}  Merged={rm}  Total={r1+r2+rm}")

        with open(ref['summary_file'], 'w') as f:
            f.write("Competitive Mapping Filter - Summary\n")
            f.write("=" * 40 + "\n")
            f.write(f"Reference:           {ref['name']}\n")
            f.write(f"Input BAM:           {bam}\n")
            f.write(f"MAPQ threshold:      {mapq}\n")
            f.write(f"Stats only:          {'yes' if stats_only else 'no'}\n\n")
            f.write(f"Total mapped reads:          {total_seen}\n")
            f.write(f"Touch target scaffold:       {touching}  ({pct:.1f}%)\n\n")
            f.write(f"Unique target reads (MAPQ > {mapq}): {st['unique']}\n\n")
            f.write(f"Ambiguous reads (MAPQ=0, target-touching):\n")
            f.write(f"  Target-only (no competing alignments): {st['target_only']}\n")
            f.write(f"  Recovered (target strictly better NM): {st['recovered']}\n")
            f.write(f"  Ties (equal NM, kept):                 {st['ties']}\n")
            f.write(f"  Discarded (non-target better):         {st['discarded']}\n")
            f.write(f"  Not touching target:                   {st['not_target']}\n")
            if st['skipped_nm']:
                f.write(f"  Skipped (no NM tag):                   {st['skipped_nm']}\n")
            if not stats_only:
                f.write(f"\nRecovered FASTQ:\n")
                f.write(f"  R1 reads:        {r1}\n")
                f.write(f"  R2 reads:        {r2}\n")
                f.write(f"  Merged/SE reads: {rm}\n")
                f.write(f"  Total:           {r1+r2+rm}\n")

    print(f"\n{'=' * 45}")
    print(f"Done. Output in: {outdir}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            'Competitive Mapping Filter\n'
            '==========================\n'
            'Recovers mitochondrial (or other target-scaffold) reads from a BAM\n'
            'that was aligned to a combined reference containing both target and\n'
            'non-target scaffolds.  Reads are split into two pools based on MAPQ:\n'
            '\n'
            '  UNIQUE   MAPQ > threshold AND maps to a target scaffold\n'
            '           → written directly to a new BAM\n'
            '\n'
            '  AMBIGUOUS  MAPQ = 0 AND touches a target scaffold\n'
            '             (primary alignment or any BWA XA secondary entry)\n'
            '           → classified by comparing NM (edit distance) of the best\n'
            '             target alignment against the best non-target alignment:\n'
            '               target NM <  non-target NM  → recovered  (target better)\n'
            '               target NM == non-target NM  → tie        (kept)\n'
            '               target NM >  non-target NM  → discarded  (non-target better)\n'
            '           → recovered / tie reads written as gzip FASTQ\n'
            '\n'
            'Run modes:\n'
            '  Single-reference  -s scaffolds.txt -n name    (one scaffold set)\n'
            '  Multi-reference   -r references.tsv           (two-column: path  nickname)\n'
            '                    The BAM is read ONCE for all scaffold sets.\n'
            '\n'
            'IMPORTANT: The input BAM must have been aligned with BWA and must\n'
            'retain MAPQ=0 reads (do not filter with -q before running this script).\n'
            'The BAM must be coordinate-sorted and indexed.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  # Single reference, full output:\n'
            '  %(prog)s -i aligned.bam -s mt_scaffolds.txt -o output/ -n quagga -q 1 -t 16\n'
            '\n'
            '  # Multi-reference, full output (one subfolder per reference):\n'
            '  %(prog)s -i aligned.bam -r refs.tsv -o output/\n'
            '\n'
            '  # Multi-reference, stats only (flat summary txts, no BAM/FASTQ):\n'
            '  %(prog)s -i aligned.bam -r refs.tsv -o output/ --stats-only\n'
            '\n'
            'refs.tsv format (tab or space separated, # = comment):\n'
            '  /path/to/mt_scaffolds.txt    mitochondria\n'
            '  /path/to/y_scaffolds.txt     y_chromosome\n'
            '\n'
            'Full output (per reference, subfolders when multi-ref):\n'
            '  <name>_unique.bam              Unique target reads,\n'
            '                                 indexed automatically.\n'
            '  <name>_recovered_R1.fq.gz      Recovered ambiguous reads — R1 of a pair.\n'
            '  <name>_recovered_R2.fq.gz      Recovered ambiguous reads — R2 of a pair.\n'
            '  <name>_recovered_merged.fq.gz  Recovered ambiguous reads — merged or SE.\n'
            '  <name>_summary.txt             Run statistics.\n'
            '\n'
            'Stats-only output (flat, directly in output dir):\n'
            '  <name>_summary.txt             One per reference.\n'
        ),
    )

    parser.add_argument(
        '-i', '--input', required=True, metavar='BAM',
        help=(
            'Input BAM aligned to a COMBINED reference (target + non-target scaffolds). '
            'Must be coordinate-sorted, indexed (.bai), and must RETAIN MAPQ=0 reads — '
            'do not pre-filter with samtools view -q. Required.'
        ),
    )

    scaffold_group = parser.add_mutually_exclusive_group(required=True)
    scaffold_group.add_argument(
        '-s', '--target-scaffolds', metavar='FILE',
        help=(
            'Plain-text file listing target scaffold names, one per line '
            '(e.g. mitochondrial chromosome names). '
            'Names must match sequence names in the BAM header exactly. '
            'Use -n to set the output name prefix. '
            'Mutually exclusive with -r.'
        ),
    )
    scaffold_group.add_argument(
        '-r', '--references', metavar='FILE',
        help=(
            'Two-column file: scaffold_file_path  nickname. '
            'Each row defines one scaffold set to run in a single BAM pass. '
            'Lines starting with # are ignored. '
            'Mutually exclusive with -s.'
        ),
    )

    parser.add_argument(
        '-o', '--output-dir', required=True, metavar='DIR',
        help=(
            'Output directory. With -r and full output, a subfolder is created '
            'per reference inside this directory. Created if it does not exist. Required.'
        ),
    )
    parser.add_argument(
        '-n', '--sample-name', default='sample', metavar='STR',
        help='Output filename prefix when using -s (default: sample). Ignored with -r.',
    )
    parser.add_argument(
        '--stats-only', action='store_true',
        help=(
            'Skip all BAM/FASTQ output; write only per-reference summary.txt files. '
            'Useful for a quick stats run before committing to full output.'
        ),
    )
    parser.add_argument(
        '-q', '--mapq-threshold', type=int, default=0, metavar='INT',
        help=(
            'MAPQ threshold for classifying a read as "unique". '
            'Reads with MAPQ STRICTLY GREATER THAN this value that map to a target '
            'scaffold are written to the unique BAM. '
            'All MAPQ=0 target-touching reads go through competitive NM filtering. '
            'Default: 0 (MAPQ ≥ 1 → unique). Common choices: 0, 1, 20, 25.'
        ),
    )
    parser.add_argument(
        '-t', '--threads', type=int, default=4, metavar='INT',
        help='Threads for pysam and samtools (default: 4).',
    )
    parser.add_argument(
        '--samtools', default='samtools', metavar='PATH',
        help='Path to samtools, used only for indexing output BAMs (default: samtools).',
    )

    args = parser.parse_args()

    bam    = os.path.abspath(args.input)
    outdir = os.path.abspath(args.output_dir)
    mapq   = args.mapq_threshold

    # --- Validate BAM ---
    if not os.path.exists(bam):
        print(f"ERROR: Input BAM not found: {bam}", file=sys.stderr)
        sys.exit(1)
    if not bam.endswith('.bam'):
        print("ERROR: Input must be a .bam file", file=sys.stderr)
        sys.exit(1)
    has_index = os.path.exists(bam + '.bai') or os.path.exists(bam.replace('.bam', '.bai'))
    if not has_index:
        print(f"ERROR: BAM index not found. Run: {args.samtools} index {bam}", file=sys.stderr)
        sys.exit(1)

    # --- Load scaffold sets ---
    if args.target_scaffolds:
        refs = load_single_ref(os.path.abspath(args.target_scaffolds), args.sample_name)
    else:
        refs = load_refs_file(os.path.abspath(args.references))

    mode = 'stats-only' if args.stats_only else 'full output'
    print(f"Input BAM:      {bam}")
    print(f"References:     {len(refs)} ({',  '.join(r['name'] for r in refs)})")
    print(f"MAPQ threshold: MAPQ > {mapq} → unique")
    print(f"Mode:           {mode}")

    run_pipeline(bam, refs, outdir, args.stats_only, mapq, args.threads, args.samtools)


if __name__ == '__main__':
    main()
