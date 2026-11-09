#!/usr/bin/env python3
"""
COMFRT: COmpetitively Mapped Fragment Recovery Tool (Competitive Mapping Filter)

Recovers reads from a BAM aligned to a combined reference (target + non-target
scaffolds).  Single pass over the input BAM using pysam.  Only primary
alignments are classified; secondary and supplementary records are skipped.

  Unique reads     MAPQ > threshold and primary alignment on a target scaffold
                   → written directly to a new BAM with a target-only header

  Ambiguous reads  MAPQ <= threshold and touches a target scaffold
                   (primary alignment or XA tag)
                   → classified by comparing best NM across target vs non-target
                   → kept reads (target strictly better, or equal unless
                     --drop-ties) written as FASTQ for remapping to the target
                   → discarded reads (non-target strictly better) dropped

No temp files.  samtools is only used to index the output BAM.

Originally by Kirstine Tersbøl Melsen
Generalized and rewritten with pysam
"""

import argparse
import contextlib
import gzip
import os
import shlex
import subprocess
import sys

import pysam

__version__ = '1.0.0'

FLAG_PROPER_PAIR = 0x2
FLAG_MATE_UNMAPPED = 0x8
FLAG_MATE_REVERSE = 0x20


# ============================================================================
# CLASSIFICATION HELPERS
# ============================================================================

def parse_xa(xa_value):
    """
    Parse XA tag payload (the value returned by pysam get_tag('XA'), without
    the 'XA:Z:' prefix) into a list of {'scaffold': ..., 'mismatches': ...}.
    """
    xa_entries = [entry_text for entry_text in xa_value.rstrip(';').split(';') if entry_text]
    parsed_alignments = []
    for entry_text in xa_entries:
        entry_fields = entry_text.split(',')
        if len(entry_fields) != 4:
            continue
        try:
            parsed_alignments.append({'scaffold': entry_fields[0], 'mismatches': int(entry_fields[3])})
        except ValueError:
            continue
    return parsed_alignments


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
    Classify an ambiguous read (MAPQ <= threshold) that touches a target scaffold.

    Uses NM (edit distance: mismatches + indels) to compare alignments.
    AS (alignment score) would be more accurate — it reflects BWA's actual
    gap-open/extend penalties — but the XA tag only reports NM for secondary
    alignments, so NM is used for consistency across all alignment comparisons.

    Returns:
        'keep'    — target alignment is strictly best, or no non-target
                    alignment is listed (primary or XA)
        'tie'     — best target NM equals best non-target NM
        'discard' — non-target alignment is strictly better
        None      — no NM tag, skip
    """
    if not read.has_tag('NM'):
        return None

    primary_nm = read.get_tag('NM')
    secondary_alignments = parse_xa(read.get_tag('XA')) if read.has_tag('XA') else []

    all_alignments = [{'scaffold': read.reference_name, 'mismatches': primary_nm}]
    all_alignments.extend(secondary_alignments)

    target_mismatch_counts = [alignment['mismatches'] for alignment in all_alignments if alignment['scaffold'] in target_scaffolds]
    nontarget_mismatch_counts = [alignment['mismatches'] for alignment in all_alignments if alignment['scaffold'] not in target_scaffolds]

    if not target_mismatch_counts:
        return 'discard'
    if not nontarget_mismatch_counts:
        return 'keep'
    minimum_target_mismatch = min(target_mismatch_counts)
    minimum_nontarget_mismatch = min(nontarget_mismatch_counts)
    if minimum_target_mismatch < minimum_nontarget_mismatch:
        return 'keep'
    elif minimum_target_mismatch == minimum_nontarget_mismatch:
        return 'tie'
    else:
        return 'discard'


_REVCOMP_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def _reverse_complement(seq):
    return seq.translate(_REVCOMP_TABLE)[::-1]


def read_to_fastq_str(read):
    """
    Convert a pysam AlignedSegment to a FASTQ-formatted string.

    pysam's query_sequence/query_qualities are in BAM/alignment orientation:
    for a reverse-strand read that's the reverse complement of what the
    sequencer actually produced. Flip it back so the exported FASTQ matches
    the original read orientation (required for correct remapping/merging).
    """
    seq = read.query_sequence
    if seq is None:
        return None
    query_qualities = read.query_qualities
    if read.is_reverse:
        seq = _reverse_complement(seq)
        if query_qualities is not None:
            query_qualities = query_qualities[::-1]
    qual_str = ''.join(chr(quality_score + 33) for quality_score in query_qualities) if query_qualities is not None else 'I' * len(seq)
    return f"@{read.query_name}\n{seq}\n+\n{qual_str}\n"


def write_recovered_read_to_bucket(read, recovered_fastq_record, reference_config):
    """
    Write one recovered read to the correct FASTQ bucket.

    Only true read pairs are written to R1/R2. If a mate is missing, the
    singleton read is written to merged/SE during final flush.
    """
    if not read.is_paired:
        reference_config['fm'].write(recovered_fastq_record)
        reference_config['merged_count'] += 1
        return

    if not (read.is_read1 or read.is_read2):
        reference_config['fm'].write(recovered_fastq_record)
        reference_config['merged_count'] += 1
        return

    pending_pairs_by_query_name = reference_config['pending_pairs_by_query_name']
    query_name = read.query_name
    pair_entry = pending_pairs_by_query_name.setdefault(query_name, {'r1': None, 'r2': None})

    if read.is_read1:
        if pair_entry['r1'] is not None:
            reference_config['fm'].write(pair_entry['r1'])
            reference_config['merged_count'] += 1
        pair_entry['r1'] = recovered_fastq_record
    elif read.is_read2:
        if pair_entry['r2'] is not None:
            reference_config['fm'].write(pair_entry['r2'])
            reference_config['merged_count'] += 1
        pair_entry['r2'] = recovered_fastq_record

    if pair_entry['r1'] is not None and pair_entry['r2'] is not None:
        reference_config['f1'].write(pair_entry['r1'])
        reference_config['f2'].write(pair_entry['r2'])
        reference_config['r1_count'] += 1
        reference_config['r2_count'] += 1
        del pending_pairs_by_query_name[query_name]


def flush_unpaired_recovered_reads(reference_config):
    """Move any leftover recovered singleton mates into merged/SE output."""
    pending_pairs_by_query_name = reference_config['pending_pairs_by_query_name']
    for pair_entry in pending_pairs_by_query_name.values():
        if pair_entry['r1'] is not None:
            reference_config['fm'].write(pair_entry['r1'])
            reference_config['merged_count'] += 1
        if pair_entry['r2'] is not None:
            reference_config['fm'].write(pair_entry['r2'])
            reference_config['merged_count'] += 1
    pending_pairs_by_query_name.clear()


# ============================================================================
# UNIQUE BAM OUTPUT
# ============================================================================

def build_target_header(input_header, target_scaffolds, command_line):
    """
    Build the unique BAM header: only the target @SQ lines (input order kept),
    @HD/@RG/@PG/@CO copied from the input, and a @PG line for this run.

    Returns (header_dict, new_tid_by_old_tid).  The list maps every input
    reference id to its id in the new header, or -1 for non-target scaffolds.
    """
    header_dict = input_header.to_dict()
    target_sequences = [
        sequence_entry for sequence_entry in header_dict.get('SQ', [])
        if sequence_entry.get('SN') in target_scaffolds
    ]
    new_tid_by_name = {sequence_entry['SN']: new_tid for new_tid, sequence_entry in enumerate(target_sequences)}
    new_tid_by_old_tid = [new_tid_by_name.get(name, -1) for name in input_header.references]
    header_dict['SQ'] = target_sequences

    program_entries = list(header_dict.get('PG', []))
    existing_program_ids = {program_entry.get('ID') for program_entry in program_entries}
    referenced_program_ids = {program_entry.get('PP') for program_entry in program_entries}
    chain_end_ids = [
        program_entry['ID'] for program_entry in program_entries
        if 'ID' in program_entry and program_entry['ID'] not in referenced_program_ids
    ]
    program_id, suffix = 'comfrt', 1
    while program_id in existing_program_ids:
        program_id = f'comfrt.{suffix}'
        suffix += 1
    comfrt_program_entry = {'ID': program_id, 'PN': 'comfrt', 'VN': __version__, 'CL': command_line}
    if chain_end_ids:
        comfrt_program_entry['PP'] = chain_end_ids[-1]
    header_dict['PG'] = program_entries + [comfrt_program_entry]
    return header_dict, new_tid_by_old_tid


def write_unique_read(read, output_bam, new_tid_by_old_tid):
    """
    Write a unique target read to a BAM whose header holds only target scaffolds.

    Reference ids are translated to the new header.  A mate on a non-target
    scaffold has no @SQ line there, so that mate is recorded as unmapped
    (RNEXT '*', PNEXT 0, TLEN 0, mate-unmapped flag, MC/MQ tags removed).
    The input read is left unchanged for the other references in the pass.
    """
    old_tid = read.reference_id
    old_mate_tid = read.next_reference_id
    new_tid = new_tid_by_old_tid[old_tid]
    new_mate_tid = new_tid_by_old_tid[old_mate_tid] if old_mate_tid >= 0 else -1

    if old_mate_tid >= 0 and new_mate_tid < 0:
        read_fields = read.to_dict()
        flag = read.flag
        if read.is_paired:
            flag = (flag | FLAG_MATE_UNMAPPED) & ~(FLAG_PROPER_PAIR | FLAG_MATE_REVERSE)
        read_fields.update(flag=str(flag), next_ref_name='*', next_ref_pos='0', length='0')
        read_fields['tags'] = [tag for tag in read_fields['tags'] if not tag.startswith(('MC:', 'MQ:'))]
        output_bam.write(pysam.AlignedSegment.from_dict(read_fields, output_bam.header))
        return

    read.reference_id = new_tid
    read.next_reference_id = new_mate_tid
    try:
        output_bam.write(read)
    finally:
        read.reference_id = old_tid
        read.next_reference_id = old_mate_tid


# ============================================================================
# REFERENCE LOADING AND CHECKS
# ============================================================================

def _load_scaffold_set(path, label='Scaffolds file'):
    """
    Read scaffold names, one per line.  Only the first whitespace-separated
    field is used and a leading '>' is dropped, so a samtools faidx .fai file
    or FASTA header lines work as well as a plain list of names.
    """
    if not os.path.exists(path):
        print(f"ERROR: {label} not found: {path}", file=sys.stderr)
        sys.exit(1)
    scaffolds = set()
    with open(path) as input_file_handle:
        for line in input_file_handle:
            line_fields = line.split()
            if line_fields and line_fields[0].lstrip('>'):
                scaffolds.add(line_fields[0].lstrip('>'))
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

    Lines starting with # are ignored.  Relative paths are resolved from the
    current working directory.
    Returns list of {'name': str, 'scaffolds': set}.
    """
    references = []
    with open(refs_file) as input_file_handle:
        for line_number, raw_line in enumerate(input_file_handle, 1):
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            columns = line.split()
            if len(columns) < 2:
                print(
                    f"ERROR: {refs_file}:{line_number}: expected 2 columns "
                    f"(scaffold_file  nickname), got: {line!r}",
                    file=sys.stderr,
                )
                sys.exit(1)
            scaffolds_path, nickname = columns[0], columns[1]
            scaffolds = _load_scaffold_set(scaffolds_path, label=f"scaffold file '{nickname}'")
            references.append({'name': nickname, 'scaffolds': scaffolds})
    if not references:
        print(f"ERROR: No entries found in {refs_file}", file=sys.stderr)
        sys.exit(1)
    return references


def check_reference_names(refs):
    """Output names become file and folder names: they must be unique and contain no '/'."""
    seen_names = set()
    for reference_config in refs:
        name = reference_config['name']
        if not name or os.sep in name or name in ('.', '..'):
            print(f"ERROR: Invalid output name {name!r} (used as a file/folder name)", file=sys.stderr)
            sys.exit(1)
        if name in seen_names:
            print(f"ERROR: Nickname {name!r} is used more than once; outputs would overwrite each other", file=sys.stderr)
            sys.exit(1)
        seen_names.add(name)


def check_scaffolds_in_header(refs, header_scaffolds):
    """Stop if no scaffold of a reference is in the BAM header; warn if some are missing."""
    for reference_config in refs:
        scaffolds = reference_config['scaffolds']
        missing_scaffolds = sorted(scaffolds - header_scaffolds)
        if len(missing_scaffolds) == len(scaffolds):
            print(
                f"ERROR: None of the {len(scaffolds)} scaffold names for '{reference_config['name']}' "
                f"are in the BAM header (first: {missing_scaffolds[0]!r}). Names must match the "
                f"@SQ SN: names exactly.",
                file=sys.stderr,
            )
            sys.exit(1)
        if missing_scaffolds:
            print(
                f"WARNING: {len(missing_scaffolds)} of {len(scaffolds)} scaffold names for "
                f"'{reference_config['name']}' are not in the BAM header and are ignored "
                f"(e.g. {', '.join(missing_scaffolds[:5])})",
                file=sys.stderr,
            )


# ============================================================================
# SUMMARY
# ============================================================================

def summary_lines(reference_config, bam, mapq, stats_only, drop_ties, total_mapped_reads, skipped_non_primary):
    """Lines of the per-reference summary (printed and written to <name>_summary.txt)."""
    reference_stats = reference_config['stats']
    target_touching_ambiguous = (
        reference_stats['recovered']
        + reference_stats['ties']
        + reference_stats['target_only_xa']
        + reference_stats['target_only_no_xa']
        + reference_stats['discarded']
        + reference_stats['skipped_nm']
    )
    total_ambiguous = target_touching_ambiguous + reference_stats['not_target']
    touch_target = reference_stats['unique'] + target_touching_ambiguous
    touch_target_percent = 100.0 * touch_target / total_mapped_reads if total_mapped_reads else 0.0
    tie_handling = 'dropped' if drop_ties else 'kept'

    def row(label, value):
        return f"{label:<48}{value}"

    lines = [
        "COMFRT (Competitive Mapping Filter) - Summary",
        "=" * 48,
        row("COMFRT version:", __version__),
        row("Reference:", reference_config['name']),
        row("Input BAM:", bam),
        row("MAPQ threshold:", f"{mapq}  (unique: MAPQ > {mapq}; ambiguous: MAPQ <= {mapq})"),
        row("Ties:", tie_handling),
        row("Stats only:", 'yes' if stats_only else 'no'),
        "",
        row("Total mapped reads (primary alignments):", total_mapped_reads),
        row(f"  Unique target reads (MAPQ > {mapq}):", reference_stats['unique']),
        row(f"  Non-target reads (MAPQ > {mapq}):", reference_stats['nontarget_confident']),
        row(f"  Ambiguous reads (MAPQ <= {mapq}):", total_ambiguous),
        row("    Not touching target:", reference_stats['not_target']),
        row("    Target-touching:", target_touching_ambiguous),
    ]
    if mapq > 0:
        lines.append(row(f"      of which MAPQ 1-{mapq}:", reference_stats['touching_mapq_above_0']))
    if skipped_non_primary:
        lines.append(row("Secondary/supplementary records (skipped):", skipped_non_primary))
    lines += [
        "",
        row("Touch target (unique + target-touching):", f"{touch_target}  ({touch_target_percent:.1f}%)"),
        "",
        "Target-touching ambiguous breakdown:",
        row("  Recovered (target strictly better NM):", reference_stats['recovered']),
        row(f"  Ties (equal NM, {tie_handling}):", reference_stats['ties']),
        row("  Target-only, XA lists only target hits:", reference_stats['target_only_xa']),
        row("  Target-only, no XA tag (not compared):", reference_stats['target_only_no_xa']),
        row("  Discarded (non-target strictly better NM):", reference_stats['discarded']),
    ]
    if reference_stats['skipped_nm']:
        lines.append(row("  Skipped (no NM tag):", reference_stats['skipped_nm']))
    if not stats_only:
        recovered_r1_count = reference_config['r1_count']
        recovered_r2_count = reference_config['r2_count']
        recovered_merged_count = reference_config['merged_count']
        lines += [
            "",
            "Recovered FASTQ (for remapping to the target):",
            row("  R1 reads:", recovered_r1_count),
            row("  R2 reads:", recovered_r2_count),
            row("  Merged/SE reads:", recovered_merged_count),
            row("  Total:", recovered_r1_count + recovered_r2_count + recovered_merged_count),
        ]
    return lines


# ============================================================================
# PIPELINE
# ============================================================================

def run_pipeline(bam, refs, outdir, stats_only, mapq, threads, samtools_path,
                 drop_ties=False, command_line='comfrt.py'):
    """
    Single pass over `bam`, classifying reads against all scaffold sets simultaneously.

    refs: list of {'name': str, 'scaffolds': set}
    """
    reference_count = len(refs)
    os.makedirs(outdir, exist_ok=True)

    # -t is a total: half for reading the input BAM, half shared by the output BAMs
    if stats_only:
        input_threads, output_threads = threads, 0
    else:
        input_threads = max(1, threads // 2)
        output_threads = max(1, (threads - input_threads) // reference_count)

    # Build per-ref output paths and stats containers
    for reference_config in refs:
        # Full output → subdir per ref when multi-ref.  Stats-only → flat.
        if reference_count > 1 and not stats_only:
            reference_output_dir = os.path.join(outdir, reference_config['name'])
            os.makedirs(reference_output_dir, exist_ok=True)
        else:
            reference_output_dir = outdir

        reference_config['summary_file'] = os.path.join(
            reference_output_dir,
            f"{reference_config['name']}_summary.txt",
        )

        if not stats_only:
            reference_config['unique_bam'] = os.path.join(reference_output_dir, f"{reference_config['name']}_unique.bam")
            reference_config['r1_out'] = os.path.join(reference_output_dir, f"{reference_config['name']}_recovered_R1.fq.gz")
            reference_config['r2_out'] = os.path.join(reference_output_dir, f"{reference_config['name']}_recovered_R2.fq.gz")
            reference_config['merged_out'] = os.path.join(reference_output_dir, f"{reference_config['name']}_recovered_merged.fq.gz")

        reference_config['stats'] = {
            'unique':                0,
            'nontarget_confident':   0,
            'not_target':            0,
            'recovered':             0,
            'ties':                  0,
            'target_only_xa':        0,
            'target_only_no_xa':     0,
            'discarded':             0,
            'skipped_nm':            0,
            'touching_mapq_above_0': 0,
        }
        reference_config['r1_count'] = 0
        reference_config['r2_count'] = 0
        reference_config['merged_count'] = 0
        reference_config['pending_pairs_by_query_name'] = {}

    total_mapped_reads = 0
    skipped_non_primary = 0
    print(f"\n[1/2] Scanning BAM ({reference_count} reference{'s' if reference_count > 1 else ''})...")

    with contextlib.ExitStack() as stack:
        input_bam = stack.enter_context(pysam.AlignmentFile(bam, 'rb', threads=input_threads))

        if not stats_only:
            for reference_config in refs:
                target_header, reference_config['new_tid_by_old_tid'] = build_target_header(
                    input_bam.header, reference_config['scaffolds'], command_line,
                )
                reference_config['bam_out'] = stack.enter_context(
                    pysam.AlignmentFile(reference_config['unique_bam'], 'wb', header=target_header, threads=output_threads)
                )
                reference_config['f1'] = stack.enter_context(gzip.open(reference_config['r1_out'], 'wt', compresslevel=6))
                reference_config['f2'] = stack.enter_context(gzip.open(reference_config['r2_out'], 'wt', compresslevel=6))
                reference_config['fm'] = stack.enter_context(gzip.open(reference_config['merged_out'], 'wt', compresslevel=6))

        for read in input_bam:
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                skipped_non_primary += 1
                continue
            total_mapped_reads += 1
            primary_scaffold = read.reference_name
            mapping_quality = read.mapping_quality

            for reference_config in refs:
                target_scaffolds = reference_config['scaffolds']
                reference_stats = reference_config['stats']
                primary_on_target = primary_scaffold in target_scaffolds

                # Confidently placed: unique if on target, otherwise not considered
                if mapping_quality > mapq:
                    if primary_on_target:
                        if not stats_only:
                            write_unique_read(read, reference_config['bam_out'], reference_config['new_tid_by_old_tid'])
                        reference_stats['unique'] += 1
                    else:
                        reference_stats['nontarget_confident'] += 1
                    continue

                # Ambiguous: MAPQ <= threshold
                if not touches_target(read, target_scaffolds):
                    reference_stats['not_target'] += 1
                    continue
                if mapping_quality > 0:
                    reference_stats['touching_mapq_above_0'] += 1

                result = classify(read, target_scaffolds)

                if result is None:
                    reference_stats['skipped_nm'] += 1
                    continue
                if result == 'discard':
                    reference_stats['discarded'] += 1
                    continue
                if result == 'tie':
                    reference_stats['ties'] += 1
                    if drop_ties:
                        continue
                else:
                    # 'keep': split by whether a non-target hit was listed at all
                    secondary_alignments = parse_xa(read.get_tag('XA')) if read.has_tag('XA') else []
                    if not primary_on_target or any(
                        alignment['scaffold'] not in target_scaffolds for alignment in secondary_alignments
                    ):
                        reference_stats['recovered'] += 1
                    elif secondary_alignments:
                        reference_stats['target_only_xa'] += 1
                    else:
                        reference_stats['target_only_no_xa'] += 1

                if not stats_only:
                    recovered_fastq_record = read_to_fastq_str(read)
                    if recovered_fastq_record is not None:
                        write_recovered_read_to_bucket(read, recovered_fastq_record, reference_config)

        if not stats_only:
            for reference_config in refs:
                flush_unpaired_recovered_reads(reference_config)

    # Index unique BAMs (written with target-only headers)
    if not stats_only:
        print(f"[2/2] Indexing {reference_count} unique BAM{'s' if reference_count > 1 else ''}...")
        for reference_config in refs:
            subprocess.run([samtools_path, 'index', reference_config['unique_bam']], check=True)

    # Print and write summaries
    for reference_config in refs:
        lines = summary_lines(
            reference_config, bam, mapq, stats_only, drop_ties, total_mapped_reads, skipped_non_primary,
        )
        print("\n" + "\n".join(lines))
        with open(reference_config['summary_file'], 'w') as summary_file_handle:
            summary_file_handle.write("\n".join(lines) + "\n")

    print(f"\n{'=' * 48}")
    print(f"Done. Output in: {outdir}")


# ============================================================================
# MAIN
# ============================================================================

def _non_negative_int(text):
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got {text!r}")
    if value < 0:
        raise argparse.ArgumentTypeError(f"must be 0 or more, got {value}")
    return value


def _positive_int(text):
    value = _non_negative_int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be 1 or more")
    return value


def main():
    parser = argparse.ArgumentParser(
        description=(
            f'COMFRT {__version__}: Competitive Mapping Filter\n'
            '==========================================\n'
            'Recovers mitochondrial (or other target-scaffold) reads from a BAM\n'
            'that was aligned to a combined reference containing both target and\n'
            'non-target scaffolds.  Only primary alignments are used.  Reads are\n'
            'split into two pools based on MAPQ (threshold q, set with -q):\n'
            '\n'
            '  UNIQUE     MAPQ > q AND primary alignment on a target scaffold\n'
            '           → written directly to a new BAM (target-only header)\n'
            '\n'
            '  AMBIGUOUS  MAPQ <= q AND touches a target scaffold\n'
            '             (primary alignment or any BWA XA secondary entry)\n'
            '           → classified by comparing NM (edit distance) of the best\n'
            '             target alignment against the best non-target alignment:\n'
            '               target NM <  non-target NM  → recovered  (target better)\n'
            '               target NM == non-target NM  → tie        (kept unless --drop-ties)\n'
            '               target NM >  non-target NM  → discarded  (non-target better)\n'
            '               no non-target hit listed    → target-only (kept, not compared)\n'
            '           → kept reads written as gzip FASTQ for remapping to the target\n'
            '\n'
            'Run modes:\n'
            '  Single-reference  -s scaffolds.txt -n name    (one scaffold set)\n'
            '  Multi-reference   -r references.tsv           (two-column: path  nickname)\n'
            '                    The BAM is read ONCE for all scaffold sets.\n'
            '\n'
            'IMPORTANT: The input BAM must have been aligned with BWA and must\n'
            'retain MAPQ=0 reads (do not filter with -q before running this script).\n'
            'For full output it must be coordinate-sorted; an index is not needed.\n'
            'BWA leaves out the XA tag when a read has more hits than samse/sampe -n\n'
            '(default 3); those reads cannot be compared (see the README).'
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
            '  <name>_unique.bam              Unique target reads, target-only header,\n'
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

    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument(
        '-i', '--input', required=True, metavar='BAM',
        help=(
            'Input BAM aligned with BWA to a COMBINED reference (target + non-target '
            'scaffolds). Must RETAIN MAPQ=0 reads — do not pre-filter with samtools view -q. '
            'Must be coordinate-sorted for full output (not needed with --stats-only). Required.'
        ),
    )

    scaffold_group = parser.add_mutually_exclusive_group(required=True)
    scaffold_group.add_argument(
        '-s', '--target-scaffolds', metavar='FILE',
        help=(
            'Text file of target scaffold names, one per line '
            '(e.g. mitochondrial chromosome names). Only the first column is used and a '
            'leading ">" is dropped, so a .fai file works too. '
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
            'Nicknames must be unique. Relative paths are read from the current directory. '
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
        '-q', '--mapq-threshold', type=_non_negative_int, default=0, metavar='INT',
        help=(
            'MAPQ threshold q. Reads with MAPQ > q whose primary alignment is on a target '
            'scaffold are written to the unique BAM. Reads with MAPQ <= q that touch a '
            'target scaffold go through the competitive NM comparison. '
            'Default: 0 (MAPQ >= 1 → unique, MAPQ 0 → compared). Common choices: 0, 1, 20, 25.'
        ),
    )
    parser.add_argument(
        '--drop-ties', action='store_true',
        help=(
            'Do not export ties (best target NM equal to best non-target NM). '
            'Ties are still counted in the summary. Default: ties are kept.'
        ),
    )
    parser.add_argument(
        '-t', '--threads', type=_positive_int, default=4, metavar='INT',
        help=(
            'Total threads for BAM (de)compression, split between the input BAM and '
            'the output BAMs (at least one per output BAM) (default: 4).'
        ),
    )
    parser.add_argument(
        '--samtools', default='samtools', metavar='PATH',
        help='Path to samtools, used only for indexing output BAMs (default: samtools).',
    )

    args = parser.parse_args()

    input_bam_path = os.path.abspath(args.input)
    output_directory = os.path.abspath(args.output_dir)
    mapq_threshold = args.mapq_threshold
    command_line = shlex.join(sys.argv).replace('\t', ' ').replace('\n', ' ')

    # --- Validate BAM ---
    if not os.path.exists(input_bam_path):
        print(f"ERROR: Input BAM not found: {input_bam_path}", file=sys.stderr)
        sys.exit(1)
    if not input_bam_path.endswith('.bam'):
        print("ERROR: Input must be a .bam file", file=sys.stderr)
        sys.exit(1)
    with pysam.AlignmentFile(input_bam_path, 'rb') as input_bam_file:
        header_scaffolds = set(input_bam_file.references)
        sort_order = input_bam_file.header.to_dict().get('HD', {}).get('SO')
    if not args.stats_only and sort_order in ('queryname', 'unsorted'):
        print(
            f"ERROR: The input BAM header says SO:{sort_order}. Full output needs a "
            f"coordinate-sorted BAM so the unique BAM can be indexed. "
            f"Run: {args.samtools} sort, or use --stats-only.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --- Load scaffold sets ---
    if args.target_scaffolds:
        refs = load_single_ref(os.path.abspath(args.target_scaffolds), args.sample_name)
    else:
        refs = load_refs_file(os.path.abspath(args.references))
    check_reference_names(refs)
    check_scaffolds_in_header(refs, header_scaffolds)

    mode = 'stats-only' if args.stats_only else 'full output'
    print(f"COMFRT {__version__}")
    print(f"Input BAM:      {input_bam_path}")
    print(f"References:     {len(refs)} ({',  '.join(reference_entry['name'] for reference_entry in refs)})")
    print(f"MAPQ threshold: MAPQ > {mapq_threshold} → unique; MAPQ <= {mapq_threshold} → compared")
    print(f"Ties:           {'dropped' if args.drop_ties else 'kept'}")
    print(f"Mode:           {mode}")

    run_pipeline(
        input_bam_path,
        refs,
        output_directory,
        args.stats_only,
        mapq_threshold,
        args.threads,
        args.samtools,
        drop_ties=args.drop_ties,
        command_line=command_line,
    )


if __name__ == '__main__':
    main()
