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


def _header_line(tag, fields):
    """Format a SAM header line from a tag (HD/SQ/RG/PG) and key/value fields."""
    return "@" + tag + "\t" + "\t".join(
        f"{field_key}:{field_value}" for field_key, field_value in fields.items()
    )


def _write_target_only_header(header_path, bam_path, target_scaffolds):
    """
    Build a SAM header file that keeps only target scaffold @SQ entries.

    Preserves @HD, @RG, @PG, and @CO records from the input BAM header so
    provenance and read-group metadata remain intact.
    """
    with pysam.AlignmentFile(bam_path, 'rb') as input_bam_file:
        input_header = input_bam_file.header.to_dict()

    with open(header_path, 'w') as output_header_file:
        if 'HD' in input_header:
            output_header_file.write(_header_line('HD', input_header['HD']) + '\n')

        if 'SQ' in input_header:
            for sequence_entry in input_header['SQ']:
                sequence_name = sequence_entry.get('SN')
                if sequence_name in target_scaffolds:
                    output_header_file.write(_header_line('SQ', sequence_entry) + '\n')

        for read_group_entry in input_header.get('RG', []):
            output_header_file.write(_header_line('RG', read_group_entry) + '\n')
        for program_entry in input_header.get('PG', []):
            output_header_file.write(_header_line('PG', program_entry) + '\n')
        for comment_entry in input_header.get('CO', []):
            output_header_file.write(f"@CO\t{comment_entry}\n")


def _rewrite_bam_with_header(unique_bam, header_path, samtools_path):
    """
    Re-encode an existing BAM with a replacement header.

    Uses streaming `samtools view` to avoid loading the full BAM into memory,
    then atomically replaces the original BAM on success.
    """
    temporary_bam = unique_bam + '.tmp.bam'

    with open(header_path, 'rb') as header_file_handle:
        replacement_header_bytes = header_file_handle.read()

    samtools_stream_reader = subprocess.Popen(
        [samtools_path, 'view', unique_bam],
        stdout=subprocess.PIPE,
    )
    samtools_stream_writer = subprocess.Popen(
        [samtools_path, 'view', '-b', '-o', temporary_bam, '-'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    samtools_stream_writer.stdin.write(replacement_header_bytes)
    while True:
        sam_chunk = samtools_stream_reader.stdout.read(1024 * 1024)
        if not sam_chunk:
            break
        samtools_stream_writer.stdin.write(sam_chunk)
    samtools_stream_writer.stdin.close()

    reader_return_code = samtools_stream_reader.wait()
    writer_return_code = samtools_stream_writer.wait()
    writer_stderr_bytes = (
        samtools_stream_writer.stderr.read()
        if samtools_stream_writer.stderr is not None
        else b''
    )

    if reader_return_code != 0 or writer_return_code != 0:
        if os.path.exists(temporary_bam):
            os.remove(temporary_bam)
        writer_stderr_message = writer_stderr_bytes.decode(errors='replace').strip()
        raise RuntimeError(
            f"Failed to rebuild BAM with target-only header ({unique_bam}). "
            f"samtools view return codes: {reader_return_code}, {writer_return_code}. "
            f"{writer_stderr_message}"
        )

    os.replace(temporary_bam, unique_bam)
    os.remove(header_path)


# ============================================================================
# REFERENCE LOADING
# ============================================================================

def _load_scaffold_set(path, label='Scaffolds file'):
    if not os.path.exists(path):
        print(f"ERROR: {label} not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as input_file_handle:
        scaffolds = set(line.strip() for line in input_file_handle if line.strip())
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


# ============================================================================
# PIPELINE
# ============================================================================

def run_pipeline(bam, refs, outdir, stats_only, mapq, threads, samtools_path):
    """
    Single pass over `bam`, classifying reads against all scaffold sets simultaneously.

    refs: list of {'name': str, 'scaffolds': set}
    """
    reference_count = len(refs)
    os.makedirs(outdir, exist_ok=True)

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
            'unique':      0,
            'not_target':  0,
            'target_only': 0,
            'recovered':   0,
            'ties':        0,
            'discarded':   0,
            'skipped_nm':  0,
        }
        reference_config['r1_count'] = 0
        reference_config['r2_count'] = 0
        reference_config['merged_count'] = 0
        reference_config['pending_pairs_by_query_name'] = {}

    total_mapped_reads = 0
    print(f"\n[1/2] Scanning BAM ({reference_count} reference{'s' if reference_count > 1 else ''})...")

    with contextlib.ExitStack() as stack:
        input_bam = stack.enter_context(pysam.AlignmentFile(bam, 'rb', threads=threads))

        if not stats_only:
            for reference_config in refs:
                reference_config['bam_out'] = stack.enter_context(
                    pysam.AlignmentFile(reference_config['unique_bam'], 'wb', template=input_bam, threads=threads)
                )
                reference_config['f1'] = stack.enter_context(gzip.open(reference_config['r1_out'], 'wt'))
                reference_config['f2'] = stack.enter_context(gzip.open(reference_config['r2_out'], 'wt'))
                reference_config['fm'] = stack.enter_context(gzip.open(reference_config['merged_out'], 'wt'))

        for read in input_bam:
            if read.is_unmapped:
                continue
            total_mapped_reads += 1

            for reference_config in refs:
                target_scaffolds = reference_config['scaffolds']
                reference_stats = reference_config['stats']

                # Unique: MAPQ > threshold AND primary alignment on a target scaffold
                if read.mapping_quality > mapq and read.reference_name in target_scaffolds:
                    if not stats_only:
                        reference_config['bam_out'].write(read)
                    reference_stats['unique'] += 1
                    continue  # skip ambiguous classification for this ref

                # Ambiguous: MAPQ = 0
                if read.mapping_quality == 0:
                    if not touches_target(read, target_scaffolds):
                        reference_stats['not_target'] += 1
                        continue

                    result = classify(read, target_scaffolds)

                    if result is None:
                        reference_stats['skipped_nm'] += 1
                        continue

                    if result in ('keep', 'tie'):
                        if not stats_only:
                            recovered_fastq_record = read_to_fastq_str(read)
                            if recovered_fastq_record is not None:
                                write_recovered_read_to_bucket(read, recovered_fastq_record, reference_config)

                        if result == 'tie':
                            reference_stats['ties'] += 1
                        else:
                            xa_value = read.get_tag('XA') if read.has_tag('XA') else None
                            secondary_alignments = parse_xa(xa_value) if xa_value else []
                            if any(alignment['scaffold'] not in target_scaffolds for alignment in secondary_alignments):
                                reference_stats['recovered'] += 1
                            else:
                                reference_stats['target_only'] += 1
                    else:
                        reference_stats['discarded'] += 1

        if not stats_only:
            for reference_config in refs:
                flush_unpaired_recovered_reads(reference_config)

    # Rebuild and index unique BAMs with target-only headers
    if not stats_only:
        print(f"[2/2] Rebuilding headers and indexing {reference_count} unique BAM{'s' if reference_count > 1 else ''}...")
        for reference_config in refs:
            target_header_path = reference_config['unique_bam'] + '.target_header.sam'
            _write_target_only_header(target_header_path, bam, reference_config['scaffolds'])
            _rewrite_bam_with_header(reference_config['unique_bam'], target_header_path, samtools_path)
            subprocess.run([samtools_path, 'index', reference_config['unique_bam']], check=True)

    # Print and write summaries
    for reference_config in refs:
        reference_stats = reference_config['stats']
        recovered_r1_count = reference_config['r1_count']
        recovered_r2_count = reference_config['r2_count']
        recovered_merged_count = reference_config['merged_count']
        ambiguous_target_touching_reads = (
            reference_stats['target_only']
            + reference_stats['recovered']
            + reference_stats['ties']
            + reference_stats['discarded']
            + reference_stats['skipped_nm']
        )
        total_ambiguous_reads = ambiguous_target_touching_reads + reference_stats['not_target']
        target_touching_reads = (
            reference_stats['unique']
            + reference_stats['target_only']
            + reference_stats['recovered']
            + reference_stats['ties']
            + reference_stats['discarded']
            + reference_stats['skipped_nm']
        )
        target_touching_percent = 100.0 * target_touching_reads / total_mapped_reads if total_mapped_reads else 0.0

        print(f"\n{'=' * 45}")
        print(f"Reference: {reference_config['name']}")
        print(f"  Total mapped reads:            {total_mapped_reads}")
        print(f"  Touch target scaffold:         {target_touching_reads}  ({target_touching_percent:.1f}%)")
        print(f"  Total ambiguous (MAPQ=0):      {total_ambiguous_reads}")
        print(f"    └─ target-touching:          {ambiguous_target_touching_reads}")
        print(f"    └─ not touching target:      {reference_stats['not_target']}")
        print(f"  Unique target reads:           {reference_stats['unique']}")
        print(f"  Ambiguous → target-only:       {reference_stats['target_only']}")
        print(f"  Ambiguous → recovered:         {reference_stats['recovered']}")
        print(f"  Ambiguous → ties (kept):       {reference_stats['ties']}")
        print(f"  Ambiguous → discarded:         {reference_stats['discarded']}")
        if reference_stats['skipped_nm']:
            print(f"  Skipped (no NM tag):           {reference_stats['skipped_nm']}")
        if not stats_only:
            print(
                f"  Recovered FASTQ: R1={recovered_r1_count}  "
                f"R2={recovered_r2_count}  "
                f"Merged={recovered_merged_count}  "
                f"Total={recovered_r1_count+recovered_r2_count+recovered_merged_count}"
            )

        with open(reference_config['summary_file'], 'w') as summary_file_handle:
            summary_file_handle.write("Competitive Mapping Filter - Summary\n")
            summary_file_handle.write("=" * 40 + "\n")
            summary_file_handle.write(f"Reference:           {reference_config['name']}\n")
            summary_file_handle.write(f"Input BAM:           {bam}\n")
            summary_file_handle.write(f"MAPQ threshold:      {mapq}\n")
            summary_file_handle.write(f"Stats only:          {'yes' if stats_only else 'no'}\n\n")
            summary_file_handle.write(f"Total mapped reads:          {total_mapped_reads}\n")
            summary_file_handle.write(
                f"Touch target scaffold:       {target_touching_reads}  "
                f"({target_touching_percent:.1f}%)\n\n"
            )
            summary_file_handle.write(
                f"Unique target reads (MAPQ > {mapq}): {reference_stats['unique']}\n\n"
            )
            summary_file_handle.write(f"Total ambiguous reads (MAPQ=0): {total_ambiguous_reads}\n")
            summary_file_handle.write(
                f"  Target-touching ambiguous:     {ambiguous_target_touching_reads}\n"
            )
            summary_file_handle.write(
                f"  Not touching target:           {reference_stats['not_target']}\n\n"
            )
            summary_file_handle.write("Target-touching ambiguous breakdown:\n")
            summary_file_handle.write(
                f"  Target-only (multiple alignments): {reference_stats['target_only']}\n"
            )
            summary_file_handle.write(
                f"  Recovered (target strictly better NM): {reference_stats['recovered']}\n"
            )
            summary_file_handle.write(
                f"  Ties (equal NM, kept):                 {reference_stats['ties']}\n"
            )
            summary_file_handle.write(
                f"  Discarded (non-target better):         {reference_stats['discarded']}\n"
            )
            if reference_stats['skipped_nm']:
                summary_file_handle.write(
                    f"  Skipped (no NM tag):                   {reference_stats['skipped_nm']}\n"
                )
            if not stats_only:
                summary_file_handle.write("\nRecovered FASTQ:\n")
                summary_file_handle.write(f"  R1 reads:        {recovered_r1_count}\n")
                summary_file_handle.write(f"  R2 reads:        {recovered_r2_count}\n")
                summary_file_handle.write(f"  Merged/SE reads: {recovered_merged_count}\n")
                summary_file_handle.write(
                    f"  Total:           {recovered_r1_count+recovered_r2_count+recovered_merged_count}\n"
                )

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

    input_bam_path = os.path.abspath(args.input)
    output_directory = os.path.abspath(args.output_dir)
    mapq_threshold = args.mapq_threshold

    # --- Validate BAM ---
    if not os.path.exists(input_bam_path):
        print(f"ERROR: Input BAM not found: {input_bam_path}", file=sys.stderr)
        sys.exit(1)
    if not input_bam_path.endswith('.bam'):
        print("ERROR: Input must be a .bam file", file=sys.stderr)
        sys.exit(1)
    has_index = (
        os.path.exists(input_bam_path + '.bai')
        or os.path.exists(input_bam_path.replace('.bam', '.bai'))
    )
    if not has_index:
        print(f"ERROR: BAM index not found. Run: {args.samtools} index {input_bam_path}", file=sys.stderr)
        sys.exit(1)

    # --- Load scaffold sets ---
    if args.target_scaffolds:
            refs = load_single_ref(os.path.abspath(args.target_scaffolds), args.sample_name)
    else:
        refs = load_refs_file(os.path.abspath(args.references))

    mode = 'stats-only' if args.stats_only else 'full output'
    print(f"Input BAM:      {input_bam_path}")
    print(f"References:     {len(refs)} ({',  '.join(reference_entry['name'] for reference_entry in refs)})")
    print(f"MAPQ threshold: MAPQ > {mapq_threshold} → unique")
    print(f"Mode:           {mode}")

    run_pipeline(
        input_bam_path,
        refs,
        output_directory,
        args.stats_only,
        mapq_threshold,
        args.threads,
        args.samtools,
    )


if __name__ == '__main__':
    main()
