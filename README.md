# COmpetitively Mapped Fragment Recovery Tool (COMFRT)
## Competitive Mapping Filter

Recovers target-scaffold reads from a BAM aligned to a combined reference, using a single-pass competitive NM comparison.

Originally written for NUMT (Nuclear Mitochondrial DNA) recovery from ancient DNA
sequencing data, by Kirstine Tersbøl Melsen.

---

## How it works

Reads are classified into two pools based on MAPQ:

| Pool | Condition | Outcome |
|------|-----------|---------|
| **Unique** | MAPQ > threshold AND primary alignment on a target scaffold | Written to a cleaned BAM (non-target `@SQ` lines stripped) |
| **Ambiguous** | MAPQ = 0 AND touches a target scaffold (primary or BWA `XA` tag) | Classified by competitive NM comparison (see below) |
| Ignored | MAPQ = 0 AND does not touch any target scaffold | Dropped |

**Competitive NM comparison** — for each ambiguous read, the best edit distance (NM) on
target scaffolds is compared against the best NM on non-target scaffolds:

```
target NM  <  non-target NM  →  recovered  (target alignment is better)
target NM  == non-target NM  →  tie        (kept — target is equally good)
target NM  >  non-target NM  →  discarded  (non-target alignment is better)
```


> **IMPORTANT**: The input BAM must be aligned with BWA to a **combined reference**
> containing both target and non-target scaffolds, and must **retain MAPQ=0 reads**.
> Do not pre-filter with `samtools view -q` before running this script.
> The BAM must be coordinate-sorted and indexed.

---

## Requirements

- Python 3.8+
- [pysam](https://pysam.readthedocs.io/) >= 0.22
- samtools (only used for indexing the output unique BAM)

```bash
pip install pysam
```

---

## Usage

```
comfrt.py -i BAM (-s FILE | -r FILE) -o DIR [options]
```

### Single-reference mode

One scaffold set, one output name prefix:

```bash
comfrt.py \
    -i aligned.bam \
    -s scaffolds.txt \
    -o output/ \
    -n OUTPUT_PREFIX \
    -q 1 \
    -t 16
```

### Multi-reference mode

Multiple scaffold sets in a single BAM pass. The BAM is read **once** for all references:

```bash
comfrt.py \
    -i aligned.bam \
    -r refs.tsv \
    -o output/ \
    -t 16
```

### Stats-only mode

Skip all BAM/FASTQ output — only write per-reference `_summary.txt` files.
Useful before committing to a full run:

```bash
comfrt.py \
    -i aligned.bam \
    -r refs.tsv \
    -o output/ \
    --stats-only
```

---

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `-i / --input` | yes | Input BAM. Combined reference, coordinate-sorted, indexed, MAPQ=0 retained. |
| `-s / --target-scaffolds` | one of `-s`/`-r` | Plain-text file of target scaffold names, one per line. Use `-n` to set output prefix. |
| `-r / --references` | one of `-s`/`-r` | Two-column TSV: `scaffold_file_path  nickname`. One row per scaffold set. `#` lines ignored. |
| `-o / --output-dir` | yes | Output directory. Created if absent. With `-r` (full output), one subfolder per reference. |
| `-n / --sample-name` | no | Output filename prefix when using `-s` (default: `sample`). Ignored with `-r`. |
| `--stats-only` | no | Write only `_summary.txt` files; skip BAM/FASTQ output entirely. |
| `-q / --mapq-threshold` | no | Reads with MAPQ **strictly greater than** this go to the unique BAM. Default: `0` (MAPQ ≥ 1 → unique). |
| `-t / --threads` | no | Threads for pysam and samtools (default: `4`). |
| `--samtools` | no | Path to samtools executable (default: `samtools`). |

### `-r` / `--references` file format

Tab or space separated, `#` comments allowed:

```
# scaffold_file               nickname
/data/refs/mt_scaffolds.txt   mitochondria
/data/refs/y_scaffolds.txt    y_chromosome
```

Scaffold files are plain text, one scaffold name per line, matching sequence names in the
BAM header exactly.

The textfile for the scaffold list can be made by using faidx, taking the first column e.g: `cut -f1 > SCAFFOLDS.txt `

---

## Output files

### Full output (default)

With `-s` (single ref), files go directly in `output/`.
With `-r` (multi-ref), files go in `output/<nickname>/`.

| File | Description |
|------|-------------|
| `<name>_unique.bam` | Unique target reads (MAPQ > threshold). BAM header stripped of non-target `@SQ` lines. Indexed automatically. |
| `<name>_unique.bam.bai` | BAI index for unique BAM. |
| `<name>_recovered_R1.fq.gz` | Recovered ambiguous reads — R1 of a pair. |
| `<name>_recovered_R2.fq.gz` | Recovered ambiguous reads — R2 of a pair. |
| `<name>_recovered_merged.fq.gz` | Recovered ambiguous reads — merged or single-end. |
| `<name>_summary.txt` | Run statistics (see below). |

### Stats-only output

With `--stats-only`, only `<name>_summary.txt` files are written, directly in `output/`
regardless of whether `-s` or `-r` is used.

### Summary file contents

```
Competitive Mapping Filter - Summary
========================================
Reference:           Alien
Input BAM:           /path/to/aligned.bam
MAPQ threshold:      20
Stats only:          no

Total mapped reads:          4103781
Touch target scaffold:       3347917  (81.6%)

Unique target reads (MAPQ > 20): 2981495

Total ambiguous reads (MAPQ=0): 462012
  Target-touching ambiguous:     366422
  Not touching target:           95590

Target-touching ambiguous breakdown:
  Target-only (multiple alignments): 321562
  Recovered (target strictly better NM): 1193
  Ties (equal NM, kept):                 42471
  Discarded (non-target better):         1196

Recovered FASTQ:
  R1 reads:        31555
  R2 reads:        31555
  Merged/SE reads: 302116
  Total:           365226
```

---

## Workflow

### Overview

```
reads.fq
    │
    ├─── Step 1 ──► Map to COMBINED reference       (target + non-target)
    │                MAPQ=0 must be retained                │
    │                                                       ▼
    │                           Step 2 ──► comfrt.py
    │                                   │
    │                      ┌────────────┴────────────┐
    │                      ▼                         ▼
    │              <name>_unique.bam        recovered_R1/R2/merged.fq.gz
    │              (MAPQ > threshold,       (ambiguous reads where target
    │               header = target only)    alignment wins or ties)
    │                      │                       │
    ├─── Step 3  ──► Map to TARGET reference only  │          Step 4
    │                (used for final BAM)          └──► Remap recovered reads
    │                       │                           to TARGET reference
    │                       │                                    │
    │                       └──────────────┬─────────────────────┘
    │                                      ▼
    └─────────────────────────────► Step 5: Merge BAMs
                                    samtools merge final.bam \
                                        target_only.bam \
                                        remapped_recovered.bam
```

---

## Notes

- **MAPQ threshold** (`-q`): default `0` means any read with MAPQ ≥ 1 on a target scaffold
  is "unique". For stricter uniqueness, use `-q 20` or `-q 25`. Common BWA behaviour:
  MAPQ=0 means the read maps equally well to ≥ 2 places; MAPQ=25/37 indicates a unique hit.

- **Multi-reference single pass**: when using `-r`, all scaffold sets are evaluated
  simultaneously in one pass. A read can be classified differently for each reference —
  e.g. "unique" for one target but "not touching" for another. Each reference
  is fully independent.

- **Ties are kept**: reads where the best target NM equals the best non-target NM are
  written to the recovered FASTQ (counted separately as "ties" in the summary).

- **PE awareness**: recovered reads are split into R1/R2/merged based on BAM pair flags.
  All three FASTQ files are always written, even if empty.

- **Header compatibility for merging (Steps 3–5)**: `samtools merge` requires all input
  BAMs to share the same `@SQ` header lines. The unique BAM from Step 2 is already stripped
  to target-scaffold headers only. Steps 3 and 4 must therefore align against the same
  isolated target reference — not the combined one from Step 1 — to produce a compatible header.
  Attempting to merge the combined-reference BAM directly would fail.
