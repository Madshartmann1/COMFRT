# COmpetitively Mapped Fragment Recovery Tool (COMFRT)
## Competitive Mapping Filter

Recovers target-scaffold reads from a BAM aligned to a combined reference, using a single-pass competitive NM comparison.

Originally written for NUMT (Nuclear Mitochondrial DNA) recovery from ancient DNA
sequencing data, by Kirstine Tersbøl Melsen.

Current version: **1.0.0** (`comfrt.py --version`). The version is written into every
summary file and into the `@PG` line of the unique BAM. See [Version history](#version-history).

---

## How it works

Only primary alignments are used; secondary and supplementary records are skipped.
With a MAPQ threshold `q` (`-q`, default `0`), every mapped read falls into exactly one of these rows:

| Pool | Condition | Outcome |
|------|-----------|---------|
| **Unique** | MAPQ > q AND primary alignment on a target scaffold | Written to the unique BAM (non-target `@SQ` lines stripped) |
| Non-target | MAPQ > q AND primary alignment on a non-target scaffold | Ignored |
| **Ambiguous** | MAPQ ≤ q AND touches a target scaffold (primary alignment or BWA `XA` tag) | Classified by competitive NM comparison (see below) |
| Not touching | MAPQ ≤ q AND does not touch any target scaffold | Ignored |

With the default `-q 0`, MAPQ ≥ 1 on a target scaffold is unique and only MAPQ 0 reads are compared.

**Competitive NM comparison**: for each ambiguous read, the best edit distance (NM) on
target scaffolds is compared against the best NM on non-target scaffolds, using the
primary alignment and every hit listed in the `XA` tag:

```
target NM  <  non-target NM   →  recovered    (target alignment is better; exported)
target NM  == non-target NM   →  tie          (exported, unless --drop-ties)
target NM  >  non-target NM   →  discarded    (non-target alignment is better)
no non-target hit listed      →  target-only  (exported; nothing to compare against)
```

Exported reads are written as FASTQ, to be remapped to the target reference alone
(see [Workflow](#workflow)).

> **IMPORTANT**: The input BAM must be aligned with BWA to a **combined reference**
> containing both target and non-target scaffolds, and must **retain MAPQ=0 reads**.
> Do not pre-filter with `samtools view -q` before running this script.
> For full output the BAM must be coordinate-sorted (the unique BAM is indexed). An index is not needed.

### The XA tag limit: reads that are never compared

BWA only writes the `XA` tag when a read has few enough alternative hits:

| Aligner | Option | Default |
|---------|--------|---------|
| `bwa samse` | `-n` max hits listed in `XA` | 3 |
| `bwa sampe` | `-n` (proper pairs) / `-N` (discordant pairs) | 3 / 10 |
| `bwa mem` | `-h` | 5 |

A read with more hits than this gets **no `XA` tag at all**, so COMFRT only sees its primary
alignment. BWA picks that primary alignment at random among equally good hits. If it lands on
a target scaffold, the read is kept as "target-only" without any comparison, even if some of
its unlisted hits are on non-target scaffolds. If it lands on a non-target scaffold, the read
counts as "not touching".

The summary reports these reads on their own line: `Target-only, no XA tag (not compared)`.
In repetitive targets this can be the largest group of exported reads. To compare more of
them, raise `-n` (and `-N` for `sampe`) in the Step 1 mapping so that more alternative hits
are listed.

### bwa aln vs bwa mem

COMFRT was developed and tested with `bwa aln` + `samse`/`sampe`. `bwa mem` output also works
(secondary and supplementary records are skipped). But NM does not count soft-clipped bases,
so with `bwa mem` a heavily clipped hit can look better than a full-length one.

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

Skip all BAM/FASTQ output and write only per-reference `_summary.txt` files.
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
| `-i / --input` | yes | Input BAM. Combined reference, MAPQ=0 retained, coordinate-sorted for full output. |
| `-s / --target-scaffolds` | one of `-s`/`-r` | Text file of target scaffold names, one per line. Use `-n` to set output prefix. |
| `-r / --references` | one of `-s`/`-r` | Two-column TSV: `scaffold_file_path  nickname`. One row per scaffold set. `#` lines ignored. |
| `-o / --output-dir` | yes | Output directory. Created if absent. With `-r` (full output), one subfolder per reference. |
| `-n / --sample-name` | no | Output filename prefix when using `-s` (default: `sample`). Ignored with `-r`. |
| `--stats-only` | no | Write only `_summary.txt` files; skip BAM/FASTQ output entirely. |
| `-q / --mapq-threshold` | no | Threshold `q`: MAPQ **> q** on a target scaffold → unique BAM; MAPQ **≤ q** touching a target → compared. Default: `0`. |
| `--drop-ties` | no | Do not export ties. They are still counted in the summary. Default: ties are kept. |
| `-t / --threads` | no | Total threads for BAM (de)compression, split between the input and output BAMs (default: `4`). |
| `--samtools` | no | Path to samtools executable (default: `samtools`). |
| `--version` | no | Print the version and exit. |

### Scaffold files

One scaffold name per line, matching the `@SQ SN:` names in the BAM header exactly.
Only the first whitespace-separated column is read and a leading `>` is dropped, so a
`samtools faidx` index works directly, or can be cut down to names:

```bash
samtools faidx target.fa
cut -f1 target.fa.fai > SCAFFOLDS.txt
```

COMFRT stops with an error if none of the names are in the BAM header, and warns if only some are.

### `-r` / `--references` file format

Tab or space separated, `#` comments allowed. Nicknames must be unique. Relative paths are
read from the current directory, and paths cannot contain spaces:

```
# scaffold_file               nickname
/data/refs/mt_scaffolds.txt   mitochondria
/data/refs/y_scaffolds.txt    y_chromosome
```

---

## Output files

### Full output (default)

With `-s` (single ref), files go directly in `output/`.
With `-r` (multi-ref), files go in `output/<nickname>/`.

| File | Description |
|------|-------------|
| `<name>_unique.bam` | Unique target reads (MAPQ > q). Header has only the target `@SQ` lines, plus a `@PG` line for COMFRT. Indexed automatically. |
| `<name>_unique.bam.bai` | BAI index for unique BAM. |
| `<name>_recovered_R1.fq.gz` | Exported ambiguous reads: R1 of a pair where both mates were exported. |
| `<name>_recovered_R2.fq.gz` | Exported ambiguous reads: R2 of a pair where both mates were exported. |
| `<name>_recovered_merged.fq.gz` | Exported ambiguous reads: merged, single-end, or one mate of a pair. |
| `<name>_summary.txt` | Run statistics (see below). |

In the unique BAM, a read whose mate maps to a non-target scaffold has its mate recorded
as unmapped (that scaffold is not in the header).

### Stats-only output

With `--stats-only`, only `<name>_summary.txt` files are written, directly in `output/`
regardless of whether `-s` or `-r` is used.

### Summary file contents

Every mapped primary read is counted in exactly one line of the first block.
Example (`-q 20`):

```
COMFRT (Competitive Mapping Filter) - Summary
================================================
COMFRT version:                                 1.0.0
Reference:                                      Alien
Input BAM:                                      /path/to/aligned.bam
MAPQ threshold:                                 20  (unique: MAPQ > 20; ambiguous: MAPQ <= 20)
Ties:                                           kept
Stats only:                                     no

Total mapped reads (primary alignments):        4103781
  Unique target reads (MAPQ > 20):              2981495
  Non-target reads (MAPQ > 20):                 525275
  Ambiguous reads (MAPQ <= 20):                 597011
    Not touching target:                        115861
    Target-touching:                            481150
      of which MAPQ 1-20:                       114728

Touch target (unique + target-touching):        3462645  (84.4%)

Target-touching ambiguous breakdown:
  Recovered (target strictly better NM):        17634
  Ties (equal NM, kept):                        42491
  Target-only, XA lists only target hits:       70821
  Target-only, no XA tag (not compared):        345288
  Discarded (non-target strictly better NM):    4916

Recovered FASTQ (for remapping to the target):
  R1 reads:                                     41077
  R2 reads:                                     41077
  Merged/SE reads:                              394080
  Total:                                        476234
```

The three first-block counts add up to the total, and the five breakdown lines add up to
"Target-touching". Everything in the breakdown except "Discarded" is exported. Lines for
skipped records (no NM tag, or secondary/supplementary) only appear when they are non-zero.

---

## Workflow

### Overview

```
reads.fq
   │
   ▼
Step 1  Map to the COMBINED reference (target + non-target)
        keep MAPQ 0, coordinate-sort
   │
   ▼
Step 2  comfrt.py
   │
   ├────────────────────────────► <name>_unique.bam
   │                              (MAPQ > q on target, target-only header)
   ▼                                          │
<name>_recovered_{R1,R2,merged}.fq.gz         │
(ambiguous reads where the target             │
 wins, ties or has no competitor)             │
   │                                          │
   ▼                                          │
Step 3  Remap to the TARGET reference only    │
   │                                          │
   └──────────────────┬───────────────────────┘
                      ▼
Step 4  samtools merge final.bam <name>_unique.bam remapped.bam
                      │
                      ▼
Step 5  Mark duplicates again on the merged BAM
```

Only the exported reads are remapped. Do **not** merge a target-only mapping of *all* reads:
that would bring back exactly the reads the competitive mapping removed, and count the
exported reads twice.

### Example commands

Single-end/merged reads with ancient-DNA `bwa aln` settings; adapt to your data.
`combined.fa` holds target + non-target scaffolds; `target.fa` holds the target scaffolds
from the same FASTA (same names and lengths).

```bash
# Step 1: competitive mapping; raise samse -n so more alternative hits are listed in XA
bwa aln -l 999 -n 0.01 -t 8 combined.fa reads.fq.gz > reads.sai
bwa samse -n 50 -r '@RG\tID:sample\tSM:sample\tLB:sample' combined.fa reads.sai reads.fq.gz \
    | samtools sort -o competitive.bam -

# Step 2
comfrt.py -i competitive.bam -s target_scaffolds.txt -o comfrt_out -n sample -q 20

# Step 3: remap the exported reads to the target only, with the same settings.
# Same SM/LB as Step 1 so duplicate marking treats both pools as one library.
bwa aln -l 999 -n 0.01 -t 8 target.fa comfrt_out/sample_recovered_merged.fq.gz > remap.sai
bwa samse -r '@RG\tID:sample.remap\tSM:sample\tLB:sample' target.fa remap.sai \
        comfrt_out/sample_recovered_merged.fq.gz \
    | samtools view -u -F 4 - | samtools sort -o remap.bam -
# Paired reads: bwa aln on R1 and R2, then
#   bwa sampe ... | samtools fixmate -m - - | samtools sort -o remap_pe.bam -

# Step 4: merge with the unique BAM
samtools merge final.bam comfrt_out/sample_unique.bam remap.bam

# Step 5: mark duplicates again
samtools markdup final.bam final.markdup.bam
```

Do not filter the remapped reads on MAPQ: reads that still map equally well to several
places *within* the target are expected and kept by design.

---

## Notes

- **MAPQ threshold** (`-q`): default `0` means any read with MAPQ ≥ 1 on a target scaffold
  is "unique". For stricter uniqueness, use `-q 20` or `-q 25`. Reads with MAPQ ≤ q are then
  compared instead of being trusted. Common BWA behaviour: MAPQ=0 means the read maps
  equally well to ≥ 2 places; MAPQ=25/37 indicates a unique hit.

- **Multi-reference single pass**: when using `-r`, all scaffold sets are evaluated
  simultaneously in one pass. A read can be classified differently for each reference,
  e.g. "unique" for one target but "not touching" for another. Each reference
  is fully independent.

- **Ties are kept by default**: reads where the best target NM equals the best non-target NM
  are exported (counted separately as "ties" in the summary). In a contaminated sample a tie
  is a read that fits a non-target genome exactly as well as the target; use `--drop-ties`
  to leave them out.

- **Paired-end reads**: each read is classified on its own, not as a pair. Exported reads are
  split into R1/R2 when both mates were exported, otherwise they go to the merged/SE file.
  So the unique BAM can hold reads whose mate was exported or discarded; their flags still
  describe the mate as mapped. All three FASTQ files are always written, even if empty.

- **Mark duplicates again after the merge (Step 5)**: PCR duplicates of MAPQ 0 reads can be
  placed on *different* equally good hits in the competitive mapping, so they escape a first
  duplicate-marking pass. After the remap to the target alone they land on the same position.

- **Header compatibility for merging (Steps 3–4)**: the unique BAM's header has only the
  target `@SQ` lines. Remap against a target-only reference built from the same FASTA, so
  names and lengths match and the merged BAM contains only the target scaffolds.

- **Threads**: `-t` is a total. Half goes to reading the input BAM and half is shared by
  the output BAMs (at least one thread each). The Python classification loop itself is
  single-threaded.

---

## Version history

**1.0.0**: first versioned release. Compared with earlier, unversioned copies of `comfrt.py`:

- Reads with 0 < MAPQ ≤ q are now compared like MAPQ 0 reads. Before, such reads on a target
  scaffold were dropped without being counted. Results with `-q 0` are unchanged.
- The summary accounts for every mapped primary read. It splits "target-only" into reads
  whose `XA` lists only target hits and reads with no `XA` tag, records the version, and
  counts a read as "recovered" when its primary alignment is the losing non-target hit
  (before, only `XA` was checked).
- New `--drop-ties` and `--version` options.
- Secondary and supplementary records are skipped. With `bwa mem`, a supplementary record
  was exported as an extra, shorter copy of its read and counted twice.
- The unique BAM is written directly with the target-only header. This fixes a hang when
  many unique reads had their mate on a non-target scaffold, and records such mates as
  unmapped instead of leaving inconsistent mate fields. A `@PG` line is added.
- Scaffold names are checked against the BAM header; `.fai` files and `>` header lines are
  accepted; duplicate nicknames are rejected; a BAM index is no longer required.
- `-t` is now a total thread count.
