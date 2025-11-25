# NUMTS Recovery Pipeline

A Snakemake-based pipeline for identifying Nuclear Mitochondrial DNA segments (NUMTs) from sequencing data. This pipeline can process raw reads or pre-aligned BAM files to identify sequences that map to both mitochondrial and nuclear genomes.

## Overview

This pipeline identifies NUMTs by detecting reads that map to both the mitochondrial genome and nuclear genome. It uses a Python wrapper script to configure and execute a Snakemake pipeline with multiple processing steps.

## Features

- **Multiple input formats**: Accepts paired-end reads, single-end reads, or pre-aligned BAM files
- **Automatic input detection**: Intelligently detects input type and adapts pipeline accordingly
- **Quality control**: Integrated read preprocessing with fastp
- **Flexible alignment**: Configurable BWA alignment parameters
- **NUMTS identification**: Uses FragmentsMappingTwicePY.py to identify NUMTS sequences
- **Optional remapping**: Can generate remapped BAM files of recovered (non-NUMTS) sequences
- **Comprehensive logging**: Detailed logs for each processing step

## Requirements

### Software Dependencies

- Python 3.6+
- Snakemake 5.0+
- BWA
- samtools
- fastp (for read preprocessing)

### Python Packages

```bash
pip install pyyaml numpy pandas snakemake
```

## Installation

1. Clone or download this repository
2. Ensure all dependencies are installed and in your PATH
3. Make the wrapper script executable:

```bash
chmod +x numts_pipeline.py
```

## Usage

### Basic Usage

#### For paired-end reads:
```bash
./numts_pipeline.py \
    -i /path/to/reads/directory/ \
    -r /path/to/reference.fasta \
    -m NC_010642.1 \
    -o /path/to/output/
```

#### For pre-aligned BAM file:
```bash
./numts_pipeline.py \
    -i /path/to/aligned.bam \
    -r /path/to/reference.fasta \
    -m NC_010642.1 \
    -o /path/to/output/
```

### Advanced Usage with Custom Parameters

```bash
./numts_pipeline.py \
    -i /path/to/reads/ \
    -r /path/to/reference.fasta \
    -m NC_010642.1 \
    -o /path/to/output/ \
    -s tiger_sample \
    -t 16 \
    --min-quality 30 \
    --min-length 35 \
    --bwa-mismatch-penalty 0.03 \
    --remap \
    --remap-quality 25
```

### With Multiple Mitochondrial Scaffolds

```bash
./numts_pipeline.py \
    -i /path/to/reads/ \
    -r /path/to/reference.fasta \
    -m "NC_010642.1,MT_scaffold2" \
    -o /path/to/output/
```

## Command-Line Arguments

### Required Arguments

| Argument | Description |
|----------|-------------|
| `-i, --input` | Input: BAM file, directory with FASTQ reads, or comma-separated R1,R2 files |
| `-r, --reference` | Reference genome FASTA file (must contain mitochondrial scaffold) |
| `-m, --mitochondrial-scaffold` | Name(s) of mitochondrial scaffold(s) in reference (comma-separated) |
| `-o, --output-dir` | Output directory for results |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-s, --sample-name` | "sample" | Sample name for output files |
| `-t, --threads` | 4 | Total number of threads/cores to use |

### Preprocessing Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-quality` | 25 | Minimum quality score for fastp (-q parameter) |
| `--min-length` | 30 | Minimum read length for fastp (-l parameter) |

### BWA Alignment Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--bwa-seed-length` | 999 | BWA seed length (-l parameter) |
| `--bwa-max-gap-opens` | 2 | BWA maximum gap opens (-o parameter) |
| `--bwa-mismatch-penalty` | 0.04 | BWA mismatch penalty (-n parameter) |

### Pipeline Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--remap` | False | Generate remapped BAM file of recovered reads |
| `--remap-quality` | 20 | Minimum mapping quality for remapping step |
| `--keep-temp` | False | Keep temporary files |
| `--dry-run` | False | Perform Snakemake dry run only |

### Tool Paths

| Argument | Default | Description |
|----------|---------|-------------|
| `--fastp-path` | "fastp" | Path to fastp executable |
| `--bwa-path` | "bwa" | Path to BWA executable |
| `--samtools-path` | "samtools" | Path to samtools executable |

## Pipeline Workflow

### For FASTQ Input (Paired-End)

1. **Input Detection**: Identify R1 and R2 files
2. **Preprocessing**: fastp deduplication, quality filtering, and merging
3. **Reference Indexing**: Index reference genome for BWA (if needed)
4. **Alignment**: BWA aln with mapq=0 to retain multi-mapping reads
5. **BAM Cleanup**: Remove duplicates, sort, and index
6. **Statistics**: Calculate mapping statistics
7. **MT Filtering**: Extract reads mapping to mitochondrial scaffold
8. **NUMTS Identification**: Run FragmentsMappingTwicePY.py
9. **[Optional] Remapping**: Remap recovered sequences to MT genome

### For FASTQ Input (Single-End)

Same as paired-end, but skips merging step in preprocessing.

### For BAM Input

1. **Input Validation**: Check BAM file and warn about mapq=0 requirement
2. **BAM Cleanup**: Remove duplicates, sort, and index
3. **Statistics**: Calculate mapping statistics
4. **MT Filtering**: Extract reads mapping to mitochondrial scaffold
5. **NUMTS Identification**: Run FragmentsMappingTwicePY.py
6. **[Optional] Remapping**: Remap recovered sequences to MT genome

## Output Files

### Standard Output (Always Generated)

| File | Description |
|------|-------------|
| `numts.txt` | List of read IDs identified as NUMTs (map better to nuclear than MT) |
| `records_<MT_NAME>.txt` | Summary statistics (MT-only mappings, MT+nuclear mappings, NUMTS count) |
| `config.yaml` | Configuration file used for the run |
| `logs/` | Directory containing log files for each step |
| `stats/<sample>_mapping_stats.txt` | Mapping statistics |

### Optional Output (When --remap is specified)

| File | Description |
|------|-------------|
| `<sample>_recovered.bam` | BAM file of recovered (non-NUMTS) sequences remapped to MT genome |
| `<sample>_recovered.bam.bai` | Index for recovered BAM file |

### Intermediate Files (Kept with --keep-temp)

| Directory | Contents |
|-----------|----------|
| `preprocessed/` | Merged/filtered FASTQ files from fastp |
| `alignment/` | Initial alignment BAM files |
| `processed/` | Cleaned and sorted BAM files |
| `filtered/` | MT-filtered BAM and text files |
| `remapping/` | Files for remapping step |

## Important Notes

### BAM File Requirements

**CRITICAL**: If providing a pre-aligned BAM file, it MUST have been aligned with **minimum mapping quality = 0** to retain multi-mapping reads. This is essential for NUMTS detection. The pipeline will warn you about this requirement.

### Mitochondrial Genome

The mitochondrial scaffold(s) specified with `-m` must be present in the reference genome provided with `-r`. The pipeline will extract the MT genome for remapping if `--remap` is specified.

### Single-End Reads

If only single-end reads are found, the pipeline will automatically switch to single-end mode and issue a warning. NUMTS detection may have reduced sensitivity in single-end mode.

### Memory and Performance

- Processing time depends on input size and number of threads
- For large genomes, ensure sufficient RAM (recommend 8GB+ per sample)
- Use more threads (`-t`) to speed up processing

## Examples

### Example 1: Basic NUMTS Identification from Reads

```bash
./numts_pipeline.py \
    -i /data/tiger/reads/ \
    -r /ref/tiger_genome.fasta \
    -m NC_010642.1 \
    -o /results/tiger_numts/ \
    -s tiger_001 \
    -t 8
```

**Output**: `numts.txt` and `records_NC_010642.1.txt`

### Example 2: Full Pipeline with Remapping

```bash
./numts_pipeline.py \
    -i /data/tiger/reads/ \
    -r /ref/tiger_genome.fasta \
    -m NC_010642.1 \
    -o /results/tiger_numts/ \
    -s tiger_001 \
    -t 16 \
    --remap \
    --remap-quality 30
```

**Output**: `numts.txt`, `records_NC_010642.1.txt`, and `tiger_001_recovered.bam`

### Example 3: From Pre-aligned BAM

```bash
./numts_pipeline.py \
    -i /data/tiger/aligned.bam \
    -r /ref/tiger_genome.fasta \
    -m NC_010642.1 \
    -o /results/tiger_numts/ \
    -s tiger_001
```

**Output**: `numts.txt` and `records_NC_010642.1.txt`

### Example 4: Dry Run to Check Pipeline

```bash
./numts_pipeline.py \
    -i /data/tiger/reads/ \
    -r /ref/tiger_genome.fasta \
    -m NC_010642.1 \
    -o /results/tiger_numts/ \
    --dry-run
```

This will show which rules would be executed without running them.

### Example 5: Custom Quality Settings for Ancient DNA

```bash
./numts_pipeline.py \
    -i /data/ancient_tiger/reads/ \
    -r /ref/tiger_genome.fasta \
    -m NC_010642.1 \
    -o /results/ancient_tiger/ \
    -s ancient_tiger_001 \
    -t 16 \
    --min-quality 20 \
    --min-length 25 \
    --bwa-mismatch-penalty 0.05 \
    --remap
```

## Understanding the Output

### numts.txt

Each line contains a read ID that was identified as a NUMTS. These are reads that:
- Map to both mitochondrial and nuclear genomes
- Map better to the nuclear genome than to the mitochondrial genome

### records_<MT_NAME>.txt

Contains four lines with summary statistics:
```
No. of sequences binding only to the MT genome: <count>
No. of sequences binding once to the MT and once or more to the nuclear genome: <count>
No. of numts: <count>
No. of sequences mapping better to MT than to nucl: <count>
```

### Recovered BAM (if --remap specified)

Contains only the sequences that map well to the mitochondrial genome after NUMTS removal. Useful for downstream mitochondrial genome analysis.

## Troubleshooting

### Error: "No FASTQ or BAM files found"

**Solution**: Ensure your input directory contains properly named FASTQ files (_1/_2 or _s1/_s2 or _r1/_r2 suffixes) or specify files directly.

### Error: "Reference already indexed, skipping..."

**Info**: This is normal - the pipeline detected existing BWA index files.

### Warning: "Using single-stranded read mode"

**Info**: Only one FASTQ file was found. Pipeline will proceed with single-end mode. Consider checking if R2 file is missing.

### Low NUMTS Count

**Possible causes**:
- Input BAM was filtered with high mapping quality
- Very few multi-mapping reads in the data
- Wrong mitochondrial scaffold name specified

### Pipeline Fails at Alignment Step

**Solutions**:
- Check that BWA is installed and in PATH
- Verify reference genome is valid FASTA format
- Ensure sufficient disk space for intermediate files

## Citation

If you use this pipeline, please cite:

```
FragmentsMappingTwicePY.py - Written by Kirstine Tersbøl Melsen, s215096
Pipeline wrapper - Generated for NUMTS analysis project
```

## Contact

For issues, questions, or suggestions, please contact the maintainer or open an issue on the repository.

## License

[Specify your license here]

## Version History

- **v1.0.0** (2025-11-10): Initial release with full pipeline functionality
  - Input detection for DS/SS/BAM
  - Preprocessing with fastp
  - BWA alignment with configurable parameters
  - NUMTS identification
  - Optional remapping to MT genome
