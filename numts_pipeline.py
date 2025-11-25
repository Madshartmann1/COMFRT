#!/usr/bin/env python3
"""
NUMTS Recovery Pipeline Wrapper
Identifies Nuclear Mitochondrial DNA segments (NUMTs) from sequencing data.
"""

import os
import sys
import argparse
import yaml
import subprocess
from pathlib import Path


def show_help():
    """Show detailed help with examples"""
    help_text = """
numts_pipeline - Identify Nuclear Mitochondrial DNA segments (NUMTs)

DESCRIPTION:
    This pipeline identifies NUMTs by detecting reads that map to both 
    mitochondrial and nuclear genomes. It accepts paired-end reads, single-end 
    reads, or pre-aligned BAM files (with MAPQ=0 to retain multi-mapping reads).

USAGE:
    numts_pipeline [OPTIONS]

REQUIRED ARGUMENTS:
    -i, --input PATH           Input: BAM file, directory with FASTQ reads, 
                              or comma-separated R1,R2 files
    -r, --reference FILE      Reference genome FASTA (must contain MT scaffold)
    -m, --mitochondrial-scaffold STR
                              Name(s) of mitochondrial scaffold(s) 
                              (comma-separated if multiple)
    -o, --output-dir DIR      Output directory for results

OPTIONAL ARGUMENTS:
    -s, --sample-name STR     Sample name for output files [default: sample]
    -t, --threads INT         Number of threads/cores to use [default: 4]

PREPROCESSING OPTIONS:
    --min-quality INT         Minimum quality score for fastp [default: 25]
    --min-length INT          Minimum read length for fastp [default: 30]

BWA ALIGNMENT OPTIONS:
    --bwa-seed-length INT     BWA seed length (-l parameter) [default: 999]
    --bwa-max-gap-opens INT   BWA max gap opens (-o parameter) [default: 2]
    --bwa-mismatch-penalty FLOAT
                              BWA mismatch penalty (-n parameter) [default: 0.04]

PIPELINE OPTIONS:
    --remap                   Generate remapped BAM of recovered reads
    --remap-quality INT       Min mapping quality for remapping [default: 20]
    --keep-temp               Keep temporary intermediate files
    --dry-run                 Show pipeline without executing

TOOL PATHS (if not in PATH):
    --fastp-path PATH         Path to fastp executable [default: fastp]
    --bwa-path PATH           Path to BWA executable [default: bwa]
    --samtools-path PATH      Path to samtools executable [default: samtools]

    Note: If tools are not in your PATH, run ./validate_setup.py first.
          It will generate a tool_paths.yaml template for you to fill in.

    -h, --help                Show this help message

EXAMPLES:
    # Basic usage with paired-end reads
    numts_pipeline.py -i /data/reads/ -r genome.fasta -m NC_010642.1 -o output/
    
    # With pre-aligned BAM file
    numts_pipeline.py -i aligned.bam -r genome.fasta -m NC_010642.1 -o output/
    
    # Full analysis with remapping
    numts_pipeline.py -i /data/reads/ -r genome.fasta -m NC_010642.1 \\
        -o output/ -s tiger_001 -t 16 --remap --remap-quality 30
    
    # Multiple mitochondrial scaffolds
    numts_pipeline.py -i /data/reads/ -r genome.fasta \\
        -m "NC_010642.1,MT_scaffold2" -o output/
    
    # Ancient DNA with relaxed parameters
    numts_pipeline.py -i /data/ancient_reads/ -r genome.fasta -m MT \\
        -o output/ --min-quality 20 --min-length 25 \\
        --bwa-mismatch-penalty 0.05 --remap
    
    # Dry run to check pipeline
    numts_pipeline.py -i /data/reads/ -r genome.fasta -m MT -o output/ --dry-run

INPUT FORMATS:
    Directory:  Looks for paired-end reads with suffixes:
                _1/_2, _s1/_s2, _r1/_r2 (with .fq, .fastq, .fq.gz, .fastq.gz)
                
    BAM file:   Must be aligned with MAPQ=0 to retain multi-mapping reads
                (Pipeline will warn and confirm before proceeding)
                
    Files:      Comma-separated R1,R2 files (e.g., reads_1.fq.gz,reads_2.fq.gz)

OUTPUT FILES:
    numts.txt                      List of NUMTS read IDs
    records_<MT_NAME>.txt          Summary statistics
    <sample>_recovered.bam         Remapped non-NUMTS reads (if --remap used)
    logs/                          Log files for each step
    stats/<sample>_mapping_stats.txt  Mapping statistics

IMPORTANT NOTES:
    * BAM files MUST be aligned with MAPQ=0 for NUMTS detection
    * Mitochondrial scaffold must exist in the reference genome
    * For ancient DNA, use relaxed quality parameters
    * Use --remap to get clean mitochondrial BAM for downstream analysis

For more details, see README.md
"""
    print(help_text)


def check_file_exists(filepath, name):
    """Check if a file exists and raise error if not."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{name} not found: {filepath}")
    return filepath


def detect_input_type(input_path):
    """
    Detect input type: bam file, ds (double-stranded) reads, or ss (single-stranded) reads.
    Returns: ('bam', filepath) or ('ds', [r1, r2]) or ('ss', filepath)
    """
    if os.path.isfile(input_path):
        # Single file - check if it's a BAM
        if input_path.endswith('.bam'):
            return 'bam', input_path
        else:
            raise ValueError(f"Single file provided but not a BAM file: {input_path}")
    
    elif os.path.isdir(input_path):
        # Directory - look for paired-end reads
        files = os.listdir(input_path)
        
        # Check for different paired-end naming conventions
        r1_patterns = ['_1.fq', '_1.fastq', '_s1.fq', '_s1.fastq', '_r1.fq', '_r1.fastq',
                      '_1.fq.gz', '_1.fastq.gz', '_s1.fq.gz', '_s1.fastq.gz', '_r1.fq.gz', '_r1.fastq.gz']
        r2_patterns = ['_2.fq', '_2.fastq', '_s2.fq', '_s2.fastq', '_r2.fq', '_r2.fastq',
                      '_2.fq.gz', '_2.fastq.gz', '_s2.fq.gz', '_s2.fastq.gz', '_r2.fq.gz', '_r2.fastq.gz']
        
        r1_files = []
        r2_files = []
        
        for f in files:
            for pattern in r1_patterns:
                if f.endswith(pattern):
                    r1_files.append(os.path.join(input_path, f))
                    break
        
        for f in files:
            for pattern in r2_patterns:
                if f.endswith(pattern):
                    r2_files.append(os.path.join(input_path, f))
                    break
        
        if r1_files and r2_files:
            # Found paired-end reads
            if len(r1_files) == 1 and len(r2_files) == 1:
                return 'ds', [r1_files[0], r2_files[0]]
            else:
                raise ValueError(f"Multiple R1/R2 files found in {input_path}. Please specify exact files.")
        
        # Look for single-end reads
        fastq_files = [os.path.join(input_path, f) for f in files 
                      if f.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz'))]
        
        if len(fastq_files) == 1:
            print(f"WARNING: No paired-end reads found. Using single-end mode with: {fastq_files[0]}")
            return 'ss', fastq_files[0]
        elif len(fastq_files) > 1:
            raise ValueError(f"Multiple FASTQ files found but couldn't identify pairs. Please specify exact files.")
        else:
            raise ValueError(f"No FASTQ or BAM files found in {input_path}")
    
    else:
        raise ValueError(f"Input path does not exist: {input_path}")


def parse_arguments():
    """Parse command line arguments."""
    # Custom help handling
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ['-h', '--help']):
        show_help()
        sys.exit(0)
    
    parser = argparse.ArgumentParser(
        prog='numts_pipeline',
        description='NUMTS Recovery Pipeline - Identifies Nuclear Mitochondrial DNA segments',
        add_help=False  # We handle help ourselves
    )
    
    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='PATH',
                       help='Input: BAM file, directory with FASTQ reads, or comma-separated R1,R2')
    required.add_argument('-r', '--reference', required=True, metavar='FILE',
                       help='Reference genome FASTA file (must contain mitochondrial scaffold)')
    required.add_argument('-m', '--mitochondrial-scaffold', required=True, metavar='STR',
                       help='Name(s) of mitochondrial scaffold(s) (comma-separated if multiple)')
    required.add_argument('-o', '--output-dir', required=True, metavar='DIR',
                       help='Output directory for results')
    
    # Optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-s', '--sample-name', default='sample', metavar='STR',
                       help='Sample name for output files (default: %(default)s)')
    optional.add_argument('-t', '--threads', type=int, default=4, metavar='INT',
                       help='Total number of threads/cores to use (default: %(default)s)')
    
    # Preprocessing options
    preproc = parser.add_argument_group('Preprocessing options')
    preproc.add_argument('--min-quality', type=int, default=25, metavar='INT',
                       help='Minimum quality score for fastp (default: %(default)s)')
    preproc.add_argument('--min-length', type=int, default=30, metavar='INT',
                       help='Minimum read length for fastp (default: %(default)s)')
    
    # BWA alignment options
    bwa_opts = parser.add_argument_group('BWA alignment options')
    bwa_opts.add_argument('--bwa-seed-length', type=int, default=999, metavar='INT',
                       help='BWA seed length (-l parameter) (default: %(default)s)')
    bwa_opts.add_argument('--bwa-max-gap-opens', type=int, default=2, metavar='INT',
                       help='BWA maximum gap opens (-o parameter) (default: %(default)s)')
    bwa_opts.add_argument('--bwa-mismatch-penalty', type=float, default=0.04, metavar='FLOAT',
                       help='BWA mismatch penalty (-n parameter) (default: %(default)s)')
    
    # Pipeline options
    pipeline = parser.add_argument_group('Pipeline options')
    pipeline.add_argument('--remap', action='store_true',
                       help='Generate remapped BAM file of recovered reads')
    pipeline.add_argument('--remap-quality', type=int, default=20, metavar='INT',
                       help='Minimum mapping quality for remapping step (default: %(default)s)')
    pipeline.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files')
    pipeline.add_argument('--dry-run', action='store_true',
                       help='Perform Snakemake dry run only')
    
    # Tool paths
    tools = parser.add_argument_group('Tool paths (if not in PATH)')
    tools.add_argument('--fastp-path', default='fastp', metavar='PATH',
                       help='Path to fastp executable (default: %(default)s)')
    tools.add_argument('--bwa-path', default='bwa', metavar='PATH',
                       help='Path to BWA executable (default: %(default)s)')
    tools.add_argument('--samtools-path', default='samtools', metavar='PATH',
                       help='Path to samtools executable (default: %(default)s)')
    
    # Help
    parser.add_argument('-h', '--help', action='store_true',
                       help='Show detailed help message')
    
    args = parser.parse_args()
    
    # Handle help
    if args.help:
        show_help()
        sys.exit(0)
    
    return args


def load_tool_paths(script_dir):
    """Load tool paths from tool_paths.yaml if it exists."""
    # Check in scripts directory first, then in main directory
    for location in [os.path.join(script_dir, 'scripts'), script_dir]:
        tool_paths_file = os.path.join(location, 'tool_paths.yaml')
        
        if os.path.exists(tool_paths_file):
            try:
                with open(tool_paths_file, 'r') as f:
                    data = yaml.safe_load(f)
                    if data and 'tool_paths' in data:
                        print(f"Loading tool paths from: {tool_paths_file}")
                        return data['tool_paths']
            except Exception as e:
                print(f"Warning: Could not read {tool_paths_file}: {e}")
    
    return {}


def validate_mt_scaffolds(reference_file, mt_scaffolds):
    """Check that mitochondrial scaffold names exist in the reference file."""
    print(f"Validating mitochondrial scaffold names in reference...")
    
    # Read scaffold names from reference FASTA
    ref_scaffolds = set()
    try:
        with open(reference_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract scaffold name (first word after >)
                    scaffold_name = line[1:].split()[0]
                    ref_scaffolds.add(scaffold_name)
    except Exception as e:
        raise ValueError(f"Error reading reference file: {e}")
    
    if not ref_scaffolds:
        raise ValueError(f"No sequences found in reference file: {reference_file}")
    
    # Check each MT scaffold
    missing_scaffolds = []
    for mt_scaffold in mt_scaffolds:
        if mt_scaffold not in ref_scaffolds:
            missing_scaffolds.append(mt_scaffold)
    
    if missing_scaffolds:
        print("\n" + "="*80)
        print("ERROR: Mitochondrial scaffold(s) not found in reference file!")
        print(f"Missing: {', '.join(missing_scaffolds)}")
        print("\nAvailable scaffolds in reference:")
        for scaffold in sorted(ref_scaffolds)[:20]:  # Show first 20
            print(f"  - {scaffold}")
        if len(ref_scaffolds) > 20:
            print(f"  ... and {len(ref_scaffolds) - 20} more")
        print("="*80 + "\n")
        raise ValueError(f"Mitochondrial scaffolds not found in reference: {', '.join(missing_scaffolds)}")
    
    print(f"✓ All mitochondrial scaffolds found in reference: {', '.join(mt_scaffolds)}")
    return True


def create_config(args):
    """Create configuration dictionary for Snakemake."""
    
    # Detect input type
    input_type, input_data = detect_input_type(args.input)
    
    # Parse mitochondrial scaffold names
    mt_scaffolds = [s.strip() for s in args.mitochondrial_scaffold.split(',')]
    
    # Validate that MT scaffolds exist in reference
    validate_mt_scaffolds(args.reference, mt_scaffolds)
    
    # Load tool paths from file if available
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tool_paths_from_file = load_tool_paths(script_dir)
    
    # Load tool paths from file if available
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tool_paths_from_file = load_tool_paths(script_dir)
    
    # Tool paths: prioritize command line args, then tool_paths.yaml, then defaults
    fastp_path = args.fastp_path if args.fastp_path != 'fastp' else tool_paths_from_file.get('fastp', 'fastp')
    bwa_path = args.bwa_path if args.bwa_path != 'bwa' else tool_paths_from_file.get('bwa', 'bwa')
    samtools_path = args.samtools_path if args.samtools_path != 'samtools' else tool_paths_from_file.get('samtools', 'samtools')
    
    # Show which tool paths are being used
    if tool_paths_from_file:
        print(f"  fastp: {fastp_path}")
        print(f"  bwa: {bwa_path}")
        print(f"  samtools: {samtools_path}")
    
    # Warn user if providing BAM file
    if input_type == 'bam':
        print("\n" + "="*80)
        print("WARNING: You are providing a pre-aligned BAM file.")
        print("For this pipeline to work correctly, the BAM file MUST have been aligned")
        print("with MINIMUM MAPPING QUALITY = 0 to retain multi-mapping reads.")
        print("If your BAM was filtered for high-quality mappings only, this pipeline")
        print("will NOT be able to identify NUMTs correctly.")
        print("="*80 + "\n")
        
        response = input("Continue with provided BAM file? (yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("Exiting...")
            sys.exit(0)
    
    # Warn if using single-end mode
    if input_type == 'ss':
        print("\n" + "="*80)
        print("WARNING: Using SINGLE-STRANDED read mode.")
        print("This may reduce sensitivity for NUMTS detection.")
        print("="*80 + "\n")
    
    # Create config dictionary
    config = {
        'sample_name': args.sample_name,
        'output_dir': os.path.abspath(args.output_dir),
        'threads': args.threads,
        
        # Input information
        'input_type': input_type,  # 'bam', 'ds', or 'ss'
        'reference': os.path.abspath(args.reference),
        'mt_scaffolds': mt_scaffolds,
        
        # Preprocessing parameters
        'min_quality': args.min_quality,
        'min_length': args.min_length,
        
        # BWA parameters
        'bwa_seed_length': args.bwa_seed_length,
        'bwa_max_gap_opens': args.bwa_max_gap_opens,
        'bwa_mismatch_penalty': args.bwa_mismatch_penalty,
        
        # Pipeline options
        'remap': args.remap,
        'remap_quality': args.remap_quality,
        'keep_temp': args.keep_temp,
        
        # Tool paths
        'fastp': fastp_path,
        'bwa': bwa_path,
        'samtools': samtools_path,
        
        # Script paths (point to scripts/ subdirectory where helper scripts are located)
        'script_dir': os.path.join(script_dir, 'scripts'),
    }
    
    # Add input-specific information
    if input_type == 'bam':
        config['input_bam'] = os.path.abspath(input_data)
    elif input_type == 'ds':
        config['input_r1'] = os.path.abspath(input_data[0])
        config['input_r2'] = os.path.abspath(input_data[1])
    elif input_type == 'ss':
        config['input_fastq'] = os.path.abspath(input_data)
    
    return config


def write_config_file(config, output_path):
    """Write config dictionary to YAML file."""
    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    print(f"Configuration file written to: {output_path}")


def run_snakemake(config_file, dry_run=False, threads=4):
    """Execute Snakemake pipeline."""
    
    # Get path to Snakefile (should be in scripts directory)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    snakefile = os.path.join(script_dir, 'scripts', 'Snakefile_numts')
    
    if not os.path.exists(snakefile):
        raise FileNotFoundError(f"Snakefile not found: {snakefile}")
    
    # Build Snakemake command
    cmd = [
        'snakemake',
        '--snakefile', snakefile,
        '--configfile', config_file,
        '--cores', str(threads),
        '--printshellcmds',
    ]
    
    if dry_run:
        cmd.append('--dry-run')
    
    print("\n" + "="*80)
    print("Executing Snakemake pipeline...")
    print("Command:", ' '.join(cmd))
    print("="*80 + "\n")
    
    # Run Snakemake
    try:
        subprocess.run(cmd, check=True)
        print("\n" + "="*80)
        print("Pipeline completed successfully!")
        print("="*80 + "\n")
    except subprocess.CalledProcessError as e:
        print("\n" + "="*80)
        print(f"Pipeline failed with error code {e.returncode}")
        print("="*80 + "\n")
        sys.exit(1)


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create configuration
    config = create_config(args)
    
    # Write config file
    config_file = os.path.join(args.output_dir, 'config.yaml')
    write_config_file(config, config_file)
    
    # Run Snakemake
    run_snakemake(config_file, dry_run=args.dry_run, threads=args.threads)
    
    print(f"\nResults are in: {args.output_dir}")
    print(f"NUMTS list: {os.path.join(args.output_dir, 'numts.txt')}")
    print(f"Records: {os.path.join(args.output_dir, f'records_{config['mt_scaffolds'][0]}.txt')}")
    
    if args.remap:
        print(f"Remapped BAM: {os.path.join(args.output_dir, f'{args.sample_name}_recovered.bam')}")


if __name__ == '__main__':
    main()
