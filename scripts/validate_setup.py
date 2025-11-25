#!/usr/bin/env python3
"""
Validation script to check if all dependencies are installed correctly
for the NUMTS pipeline.
"""

import subprocess
import sys
import os
from pathlib import Path


def check_command(cmd, name, version_flag='--version'):
    """Check if a command is available and print version."""
    try:
        result = subprocess.run([cmd, version_flag], 
                              capture_output=True, 
                              text=True, 
                              timeout=5)
        if result.returncode == 0 or 'version' in result.stdout.lower() or 'version' in result.stderr.lower():
            version_output = result.stdout or result.stderr
            first_line = version_output.split('\n')[0]
            print(f"✓ {name:20s} FOUND: {first_line[:60]}")
            return True
        else:
            print(f"✗ {name:20s} FOUND but version check failed")
            return False
    except FileNotFoundError:
        print(f"✗ {name:20s} NOT FOUND - please install")
        return False
    except subprocess.TimeoutExpired:
        print(f"✗ {name:20s} TIMEOUT - command hung")
        return False
    except Exception as e:
        print(f"✗ {name:20s} ERROR: {str(e)}")
        return False


def check_python_module(module_name, import_name=None):
    """Check if a Python module is available."""
    if import_name is None:
        import_name = module_name
    
    try:
        __import__(import_name)
        print(f"✓ Python: {module_name:13s} FOUND")
        return True
    except ImportError:
        print(f"✗ Python: {module_name:13s} NOT FOUND - install with: pip install {module_name}")
        return False


def check_file_exists(filepath, description):
    """Check if a file exists."""
    if os.path.exists(filepath):
        print(f"✓ {description:20s} FOUND: {filepath}")
        return True
    else:
        print(f"✗ {description:20s} NOT FOUND: {filepath}")
        return False


def generate_tool_paths_template(missing_tools, script_dir):
    """Generate a template file for tool paths."""
    template_path = os.path.join(script_dir, 'tool_paths.yaml')
    
    if os.path.exists(template_path):
        print(f"\nNote: {template_path} already exists, not overwriting.")
        return template_path
    
    template_content = {
        'tool_paths': {
            'fastp': 'fastp',
            'bwa': 'bwa',
            'samtools': 'samtools',
        }
    }
    
    # Add comments for missing tools
    import yaml
    with open(template_path, 'w') as f:
        f.write("# Tool paths configuration for NUMTS pipeline\n")
        f.write("# Edit this file to provide full paths to tools not in your PATH\n")
        f.write("# Then run the pipeline with: --fastp-path $(grep fastp tool_paths.yaml | awk '{print $2}')\n")
        f.write("#\n")
        if missing_tools:
            f.write("# The following tools were NOT FOUND during validation:\n")
            for tool in missing_tools:
                f.write(f"#   - {tool}\n")
            f.write("#\n")
        f.write("# Format: Replace 'toolname' with the full path, e.g.:\n")
        f.write("#   fastp: /usr/local/bin/fastp\n")
        f.write("#   bwa: /home/user/software/bwa/bwa\n")
        f.write("#   samtools: /opt/samtools-1.15/bin/samtools\n")
        f.write("\n")
        yaml.dump(template_content, f, default_flow_style=False)
    
    return template_path


def main():
    """Run all validation checks."""
    print("=" * 80)
    print("NUMTS Pipeline - Dependency Validation")
    print("=" * 80)
    print()
    
    all_good = True
    missing_tools = []
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Check Python version
    print("Checking Python version...")
    python_version = sys.version_info
    if python_version.major >= 3 and python_version.minor >= 6:
        print(f"✓ Python version: {python_version.major}.{python_version.minor}.{python_version.micro}")
    else:
        print(f"✗ Python version too old: {python_version.major}.{python_version.minor}.{python_version.micro}")
        print("  Required: Python 3.6 or higher")
        all_good = False
    print()
    
    # Check command-line tools
    print("Checking command-line tools...")
    tools = [
        ('snakemake', 'Snakemake', '--version'),
        ('bwa', 'BWA', None),  # BWA doesn't have --version
        ('samtools', 'samtools', '--version'),
        ('fastp', 'fastp', '--version'),
    ]
    
    for cmd, name, version_flag in tools:
        if version_flag:
            result = check_command(cmd, name, version_flag)
        else:
            # For BWA, just check if command exists
            try:
                result = subprocess.run([cmd], 
                                      capture_output=True, 
                                      text=True, 
                                      timeout=2)
                # BWA returns non-zero but that's OK if it runs
                print(f"✓ {name:20s} FOUND")
                result = True
            except FileNotFoundError:
                print(f"✗ {name:20s} NOT FOUND - please install")
                result = False
                missing_tools.append(name)
            except Exception:
                # Command exists but returns error - that's fine for BWA
                print(f"✓ {name:20s} FOUND")
                result = True
        
        if not result and version_flag:
            missing_tools.append(name)
        
        all_good = all_good and result
    print()
    
    # Check Python modules
    print("Checking Python modules...")
    modules = [
        ('yaml', 'pyyaml'),
        ('numpy', 'numpy'),
        ('pandas', 'pandas'),
    ]
    
    for import_name, pip_name in modules:
        result = check_python_module(pip_name, import_name)
        all_good = all_good and result
    print()
    
    # Check pipeline files
    print("Checking pipeline files...")
    files = [
        ('numts_pipeline.py', 'Wrapper script'),
        ('Snakefile_numts', 'Snakefile'),
        ('FragmentsMappingTwicePY.py', 'NUMTS detection script'),
        ('README.md', 'Documentation'),
    ]
    
    for filename, description in files:
        filepath = os.path.join(script_dir, filename)
        result = check_file_exists(filepath, description)
        all_good = all_good and result
    print()
    
    # Summary
    print("=" * 80)
    if all_good:
        print("✓ ALL CHECKS PASSED - Pipeline is ready to use!")
        print()
        print("Quick start:")
        print("  ./numts_pipeline.py -i INPUT -r REFERENCE.fasta -m MT_SCAFFOLD -o OUTPUT")
        print()
        print("For more information:")
        print("  ./numts_pipeline.py --help")
        print("  cat README.md")
        print("  cat QUICKSTART.md")
        return 0
    else:
        print("✗ SOME CHECKS FAILED - Please install missing dependencies")
        print()
        
        # Generate tool paths template if tools are missing
        if missing_tools:
            template_path = generate_tool_paths_template(missing_tools, script_dir)
            print(f"Generated tool paths template: {template_path}")
            print()
            print("To use custom tool paths:")
            print("  1. Edit tool_paths.yaml with full paths to your tools")
            print("  2. Run pipeline with tool path options, for example:")
            print("     ./numts_pipeline.py -i INPUT -r REF.fa -m MT -o OUT \\")
            print("       --bwa-path /path/to/bwa \\")
            print("       --samtools-path /path/to/samtools \\")
            print("       --fastp-path /path/to/fastp")
            print()
        
        print("Installation help:")
        print("  - Python packages: pip install pyyaml numpy pandas snakemake")
        print("  - BWA: https://github.com/lh3/bwa")
        print("  - samtools: http://www.htslib.org/")
        print("  - fastp: https://github.com/OpenGene/fastp")
        print()
        print("Or use conda:")
        print("  conda install -c bioconda bwa samtools fastp snakemake")
        print("  pip install pyyaml numpy pandas")
        return 1


if __name__ == '__main__':
    sys.exit(main())
