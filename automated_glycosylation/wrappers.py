#!/usr/bin/env python3
"""
...
"""

import sys
import subprocess
import argparse
from pathlib import Path

def get_scripts_dir():
    """Get the directory containing the scripts"""
    package_root = Path(__file__).parent
    return package_root / "scripts"

def get_bin_dir():
    """Get the bin directory"""
    package_root = Path(__file__).parent
    return package_root / "bin"

def run_all_analysis():
    """Wrapper for the complete pipeline"""
    bin_dir = get_bin_dir()
    script_path = bin_dir / "run_all_analysis.sh"
    
    # Make sure the script is executable
    if script_path.exists():
        script_path.chmod(0o755)
        # Pass all arguments to the script
        cmd = [str(script_path)] + sys.argv[1:]
        subprocess.run(cmd)
    else:
        print(f"Error: Pipeline script not found at {script_path}", file=sys.stderr)
        print("Please ensure bin/run_all_analysis.sh exists", file=sys.stderr)
        sys.exit(1)

...

if __name__ == "__main__":
...
