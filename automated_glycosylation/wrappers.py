#!/usr/bin/env python3
"""
Wrapper functions for the automated glycosylation pipeline.
"""

import sys
import subprocess
from pathlib import Path

def get_package_dir() -> Path:
    """Get the package installation directory."""
    return Path(__file__).parent

def run_glyco_prep():
    """Wrapper for glycosylation preparation."""
    from .cli import run_glyco_prep as _run_glyco_prep
    _run_glyco_prep()

def run_glyco_param():
    """Wrapper for parametrization."""
    from .cli import run_glyco_param as _run_glyco_param
    _run_glyco_param()

def run_glyco_orient():
    """Wrapper for carbohydrate orientation."""
    from .cli import run_glyco_orient as _run_glyco_orient
    _run_glyco_orient()

def run_all_pipeline():
    """Wrapper for complete pipeline."""
    from .cli import run_all_pipeline as _run_all_pipeline
    _run_all_pipeline()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python -m automated_glycosylation [prep|param|orient|all] [arguments]")
        sys.exit(1)
    
    command = sys.argv[1]
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    
    if command == "prep":
        run_glyco_prep()
    elif command == "param":
        run_glyco_param()
    elif command == "orient":
        run_glyco_orient()
    elif command == "all":
        run_all_pipeline()
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
