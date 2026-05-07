#!/usr/bin/env python3
"""
Wrapper functions for the automated glycosylation pipeline.
These can be imported and used directly in Python scripts.
"""

import sys
import os
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any


def run_glyco_prep(
    input_pdb: str,
    output_dir: str,
    input_tsv: Optional[str] = None,
    input_glycosylator_tsv: Optional[str] = None,
    protein_residue_start: int = 1,
    rotate_atoms: str = "OD1,CG,ND2,HD22,HD21,HB2,HB3",
    fixed_atom: str = "CB",
    center_atom: str = "CA",
    radius: float = 30.0,
    rotation_step: int = 1,
) -> str:
    """
    Run glycosylation preparation step.
    
    Args:
        input_pdb: Input PDB file
        output_dir: Output directory
        input_tsv: Input TSV file (Caselino format)
        input_glycosylator_tsv: Pre-processed glycosylator TSV
        protein_residue_start: Protein residue start number
        rotate_atoms: Atoms to rotate
        fixed_atom: Fixed atom
        center_atom: Center atom
        radius: Radius for orientation
        rotation_step: Rotation step in degrees
    
    Returns:
        Path to final PDB file
    """
    from .cli import run_glyco_prep as _run_glyco_prep
    
    # Store original argv
    original_argv = sys.argv
    
    try:
        # Create fake argv
        sys.argv = [
            "glyco-prep",
            "--input-pdb", input_pdb,
            "--output-dir", output_dir
        ]
        
        if input_tsv:
            sys.argv.extend(["--input-tsv", input_tsv])
        if input_glycosylator_tsv:
            sys.argv.extend(["--input-glycosylator-tsv", input_glycosylator_tsv])
        
        sys.argv.extend([
            "--protein-residue-start", str(protein_residue_start),
            "--rotate-atoms", rotate_atoms,
            "--fixed-atom", fixed_atom,
            "--center-atom", center_atom,
            "--radius", str(radius),
            "--rotation-step", str(rotation_step)
        ])
        
        # Run the function
        result = _run_glyco_prep()
        return result
        
    finally:
        # Restore original argv
        sys.argv = original_argv


def run_glyco_param(
    prep_output_dir: str,
    input_pdb: str,
    output_dir: str,
    skip_charmm_download: bool = False,
) -> Optional[str]:
    """
    Run parametrization step.
    
    Args:
        prep_output_dir: Output directory from prep step
        input_pdb: Input PDB file from prep step
        output_dir: Output directory for parametrization
        skip_charmm_download: Skip CHARMM force field download
    
    Returns:
        Path to final PDB file or None
    """
    from .cli import run_glyco_param as _run_glyco_param
    
    original_argv = sys.argv
    
    try:
        sys.argv = [
            "glyco-param",
            "--prep-output-dir", prep_output_dir,
            "--input-pdb", input_pdb,
            "--output-dir", output_dir
        ]
        
        if skip_charmm_download:
            sys.argv.append("--skip-charmm-download")
        
        result = _run_glyco_param()
        return result
        
    finally:
        sys.argv = original_argv


def run_glyco_orient(
    input_pdb: str,
    param_output_dir: str,
    output_dir: str,
    charmm_dir: Optional[str] = None,
    theta_step: int = 10,
    n_steps: int = 10,
    max_cycles: int = 5,
    radius: float = 300,
    use_coulomb: str = "no",
    n_workers: int = 1,
    save_individual_glycans: bool = False,
    save_before_after: bool = False,
    verbose: bool = False,
) -> str:
    """
    Run carbohydrate orientation step.
    
    Args:
        input_pdb: Input PDB file from param step
        param_output_dir: Output directory from param step
        output_dir: Output directory for orientation
        charmm_dir: CHARMM36 directory
        theta_step: Theta step for MCMC
        n_steps: Number of steps for MCMC
        max_cycles: Maximum cycles for MCMC
        radius: Radius for orientation
        use_coulomb: Use Coulomb potential
        n_workers: Number of workers
        save_individual_glycans: Save individual glycans
        save_before_after: Save before/after structures
        verbose: Verbose output
    
    Returns:
        Path to optimized PDB file
    """
    from .cli import run_glyco_orient as _run_glyco_orient
    
    original_argv = sys.argv
    
    try:
        sys.argv = [
            "glyco-orient",
            "--input-pdb", input_pdb,
            "--param-output-dir", param_output_dir,
            "--output-dir", output_dir
        ]
        
        if charmm_dir:
            sys.argv.extend(["--charmm-dir", charmm_dir])
        
        sys.argv.extend([
            "--theta-step", str(theta_step),
            "--n-steps", str(n_steps),
            "--max-cycles", str(max_cycles),
            "--radius", str(radius),
            "--use-coulomb", use_coulomb,
            "--n-workers", str(n_workers)
        ])
        
        if save_individual_glycans:
            sys.argv.append("--save-individual-glycans")
        if save_before_after:
            sys.argv.append("--save-before-after")
        if verbose:
            sys.argv.append("--verbose")
        
        result = _run_glyco_orient()
        return result
        
    finally:
        sys.argv = original_argv


def run_all_pipeline(
    input_pdb: str,
    prep_output_dir: str,
    param_output_dir: str,
    orient_output_dir: str,
    input_tsv: Optional[str] = None,
    input_glycosylator_tsv: Optional[str] = None,
    protein_residue_start: int = 1,
    rotate_atoms: str = "OD1,CG,ND2,HD22,HD21,HB2,HB3",
    fixed_atom: str = "CB",
    center_atom: str = "CA",
    radius_prep: float = 30.0,
    rotation_step: int = 1,
    skip_charmm_download: bool = False,
    charmm_dir: Optional[str] = None,
    theta_step: int = 10,
    n_steps: int = 10,
    max_cycles: int = 5,
    radius_orient: float = 300,
    use_coulomb: str = "no",
    n_workers: int = 1,
    save_individual_glycans: bool = False,
    save_before_after: bool = False,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Run complete pipeline.
    
    Args:
        input_pdb: Input PDB file
        prep_output_dir: Output directory for prep step
        param_output_dir: Output directory for param step
        orient_output_dir: Output directory for orient step
        input_tsv: Input TSV file
        input_glycosylator_tsv: Pre-processed glycosylator TSV
        protein_residue_start: Protein residue start number
        rotate_atoms: Atoms to rotate
        fixed_atom: Fixed atom
        center_atom: Center atom
        radius_prep: Radius for prep orientation
        rotation_step: Rotation step in degrees
        skip_charmm_download: Skip CHARMM force field download
        charmm_dir: CHARMM36 directory
        theta_step: Theta step for MCMC
        n_steps: Number of steps for MCMC
        max_cycles: Maximum cycles for MCMC
        radius_orient: Radius for orientation
        use_coulomb: Use Coulomb potential
        n_workers: Number of workers
        save_individual_glycans: Save individual glycans
        save_before_after: Save before/after structures
        verbose: Verbose output
    
    Returns:
        Dictionary with results from each step
    """
    from .cli import run_glyco_all as _run_glyco_all
    
    original_argv = sys.argv
    
    try:
        sys.argv = [
            "glyco-all",
            "--input-pdb", input_pdb,
            "--prep-output-dir", prep_output_dir,
            "--param-output-dir", param_output_dir,
            "--orient-output-dir", orient_output_dir
        ]
        
        if input_tsv:
            sys.argv.extend(["--input-tsv", input_tsv])
        if input_glycosylator_tsv:
            sys.argv.extend(["--input-glycosylator-tsv", input_glycosylator_tsv])
        
        sys.argv.extend([
            "--protein-residue-start", str(protein_residue_start),
            "--rotate-atoms", rotate_atoms,
            "--fixed-atom", fixed_atom,
            "--center-atom", center_atom,
            "--radius-prep", str(radius_prep),
            "--rotation-step", str(rotation_step)
        ])
        
        if skip_charmm_download:
            sys.argv.append("--skip-charmm-download")
        
        if charmm_dir:
            sys.argv.extend(["--charmm-dir", charmm_dir])
        
        sys.argv.extend([
            "--theta-step", str(theta_step),
            "--n-steps", str(n_steps),
            "--max-cycles", str(max_cycles),
            "--radius-orient", str(radius_orient),
            "--use-coulomb", use_coulomb,
            "--n-workers", str(n_workers)
        ])
        
        if save_individual_glycans:
            sys.argv.append("--save-individual-glycans")
        if save_before_after:
            sys.argv.append("--save-before-after")
        if verbose:
            sys.argv.append("--verbose")
        
        result = _run_glyco_all()
        return result
        
    finally:
        sys.argv = original_argv
