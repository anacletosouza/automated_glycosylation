#!/usr/bin/env python3
"""
Glycosylation Pipeline CLI - Unified interface for all steps
Usage: glyco-{prep,param,orient,all} [arguments]
"""

import sys
import os
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).parent.absolute()
SCRIPTS_DIR = SCRIPT_DIR / "scripts"


def find_script(step, script_name):
    """Find script in appropriate subdirectory"""
    step_dirs = {
        'prep': '1-glycosylation_preparation',
        'param': '2-parametrization_scripts',
        'orient': '3-carbohydrate_orientation'
    }
    
    if step in step_dirs:
        script_path = SCRIPTS_DIR / step_dirs[step] / script_name
        if script_path.exists():
            return str(script_path)
    
    # Search in all subdirectories
    for subdir in SCRIPTS_DIR.iterdir():
        if subdir.is_dir():
            script_path = subdir / script_name
            if script_path.exists():
                return str(script_path)
    
    raise FileNotFoundError(f"Script not found: {script_name}")


def get_abs_path(path):
    """Convert to absolute path"""
    if path is None:
        return None
    return str(Path(path).absolute())


def get_base_filename(filepath):
    """Get base filename without extension for generic naming"""
    if filepath:
        return Path(filepath).stem
    return "input"


def run_glyco_prep():
    """Run glycosylation preparation (Step 1)"""
    parser = argparse.ArgumentParser(
        description="Step 1: Glycosylation preparation",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--input-pdb", required=True, help="Input PDB file")
    parser.add_argument("--input-tsv", help="Input TSV file (Caselino format)")
    parser.add_argument("--input-glycosylator-tsv", help="Pre-processed glycosylator TSV")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--protein-residue-start", type=int, default=1, 
                       help="Protein residue start number (default: 1)")
    parser.add_argument("--rotate-atoms", default="OD1,CG,ND2,HD22,HD21,HB2,HB3",
                       help="Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)")
    parser.add_argument("--fixed-atom", default="CB", help="Fixed atom (default: CB)")
    parser.add_argument("--center-atom", default="CA", help="Center atom (default: CA)")
    parser.add_argument("--radius", type=float, default=30.0, 
                       help="Radius for orientation (default: 30.0)")
    parser.add_argument("--rotation-step", type=float, default=1.0,
                       help="Rotation step in degrees (default: 1)")
    
    args = parser.parse_args()
    
    # Find the main preparation script
    prep_script = find_script('prep', 'asn_orientation.py')
    
    # Build command
    cmd = [
        sys.executable, prep_script,
        "--input-pdb", get_abs_path(args.input_pdb),
        "--output-dir", get_abs_path(args.output_dir),
        "--protein-residue-start", str(args.protein_residue_start),
        "--rotate-atoms", args.rotate_atoms,
        "--fixed-atom", args.fixed_atom,
        "--center-atom", args.center_atom,
        "--radius", str(args.radius),
        "--rotation-step", str(args.rotation_step)
    ]
    
    if args.input_tsv:
        cmd.extend(["--input-tsv", get_abs_path(args.input_tsv)])
    if args.input_glycosylator_tsv:
        cmd.extend(["--input-glycosylator-tsv", get_abs_path(args.input_glycosylator_tsv)])
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


def run_glyco_param():
    """Run parametrization (Step 2)"""
    parser = argparse.ArgumentParser(
        description="Step 2: Parametrization",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--prep-output-dir", required=True, 
                       help="Output directory from preparation step")
    parser.add_argument("--input-pdb", required=True, help="Input PDB file from preparation")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--skip-charmm-download", action="store_true",
                       help="Skip CHARMM force field download")
    
    args = parser.parse_args()
    
    # Find the parametrization script
    param_script = find_script('param', 'run_all_2.sh')
    
    # Build command
    cmd = [
        "bash", param_script,
        "--prep-output-dir", get_abs_path(args.prep_output_dir),
        "--input-pdb", get_abs_path(args.input_pdb),
        "--output-dir", get_abs_path(args.output_dir)
    ]
    
    if args.skip_charmm_download:
        cmd.append("--skip-charmm-download")
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


def run_glyco_orient():
    """Run carbohydrate orientation (Step 3)"""
    parser = argparse.ArgumentParser(
        description="Step 3: Carbohydrate orientation optimization",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--input-pdb", required=True, 
                       help="Input PDB file from parametrization")
    parser.add_argument("--param-output-dir", required=True,
                       help="Output directory from parametrization step")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--charmm-dir", help="CHARMM36 directory (optional)")
    parser.add_argument("--theta-step", type=int, default=10,
                       help="Theta step for MCMC (default: 10)")
    parser.add_argument("--n-steps", type=int, default=10,
                       help="Number of steps for MCMC (default: 10)")
    parser.add_argument("--max-cycles", type=int, default=5,
                       help="Maximum cycles for MCMC (default: 5)")
    parser.add_argument("--radius", type=float, default=300.0,
                       help="Radius for orientation (default: 300)")
    parser.add_argument("--use-coulomb", choices=['yes', 'no'], default='no',
                       help="Use Coulomb potential (default: no)")
    parser.add_argument("--n-workers", type=int, default=1,
                       help="Number of workers (default: 1)")
    parser.add_argument("--save-individual-glycans", action="store_true",
                       help="Save individual glycans")
    parser.add_argument("--save-before-after", action="store_true",
                       help="Save before/after structures")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    # Find the orientation script
    orient_script = find_script('orient', 'run_all_3.sh')
    
    # Build command
    cmd = [
        "bash", orient_script,
        "--input-pdb", get_abs_path(args.input_pdb),
        "--param-output-dir", get_abs_path(args.param_output_dir),
        "--output-dir", get_abs_path(args.output_dir),
        "--theta-step", str(args.theta_step),
        "--n-steps", str(args.n_steps),
        "--max-cycles", str(args.max_cycles),
        "--radius", str(args.radius),
        "--use-coulomb", args.use_coulomb,
        "--n-workers", str(args.n_workers)
    ]
    
    if args.charmm_dir:
        cmd.extend(["--charmm-dir", get_abs_path(args.charmm_dir)])
    if args.save_individual_glycans:
        cmd.append("--save-individual-glycans")
    if args.save_before_after:
        cmd.append("--save-before-after")
    if args.verbose:
        cmd.append("--verbose")
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


def auto_glyco():
    """Run complete pipeline (Steps 1-3)"""
    parser = argparse.ArgumentParser(
        description="Complete glycosylation pipeline (Steps 1-3)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  auto_glyco --input-pdb protein.pdb \\
             --prep-output-dir ./prep \\
             --param-output-dir ./param \\
             --orient-output-dir ./orient

  # Run with TSV file and custom parameters
  auto_glyco --input-pdb protein.pdb \\
             --input-tsv table.tsv \\
             --prep-output-dir ./prep \\
             --param-output-dir ./param \\
             --orient-output-dir ./orient \\
             --protein-residue-start 10 \\
             --theta-step 20 \\
             --n-steps 50

  # Run only preparation
  auto_glyco --prep-only --input-pdb protein.pdb --prep-output-dir ./prep

  # Run only orientation with custom settings
  auto_glyco --orient-only \\
             --input-pdb protein_param.pdb \\
             --param-output-dir ./param \\
             --orient-output-dir ./orient \\
             --verbose
        """
    )
    
    # Required options
    parser.add_argument("--input-pdb", help="Input PDB file")
    parser.add_argument("--prep-output-dir", help="Output directory for preparation step")
    parser.add_argument("--param-output-dir", help="Output directory for parametrization step")
    parser.add_argument("--orient-output-dir", help="Output directory for orientation step")
    
    # Input options
    parser.add_argument("--input-tsv", help="Input TSV file (Caselino format)")
    parser.add_argument("--input-glycosylator-tsv", help="Pre-processed glycosylator TSV")
    
    # Preparation options
    parser.add_argument("--protein-residue-start", type=int, default=1,
                       help="Protein residue start number (default: 1)")
    parser.add_argument("--rotate-atoms", default="OD1,CG,ND2,HD22,HD21,HB2,HB3",
                       help="Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)")
    parser.add_argument("--fixed-atom", default="CB", help="Fixed atom (default: CB)")
    parser.add_argument("--center-atom", default="CA", help="Center atom (default: CA)")
    parser.add_argument("--radius-prep", type=float, default=30.0,
                       help="Radius for orientation (default: 30.0)")
    parser.add_argument("--rotation-step", type=float, default=1.0,
                       help="Rotation step in degrees (default: 1)")
    
    # Parametrization options
    parser.add_argument("--skip-charmm-download", action="store_true",
                       help="Skip CHARMM force field download")
    
    # Orientation options
    parser.add_argument("--charmm-dir", help="CHARMM36 directory (optional)")
    parser.add_argument("--theta-step", type=int, default=10,
                       help="Theta step for MCMC (default: 10)")
    parser.add_argument("--n-steps", type=int, default=10,
                       help="Number of steps for MCMC (default: 10)")
    parser.add_argument("--max-cycles", type=int, default=5,
                       help="Maximum cycles for MCMC (default: 5)")
    parser.add_argument("--radius-orient", type=float, default=300.0,
                       help="Radius for orientation (default: 300)")
    parser.add_argument("--use-coulomb", choices=['yes', 'no'], default='no',
                       help="Use Coulomb potential (default: no)")
    parser.add_argument("--n-workers", type=int, default=1,
                       help="Number of workers (default: 1)")
    parser.add_argument("--save-individual-glycans", action="store_true",
                       help="Save individual glycans")
    parser.add_argument("--save-before-after", action="store_true",
                       help="Save before/after structures")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    
    # Step selection
    parser.add_argument("--prep-only", action="store_true", help="Run only preparation step")
    parser.add_argument("--param-only", action="store_true", help="Run only parametrization step")
    parser.add_argument("--orient-only", action="store_true", help="Run only orientation step")
    
    args = parser.parse_args()
    
    # Determine which steps to run
    run_prep = not (args.param_only or args.orient_only) or args.prep_only
    run_param = not (args.prep_only or args.orient_only) or args.param_only
    run_orient = not (args.prep_only or args.param_only) or args.orient_only
    
    # Validate required arguments
    if any([run_prep, run_param, run_orient]) and not args.input_pdb:
        print("ERROR: --input-pdb is required", file=sys.stderr)
        sys.exit(1)
    
    if run_prep and not args.prep_output_dir:
        print("ERROR: --prep-output-dir is required for preparation step", file=sys.stderr)
        sys.exit(1)
    
    if run_param and not args.param_output_dir:
        print("ERROR: --param-output-dir is required for parametrization step", file=sys.stderr)
        sys.exit(1)
    
    if run_orient and not args.orient_output_dir:
        print("ERROR: --orient-output-dir is required for orientation step", file=sys.stderr)
        sys.exit(1)
    
    print("")
    print("=" * 60)
    print("AUTOMATED GLYCOSYLATION PIPELINE")
    print("=" * 60)
    print(f"Start time: {datetime.now()}")
    print("")
    
    # Track final PDBs
    final_prep_pdb = None
    final_param_pdb = None
    
    # Run preparation step if selected
    if run_prep:
        print("")
        print("=" * 60)
        print("RUNNING PREPARATION STEP")
        print("=" * 60)
        
        prep_script = find_script('prep', 'asn_orientation.py')
        cmd = [
            sys.executable, prep_script,
            "--input-pdb", get_abs_path(args.input_pdb),
            "--output-dir", get_abs_path(args.prep_output_dir),
            "--protein-residue-start", str(args.protein_residue_start),
            "--rotate-atoms", args.rotate_atoms,
            "--fixed-atom", args.fixed_atom,
            "--center-atom", args.center_atom,
            "--radius", str(args.radius_prep),
            "--rotation-step", str(args.rotation_step)
        ]
        
        if args.input_tsv:
            cmd.extend(["--input-tsv", get_abs_path(args.input_tsv)])
        if args.input_glycosylator_tsv:
            cmd.extend(["--input-glycosylator-tsv", get_abs_path(args.input_glycosylator_tsv)])
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)
        
        # Get the final PDB from preparation
        final_prep_pdb = Path(args.prep_output_dir) / "PDB_PROTEIN_GLYCOSYLATED" / "protein_renumbered.pdb"
        if not final_prep_pdb.exists():
            final_prep_pdb = Path(args.prep_output_dir) / "PDB_PROTEIN_GLYCOSYLATED" / "protein_asn_orientation.pdb"
        
        print(f"Preparation completed. Final PDB: {final_prep_pdb}")
    
    # Run parametrization step if selected
    if run_param:
        print("")
        print("=" * 60)
        print("RUNNING PARAMETRIZATION STEP")
        print("=" * 60)
        
        input_for_param = final_prep_pdb if final_prep_pdb else Path(args.input_pdb)
        
        param_script = find_script('param', 'run_all_2.sh')
        cmd = [
            "bash", param_script,
            "--prep-output-dir", get_abs_path(args.prep_output_dir) if args.prep_output_dir else ".",
            "--input-pdb", str(input_for_param),
            "--output-dir", get_abs_path(args.param_output_dir)
        ]
        
        if args.skip_charmm_download:
            cmd.append("--skip-charmm-download")
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)
        
        # Get the final PDB from parametrization
        final_param_pdb = Path(args.param_output_dir) / "PDB_GLYCOPROTEIN" / "protein_final_valence_corrected.pdb"
        
        print(f"Parametrization completed. Final PDB: {final_param_pdb}")
    
    # Run orientation step if selected
    if run_orient:
        print("")
        print("=" * 60)
        print("RUNNING ORIENTATION STEP")
        print("=" * 60)
        
        input_for_orient = final_param_pdb if final_param_pdb else Path(args.input_pdb)
        
        orient_script = find_script('orient', 'run_all_3.sh')
        cmd = [
            "bash", orient_script,
            "--input-pdb", str(input_for_orient),
            "--param-output-dir", get_abs_path(args.param_output_dir) if args.param_output_dir else ".",
            "--output-dir", get_abs_path(args.orient_output_dir),
            "--theta-step", str(args.theta_step),
            "--n-steps", str(args.n_steps),
            "--max-cycles", str(args.max_cycles),
            "--radius", str(args.radius_orient),
            "--use-coulomb", args.use_coulomb,
            "--n-workers", str(args.n_workers)
        ]
        
        if args.charmm_dir:
            cmd.extend(["--charmm-dir", get_abs_path(args.charmm_dir)])
        if args.save_individual_glycans:
            cmd.append("--save-individual-glycans")
        if args.save_before_after:
            cmd.append("--save-before-after")
        if args.verbose:
            cmd.append("--verbose")
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)
        
        final_orient_pdb = Path(args.orient_output_dir) / "PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED" / "protein_optimized.pdb"
        print(f"Orientation completed. Final PDB: {final_orient_pdb}")
    
    # Summary
    print("")
    print("=" * 60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 60)
    print(f"End time: {datetime.now()}")
    print("")
    print("Results:")
    if final_prep_pdb and final_prep_pdb.exists():
        print(f"  Preparation:    {final_prep_pdb}")
    if final_param_pdb and final_param_pdb.exists():
        print(f"  Parametrization: {final_param_pdb}")
    if run_orient and final_orient_pdb.exists():
        print(f"  Orientation:    {final_orient_pdb}")
    print("")


def main():
    """Main entry point for CLI."""
    if len(sys.argv) < 2:
        print("Usage: glyco-{prep,param,orient,all} [arguments]")
        print("\nAvailable commands:")
        print("  glyco-prep    - Run glycosylation preparation (Step 1)")
        print("  glyco-param   - Run parametrization (Step 2)")
        print("  glyco-orient  - Run carbohydrate orientation (Step 3)")
        print("  auto_glyco    - Run complete pipeline (Steps 1-3)")
        print("\nFor help on specific command: glyco-<command> -h")
        sys.exit(1)
    
    command = sys.argv[1]
    
    # Remove the command from argv
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    
    if command == "prep":
        run_glyco_prep()
    elif command == "param":
        run_glyco_param()
    elif command == "orient":
        run_glyco_orient()
    elif command == "all" or command == "auto_glyco":
        auto_glyco()
    else:
        print(f"Unknown command: {command}", file=sys.stderr)
        print("Available commands: prep, param, orient, all, auto_glyco", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
