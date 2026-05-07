#!/usr/bin/env python3
"""
Glycosylation Pipeline CLI - Unified interface for all steps
Usage: auto_glyco [arguments]
"""

import sys
import os
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).parent.absolute()
SCRIPTS_DIR = SCRIPT_DIR / "python_scripts"


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
        description='Run glycosylation preparation (Step 1)',
        usage='glyco-prep [options]'
    )
    parser.add_argument('--input-pdb', required=True, help='Input PDB file')
    parser.add_argument('--input-tsv', help='Input TSV file (Caselino format)')
    parser.add_argument('--input-glycosylator-tsv', help='Pre-processed glycosylator TSV')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--protein-residue-start', type=int, default=1, help='Protein residue start number (default: 1)')
    parser.add_argument('--rotate-atoms', default='OD1,CG,ND2,HD22,HD21,HB2,HB3', help='Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)')
    parser.add_argument('--fixed-atom', default='CB', help='Fixed atom (default: CB)')
    parser.add_argument('--center-atom', default='CA', help='Center atom (default: CA)')
    parser.add_argument('--radius', type=float, default=30.0, help='Radius for orientation (default: 30.0)')
    parser.add_argument('--rotation-step', type=float, default=1.0, help='Rotation step in degrees (default: 1)')
    
    args = parser.parse_args()
    
    # Create output directories
    output_dir = Path(args.output_dir)
    pdb_glycosylated_dir = output_dir / "PDB_PROTEIN_GLYCOSYLATED"
    pdb_glycosylated_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Run glycosylation script
    print("Step 1.1: Running glycosylation...")
    glycosylation_script = find_script('prep', '2-glycosylation_script.py')
    
    cmd = [
        sys.executable, glycosylation_script,
        args.input_pdb,
        args.input_tsv or "",
        args.input_glycosylator_tsv or "",
        str(pdb_glycosylated_dir)
    ]
    
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in glycosylation step", file=sys.stderr)
        sys.exit(1)
    
    # Step 2: Run ASN orientation
    print("Step 1.2: Running ASN orientation...")
    asn_orientation_script = find_script('prep', 'asn_orientation.py')
    
    # Find the glycosylated PDB
    glycosylated_pdb = pdb_glycosylated_dir / "protein_glycosylated.pdb"
    if not glycosylated_pdb.exists():
        # Try alternative naming
        pdb_files = list(pdb_glycosylated_dir.glob("*_glycosylated.pdb"))
        if pdb_files:
            glycosylated_pdb = pdb_files[0]
        else:
            print("Error: Could not find glycosylated PDB file", file=sys.stderr)
            sys.exit(1)
    
    cmd = [
        sys.executable, asn_orientation_script,
        str(glycosylated_pdb),
        str(output_dir),
        str(args.protein_residue_start),
        args.rotate_atoms,
        args.fixed_atom,
        args.center_atom,
        str(args.radius),
        str(args.rotation_step)
    ]
    
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in ASN orientation step", file=sys.stderr)
        sys.exit(1)
    
    # Find final output PDB
    final_pdb = pdb_glycosylated_dir / "protein_renumbered.pdb"
    if not final_pdb.exists():
        final_pdb = pdb_glycosylated_dir / "protein_asn_orientation.pdb"
    
    print(f"Preparation completed successfully!")
    print(f"Output PDB: {final_pdb}")
    
    return str(final_pdb)


def run_glyco_param():
    """Run parametrization (Step 2)"""
    parser = argparse.ArgumentParser(
        description='Run parametrization (Step 2)',
        usage='glyco-param [options]'
    )
    parser.add_argument('--prep-output-dir', required=True, help='Output directory from preparation step')
    parser.add_argument('--input-pdb', required=True, help='Input PDB file (from preparation step)')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--skip-charmm-download', action='store_true', help='Skip CHARMM force field download')
    
    args = parser.parse_args()
    
    # Create output directories
    output_dir = Path(args.output_dir)
    pdb_glycoprotein_dir = output_dir / "PDB_GLYCOPROTEIN"
    rtp_output_dir = output_dir / "RTP_FILES"
    pdb_glycoprotein_dir.mkdir(parents=True, exist_ok=True)
    rtp_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download CHARMM if needed
    if not args.skip_charmm_download:
        print("Downloading CHARMM36 force field...")
        charmm_script = find_script('param', '0-JSON_generator.py')
        # The script might include download logic
        # For now, we'll assume it's handled within the scripts
    
    # Step 1: Generate JSON from PDB
    print("Step 2.1: Generating JSON from PDB...")
    json_generator = find_script('param', '0-JSON_generator.py')
    
    cmd = [sys.executable, json_generator]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in JSON generation", file=sys.stderr)
        sys.exit(1)
    
    # Step 2: Parse PDB
    print("Step 2.2: Parsing PDB...")
    parser_pdb = find_script('param', '1-parser_pdb.py')
    
    cmd = [sys.executable, parser_pdb, args.input_pdb]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in PDB parsing", file=sys.stderr)
        sys.exit(1)
    
    # Step 3: Parse carbohydrate RTP
    print("Step 2.3: Parsing carbohydrate RTP...")
    parser_carb = find_script('param', '2-parser_carb_rtp.py')
    
    cmd = [sys.executable, parser_carb]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in carbohydrate RTP parsing", file=sys.stderr)
        sys.exit(1)
    
    # Step 4: Compare PDB and RTP
    print("Step 2.4: Comparing PDB and RTP...")
    compare_script = find_script('param', '3-comparison_pdb_rtp.py')
    
    cmd = [sys.executable, compare_script]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in comparison step", file=sys.stderr)
        sys.exit(1)
    
    # Step 5: Generate RTP
    print("Step 2.5: Generating RTP files...")
    rtp_gen1 = find_script('param', '4-rtp_generator_part1.py')
    rtp_gen2 = find_script('param', '4-rtp_generator_part2.py')
    rtp_gen3 = find_script('param', '4-rtp_generator_part3.py')
    
    for rtp_script in [rtp_gen1, rtp_gen2, rtp_gen3]:
        cmd = [sys.executable, rtp_script]
        result = subprocess.run(cmd, capture_output=False)
        if result.returncode != 0:
            print(f"Error in {Path(rtp_script).name}", file=sys.stderr)
            sys.exit(1)
    
    # Step 6: Acetylation replacement
    print("Step 2.6: Acetylation replacement...")
    acetylation_script = find_script('param', '5-acetylation_replacement.py')
    
    cmd = [sys.executable, acetylation_script]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in acetylation replacement", file=sys.stderr)
        sys.exit(1)
    
    # Step 7: Clean RTP
    print("Step 2.7: Cleaning RTP...")
    clean_rtp = find_script('param', '6-clean_rtp.py')
    
    cmd = [sys.executable, clean_rtp]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in RTP cleaning", file=sys.stderr)
        sys.exit(1)
    
    # Step 8: Run together scripts
    print("Step 2.8: Finalizing parametrization...")
    for part in range(1, 6):
        together_script = find_script('param', f'7-together_part_{part}.py')
        cmd = [sys.executable, together_script]
        result = subprocess.run(cmd, capture_output=False)
        if result.returncode != 0:
            print(f"Error in together_part_{part}", file=sys.stderr)
            sys.exit(1)
    
    # Step 9: Generate HDB
    print("Step 2.9: Generating HDB...")
    hdb_generator = find_script('param', '10-generation_hdb.py')
    
    cmd = [sys.executable, hdb_generator]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in HDB generation", file=sys.stderr)
        sys.exit(1)
    
    # Find final output PDB
    final_pdb = pdb_glycoprotein_dir / "protein_final_valence_corrected.pdb"
    
    print(f"Parametrization completed successfully!")
    print(f"Output PDB: {final_pdb}")
    
    return str(final_pdb)


def run_glyco_orient():
    """Run carbohydrate orientation (Step 3)"""
    parser = argparse.ArgumentParser(
        description='Run carbohydrate orientation (Step 3)',
        usage='glyco-orient [options]'
    )
    parser.add_argument('--input-pdb', required=True, help='Input PDB file (from parametrization step)')
    parser.add_argument('--param-output-dir', required=True, help='Parameter output directory')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--charmm-dir', help='CHARMM36 directory (optional)')
    parser.add_argument('--theta-step', type=int, default=10, help='Theta step for MCMC (default: 10)')
    parser.add_argument('--n-steps', type=int, default=10, help='Number of steps for MCMC (default: 10)')
    parser.add_argument('--max-cycles', type=int, default=5, help='Maximum cycles for MCMC (default: 5)')
    parser.add_argument('--radius', type=float, default=300.0, help='Radius for orientation (default: 300)')
    parser.add_argument('--use-coulomb', choices=['yes', 'no'], default='no', help='Use Coulomb potential (default: no)')
    parser.add_argument('--n-workers', type=int, default=1, help='Number of workers (default: 1)')
    parser.add_argument('--save-individual-glycans', action='store_true', help='Save individual glycans')
    parser.add_argument('--save-before-after', action='store_true', help='Save before/after structures')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Create output directories
    output_dir = Path(args.output_dir)
    optimized_dir = output_dir / "PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED"
    optimized_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Convert PDB to JSON
    print("Step 3.1: Converting PDB to JSON...")
    pdb_to_json = find_script('orient', '1-pdb_to_json.py')
    
    cmd = [sys.executable, pdb_to_json, args.input_pdb]
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in PDB to JSON conversion", file=sys.stderr)
        sys.exit(1)
    
    # Step 2: Run orientation optimization
    print("Step 3.2: Running MCMC orientation optimization...")
    optimize_mcmc = find_script('orient', '4-optimize_glycans_mcmc.py')
    
    cmd = [
        sys.executable, optimize_mcmc,
        "--input-pdb", args.input_pdb,
        "--output-dir", str(optimized_dir),
        "--theta-step", str(args.theta_step),
        "--n-steps", str(args.n_steps),
        "--max-cycles", str(args.max_cycles),
        "--radius", str(args.radius),
        "--use-coulomb", args.use_coulomb,
        "--n-workers", str(args.n_workers)
    ]
    
    if args.save_individual_glycans:
        cmd.append("--save-individual-glycans")
    if args.save_before_after:
        cmd.append("--save-before-after")
    if args.verbose:
        cmd.append("--verbose")
    if args.charmm_dir:
        cmd.extend(["--charmm-dir", args.charmm_dir])
    
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("Error in MCMC optimization", file=sys.stderr)
        sys.exit(1)
    
    # Find final output PDB
    final_pdb = optimized_dir / "protein_optimized.pdb"
    
    print(f"Orientation completed successfully!")
    print(f"Output PDB: {final_pdb}")
    
    return str(final_pdb)


def auto_glyco():
    """Run complete pipeline (Steps 1-3)"""
    parser = argparse.ArgumentParser(
        description='Complete automated glycosylation pipeline (Steps 1-3)',
        usage='auto_glyco [options]',
        epilog="""
Examples:
  auto_glyco --input-pdb protein.pdb --prep-output-dir ./prep --param-output-dir ./param --orient-output-dir ./orient
  
  auto_glyco --input-pdb protein.pdb --input-tsv table.tsv --prep-output-dir ./prep --param-output-dir ./param --orient-output-dir ./orient
  
  auto_glyco --prep-only --input-pdb protein.pdb --prep-output-dir ./prep
        """
    )
    
    # Required arguments
    parser.add_argument('--input-pdb', required=True, help='Input PDB file')
    parser.add_argument('--prep-output-dir', required=True, help='Output directory for preparation step')
    parser.add_argument('--param-output-dir', required=True, help='Output directory for parametrization step')
    parser.add_argument('--orient-output-dir', required=True, help='Output directory for orientation step')
    
    # Input options
    parser.add_argument('--input-tsv', help='Input TSV file (Caselino format)')
    parser.add_argument('--input-glycosylator-tsv', help='Pre-processed glycosylator TSV')
    
    # Preparation options
    parser.add_argument('--protein-residue-start', type=int, default=1, help='Protein residue start number (default: 1)')
    parser.add_argument('--rotate-atoms', default='OD1,CG,ND2,HD22,HD21,HB2,HB3', help='Atoms to rotate')
    parser.add_argument('--fixed-atom', default='CB', help='Fixed atom (default: CB)')
    parser.add_argument('--center-atom', default='CA', help='Center atom (default: CA)')
    parser.add_argument('--radius-prep', type=float, default=30.0, help='Radius for orientation (default: 30.0)')
    parser.add_argument('--rotation-step', type=float, default=1.0, help='Rotation step in degrees (default: 1)')
    
    # Parametrization options
    parser.add_argument('--skip-charmm-download', action='store_true', help='Skip CHARMM force field download')
    
    # Orientation options
    parser.add_argument('--charmm-dir', help='CHARMM36 directory (optional)')
    parser.add_argument('--theta-step', type=int, default=10, help='Theta step for MCMC (default: 10)')
    parser.add_argument('--n-steps', type=int, default=10, help='Number of steps for MCMC (default: 10)')
    parser.add_argument('--max-cycles', type=int, default=5, help='Maximum cycles for MCMC (default: 5)')
    parser.add_argument('--radius-orient', type=float, default=300.0, help='Radius for orientation (default: 300)')
    parser.add_argument('--use-coulomb', choices=['yes', 'no'], default='no', help='Use Coulomb potential (default: no)')
    parser.add_argument('--n-workers', type=int, default=1, help='Number of workers (default: 1)')
    parser.add_argument('--save-individual-glycans', action='store_true', help='Save individual glycans')
    parser.add_argument('--save-before-after', action='store_true', help='Save before/after structures')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    # Step selection
    parser.add_argument('--prep-only', action='store_true', help='Run only preparation step')
    parser.add_argument('--param-only', action='store_true', help='Run only parametrization step')
    parser.add_argument('--orient-only', action='store_true', help='Run only orientation step')
    
    args = parser.parse_args()
    
    # Create output directories
    Path(args.prep_output_dir).mkdir(parents=True, exist_ok=True)
    Path(args.param_output_dir).mkdir(parents=True, exist_ok=True)
    Path(args.orient_output_dir).mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*60)
    print("AUTOMATED GLYCOSYLATION PIPELINE")
    print("="*60)
    print(f"Start time: {datetime.now()}")
    print(f"Input PDB: {args.input_pdb}")
    print("="*60)
    
    final_prep = None
    final_param = None
    final_orient = None
    
    try:
        # Track intermediate PDBs
        current_pdb = args.input_pdb
        
        # Step 1: Preparation
        if not args.param_only and not args.orient_only:
            print("\n" + "="*60)
            print("STEP 1: GLYCOSYLATION PREPARATION")
            print("="*60)
            
            # Create a modified args for prep
            class PrepArgs:
                pass
            
            prep_args = PrepArgs()
            prep_args.input_pdb = args.input_pdb
            prep_args.input_tsv = args.input_tsv
            prep_args.input_glycosylator_tsv = args.input_glycosylator_tsv
            prep_args.output_dir = args.prep_output_dir
            prep_args.protein_residue_start = args.protein_residue_start
            prep_args.rotate_atoms = args.rotate_atoms
            prep_args.fixed_atom = args.fixed_atom
            prep_args.center_atom = args.center_atom
            prep_args.radius = args.radius_prep
            prep_args.rotation_step = args.rotation_step
            
            # Save original argv
            original_argv = sys.argv
            sys.argv = ['glyco-prep']
            for key, value in vars(prep_args).items():
                if value is not None and not isinstance(value, argparse.Namespace):
                    sys.argv.append(f'--{key.replace("_", "-")}')
                    sys.argv.append(str(value))
            
            final_prep = run_glyco_prep()
            current_pdb = final_prep
            
            # Restore original argv
            sys.argv = original_argv
        
        # Step 2: Parametrization
        if not args.prep_only and not args.orient_only:
            print("\n" + "="*60)
            print("STEP 2: PARAMETRIZATION")
            print("="*60)
            
            # Create a modified args for param
            class ParamArgs:
                pass
            
            param_args = ParamArgs()
            param_args.prep_output_dir = args.prep_output_dir
            param_args.input_pdb = current_pdb
            param_args.output_dir = args.param_output_dir
            param_args.skip_charmm_download = args.skip_charmm_download
            
            # Save original argv
            original_argv = sys.argv
            sys.argv = ['glyco-param']
            for key, value in vars(param_args).items():
                if value is not None and not isinstance(value, argparse.Namespace):
                    if key == 'skip_charmm_download' and value:
                        sys.argv.append(f'--skip-charmm-download')
                    elif value is not None:
                        sys.argv.append(f'--{key.replace("_", "-")}')
                        sys.argv.append(str(value))
            
            final_param = run_glyco_param()
            current_pdb = final_param
            
            # Restore original argv
            sys.argv = original_argv
        
        # Step 3: Orientation
        if not args.prep_only and not args.param_only:
            print("\n" + "="*60)
            print("STEP 3: CARBOHYDRATE ORIENTATION")
            print("="*60)
            
            # Create a modified args for orient
            class OrientArgs:
                pass
            
            orient_args = OrientArgs()
            orient_args.input_pdb = current_pdb
            orient_args.param_output_dir = args.param_output_dir
            orient_args.output_dir = args.orient_output_dir
            orient_args.charmm_dir = args.charmm_dir
            orient_args.theta_step = args.theta_step
            orient_args.n_steps = args.n_steps
            orient_args.max_cycles = args.max_cycles
            orient_args.radius = args.radius_orient
            orient_args.use_coulomb = args.use_coulomb
            orient_args.n_workers = args.n_workers
            orient_args.save_individual_glycans = args.save_individual_glycans
            orient_args.save_before_after = args.save_before_after
            orient_args.verbose = args.verbose
            
            # Save original argv
            original_argv = sys.argv
            sys.argv = ['glyco-orient']
            for key, value in vars(orient_args).items():
                if value is not None and not isinstance(value, argparse.Namespace):
                    if isinstance(value, bool) and value:
                        sys.argv.append(f'--{key.replace("_", "-")}')
                    elif not isinstance(value, bool):
                        sys.argv.append(f'--{key.replace("_", "-")}')
                        sys.argv.append(str(value))
            
            final_orient = run_glyco_orient()
            
            # Restore original argv
            sys.argv = original_argv
        
        # Summary
        print("\n" + "="*60)
        print("PIPELINE COMPLETED SUCCESSFULLY!")
        print("="*60)
        print(f"End time: {datetime.now()}")
        print("\nResults:")
        if final_prep:
            print(f"  Preparation:       {final_prep}")
        if final_param:
            print(f"  Parametrization:   {final_param}")
        if final_orient:
            print(f"  Orientation:       {final_orient}")
        print("="*60)
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main entry point for CLI."""
    if len(sys.argv) < 2:
        print("""Usage: glyco-{prep,param,orient,all} [arguments]
        # Run the complete pipeline
    auto_glyco \
        --input-pdb protein.pdb \
        --input-tsv table.tsv \
        --prep-output-dir ./prep \
        --param-output-dir ./param \
        --orient-output-dir ./orient

    # Or using the alternative name
    glyco-all \
        --input-pdb protein.pdb \
        --input-tsv table.tsv \
        --prep-output-dir ./prep \
        --param-output-dir ./param \
        --orient-output-dir ./orient

    # Run only preparation step
        auto_glyco \
        --prep-only \
        --input-pdb protein.pdb \
        --prep-output-dir ./prep \
        --param-output-dir ./param \
        --orient-output-dir ./orient

    # Run with custom parameters
    auto_glyco \
        --input-pdb protein.pdb \
        --prep-output-dir ./prep \
        --param-output-dir ./param \
        --orient-output-dir ./orient \
        --protein-residue-start 10 \
        --theta-step 20 \
        --n-steps 50 \
        --verbose   
        """)
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
