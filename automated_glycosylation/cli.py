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


def run_glyco_prep(args=None):
    """Run preparation step"""
    if args is None:
        parser = argparse.ArgumentParser(
            prog='glyco-prep',
            description='Run glycosylation preparation step',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        
        # Required arguments
        parser.add_argument('--input-pdb', required=True, help='Input PDB file')
        parser.add_argument('--output-dir', required=True, help='Output directory')
        
        # Optional input files
        parser.add_argument('--input-tsv', help='Input TSV file (glycosylation table)')
        parser.add_argument('--input-glycosylator-tsv', help='Pre-processed glycosylator TSV')
        
        # Optional parameters
        parser.add_argument('--protein-residue-start', type=int, default=1, 
                           help='Protein residue start number (default: 1)')
        parser.add_argument('--rotate-atoms', default="OD1,CG,ND2,HD22,HD21,HB2,HB3", 
                           help='Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)')
        parser.add_argument('--fixed-atom', default="CB", 
                           help='Fixed atom (default: CB)')
        parser.add_argument('--center-atom', default="CA", 
                           help='Center atom (default: CA)')
        parser.add_argument('--radius', type=float, default=30.0, 
                           help='Radius for orientation (default: 30.0)')
        parser.add_argument('--rotation-step', type=int, default=1, 
                           help='Rotation step in degrees (default: 1)')
        
        args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb_abs = get_abs_path(args.input_pdb)
    output_dir_abs = get_abs_path(args.output_dir)
    input_tsv_abs = get_abs_path(args.input_tsv) if args.input_tsv else None
    input_glycosylator_tsv_abs = get_abs_path(args.input_glycosylator_tsv) if args.input_glycosylator_tsv else None
    
    print("=" * 60)
    print("STEP 1: Glycosylation Preparation")
    print("=" * 60)
    print(f"Input PDB: {input_pdb_abs}")
    print(f"Output dir: {output_dir_abs}")
    if input_tsv_abs:
        print(f"Input TSV: {input_tsv_abs}")
    
    # Get base name for generic file naming
    base_name = get_base_filename(input_tsv_abs) if input_tsv_abs else "glycosylation"
    
    # Create directories
    os.makedirs(output_dir_abs, exist_ok=True)
    os.makedirs(f"{output_dir_abs}/TSV", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/PDB_PROTEIN_GLYCOSYLATED", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/EXTRACTED_CARBOHYDRATES", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/TO_TOP", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/TO_TOP/PDB", exist_ok=True)
    
    # Step 0: Correct table (if input TSV provided)
    input_tsv_glycosylator = None
    if input_tsv_abs:
        output_tsv = f"{output_dir_abs}/TSV/{base_name}_corrected.tsv"
        output_results = f"{output_dir_abs}/EXTRACTED_CARBOHYDRATES"
        
        print("\nStep 0: Correcting input table...")
        script_0 = find_script('prep', '0-correcting_caselino_table_for_variants.py')
        cmd = [
            sys.executable, script_0,
            "--input", input_tsv_abs,
            "--output", output_tsv
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        # Step 1: Convert to IUPAC notation
        print("\nStep 1: Converting to IUPAC notation...")
        output_glycosylator = f"{output_dir_abs}/TSV/{base_name}_glycosylator.tsv"
        script_1 = find_script('prep', '1-iupac_converted.py')
        cmd = [
            sys.executable, script_1,
            "--input_tsv", output_tsv,
            "--output_tsv", output_glycosylator,
            "--output_dir", output_results
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        input_tsv_glycosylator = output_glycosylator
    elif input_glycosylator_tsv_abs:
        input_tsv_glycosylator = input_glycosylator_tsv_abs
    else:
        print("WARNING: No input TSV or glycosylator TSV provided. Skipping table processing.")
    
    # Step 2a: Optimize asparagine orientations
    print("\nStep 2a: Optimizing asparagine orientations...")
    pdb_asn_optimized = f"{output_dir_abs}/PDB_PROTEIN_GLYCOSYLATED/protein_asn_orientation.pdb"
    script_asn = find_script('prep', 'asn_orientation.py')
    
    cmd = [
        sys.executable, script_asn,
        input_pdb_abs,
        "--rotate-atoms", args.rotate_atoms,
        "--fixed-atom", args.fixed_atom,
        "--center-atom", args.center_atom,
        "--radius", str(args.radius),
        "--rotation-step", str(args.rotation_step),
        "-o", pdb_asn_optimized
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # Step 2b: Run glycosylation (only if we have glycosylator data)
    final_pdb = pdb_asn_optimized
    if input_tsv_glycosylator:
        print("\nStep 2b: Running glycosylation...")
        glycosylated_pdb = f"{output_dir_abs}/PDB_PROTEIN_GLYCOSYLATED/protein_glycosylated.pdb"
        script_glyco = find_script('prep', '2-glycosylation_script.py')
        cmd = [
            sys.executable, script_glyco,
            "--input_tsv_glycosylator", input_tsv_glycosylator,
            "--input_pdb_protein", pdb_asn_optimized,
            "--protein_residue_start", str(args.protein_residue_start),
            "--output", glycosylated_pdb
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        # Step 3: Correct chain labels and residue numbers
        print("\nStep 3: Correcting chain labels and residue numbers...")
        renumbered_pdb = f"{output_dir_abs}/PDB_PROTEIN_GLYCOSYLATED/protein_renumbered.pdb"
        script_correct = find_script('prep', '3-correction_chain_labels_residue_numbers.py')
        cmd = [
            sys.executable, script_correct,
            glycosylated_pdb,
            renumbered_pdb
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        # Step 4: Extract glycans coordinates
        print("\nStep 4: Extracting glycans coordinates...")
        output_noh = f"{output_dir_abs}/PDB_PROTEIN_GLYCOSYLATED/protein_without_H.pdb"
        script_extract = find_script('prep', '4-extract_coordinates_of_glycans_from_structure.py')
        cmd = [
            sys.executable, script_extract,
            "--input_pdb", renumbered_pdb,
            "--output_noH", output_noh,
            "--output_dir", f"{output_dir_abs}/TO_TOP"
        ]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        final_pdb = renumbered_pdb
    else:
        print("\nSkipping glycosylation steps (no glycosylator data provided)")
    
    print("\n" + "=" * 60)
    print(f"Preparation completed successfully!")
    print(f"Output directory: {output_dir_abs}")
    print(f"Final PDB: {final_pdb}")
    print("=" * 60)
    
    return final_pdb


def run_glyco_param(args=None):
    """Run parametrization step"""
    if args is None:
        parser = argparse.ArgumentParser(
            prog='glyco-param',
            description='Run parametrization step',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        
        parser.add_argument('--prep-output-dir', required=True, 
                           help='Output directory from prep step')
        parser.add_argument('--input-pdb', required=True, 
                           help='Input PDB file from prep step')
        parser.add_argument('--output-dir', required=True, 
                           help='Output directory for parametrization')
        parser.add_argument('--skip-charmm-download', action='store_true',
                           help='Skip CHARMM force field download (use existing)')
        
        args = parser.parse_args()
    
    # Convert to absolute paths
    prep_output_dir_abs = get_abs_path(args.prep_output_dir)
    input_pdb_abs = get_abs_path(args.input_pdb)
    output_dir_abs = get_abs_path(args.output_dir)
    
    print("=" * 60)
    print("STEP 2: Parametrization")
    print("=" * 60)
    print(f"Prep output dir: {prep_output_dir_abs}")
    print(f"Input PDB: {input_pdb_abs}")
    print(f"Output dir: {output_dir_abs}")
    
    # Create directories
    os.makedirs(output_dir_abs, exist_ok=True)
    os.makedirs(f"{output_dir_abs}/JSON", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/PDB_GLYCOPROTEIN", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/VALENCE_GLYCAN_VARIANTS", exist_ok=True)
    
    # Setup CHARMM force field
    charmm_dir = f"{output_dir_abs}/charmm36.ff"
    charmm_rtp = f"{charmm_dir}/carb.rtp"
    charmm_hdb = f"{charmm_dir}/carb.hdb"
    charmm_rtp_backup = f"{charmm_dir}/carb.rtp.backup"
    charmm_hdb_backup = f"{charmm_dir}/carb.hdb.backup"
    
    # Create CHARMM directory if it doesn't exist
    os.makedirs(charmm_dir, exist_ok=True)
    
    # Check if carb.rtp exists, create empty if not
    if not os.path.exists(charmm_rtp):
        print(f"Creating empty {charmm_rtp}")
        with open(charmm_rtp, 'w') as f:
            f.write("; Auto-generated empty RTP file\n")
        with open(charmm_hdb, 'w') as f:
            f.write("; Auto-generated empty HDB file\n")
    
    # Download CHARMM if needed
    if not args.skip_charmm_download:
        if not os.path.exists(charmm_rtp_backup) or not os.path.exists(charmm_hdb_backup):
            print("\nDownloading CHARMM36 force field...")
            subprocess.run(["wget", "-O", "charmm36.tgz", "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz"], check=True)
            subprocess.run(["tar", "-xzf", "charmm36.tgz"], check=True)
            subprocess.run(["mv", "charmm36-jul2022.ff", charmm_dir], check=True)
            subprocess.run(["rm", "charmm36.tgz"], check=True)
            
            # Create backups
            if os.path.exists(charmm_rtp):
                subprocess.run(["cp", charmm_rtp, charmm_rtp_backup], check=True)
            if os.path.exists(charmm_hdb):
                subprocess.run(["cp", charmm_hdb, charmm_hdb_backup], check=True)
        else:
            print("\nRestoring CHARMM force field from backups...")
            subprocess.run(["cp", charmm_rtp_backup, charmm_rtp], check=True)
            subprocess.run(["cp", charmm_hdb_backup, charmm_hdb], check=True)
    
    # Generate JSONs
    print("\nGenerating JSONs for glycans...")
    script_json = find_script('param', '0-JSON_generator.py')
    cmd = [
        sys.executable, script_json,
        "--base_dir", prep_output_dir_abs,
        "--output_dir", f"{output_dir_abs}/JSON"
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # Process each glycan directory
    json_dir = f"{output_dir_abs}/JSON"
    if not Path(json_dir).exists() or not any(Path(json_dir).iterdir()):
        print(f"WARNING: No glycan directories found in {json_dir}")
        print("Skipping parametrization processing")
        return None
    
    # Process each glycan
    for glycan_dir in Path(json_dir).iterdir():
        if not glycan_dir.is_dir():
            continue
        
        basename = glycan_dir.name
        print(f"\n{'='*60}")
        print(f"Processing: {basename}")
        print(f"{'='*60}")
        
        pdb_file = glycan_dir / f"{basename}.pdb"
        json_file = glycan_dir / f"{basename}.json"
        parser_file = glycan_dir / f"{basename}_parser.pkl"
        rtp_pickle = glycan_dir / "carb_residues.pkl"
        rtp_modified = glycan_dir / "carb_modified.rtp"
        rtp_unique = glycan_dir / "carb_unique.rtp"
        
        # Check if required files exist
        if not pdb_file.exists():
            print(f"  ERROR: {pdb_file} not found, skipping")
            continue
        
        # Verify PDB has atoms
        try:
            with open(pdb_file, 'r') as f:
                content = f.read()
                if "HETATM" not in content and "ATOM" not in content:
                    print(f"  ERROR: No atoms found in {pdb_file}, skipping")
                    continue
            print(f"  PDB file OK: {pdb_file} ({pdb_file.stat().st_size} bytes)")
        except Exception as e:
            print(f"  ERROR reading PDB: {e}, skipping")
            continue
        
        # Change to glycan directory for all operations
        original_cwd = os.getcwd()
        os.chdir(str(glycan_dir))
        
        try:
            # 1. Parse PDB
            print(f"  Step 1: Parsing PDB...")
            script_1 = find_script('param', '1-parser_pdb.py')
            cmd = [sys.executable, script_1, str(pdb_file.name), "-o", str(parser_file.name)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  ERROR in parser_pdb: {result.stderr[:200]}")
                continue
            print(f"    OK: Created {parser_file.name}")
            
            # 2. Parse CHARMM RTP
            print(f"  Step 2: Parsing CHARMM RTP...")
            script_2 = find_script('param', '2-parser_carb_rtp.py')
            cmd = [sys.executable, script_2, charmm_rtp, "-o", str(rtp_pickle.name)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  ERROR in parser_carb_rtp: {result.stderr[:200]}")
                continue
            print(f"    OK: Created {rtp_pickle.name}")
            
            # 3. Compare PDB and RTP
            print(f"  Step 3: Comparing PDB and RTP...")
            script_3 = find_script('param', '3-comparison_pdb_rtp.py')
            cmd = [sys.executable, script_3, "--pdb", str(parser_file.name), "--rtp", str(rtp_pickle.name)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"    OK: Comparison complete")
            
            # 4. Generate RTP part 1
            print(f"  Step 4: Generating RTP part 1...")
            script_4a = find_script('param', '4-rtp_generator_part1.py')
            cmd = [sys.executable, script_4a, "-p", str(parser_file.name), "-r", str(rtp_pickle.name)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"    OK: RTP part 1 complete")
            
            # 5. Generate RTP part 2
            print(f"  Step 5: Generating RTP part 2...")
            script_4b = find_script('param', '4-rtp_generator_part2.py')
            cmd = [sys.executable, script_4b, "--pdb", str(parser_file.name), "--rtp", str(rtp_pickle.name), "--json", str(json_file.name)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"    OK: RTP part 2 complete")
            
            # 6. Generate RTP part 3
            print(f"  Step 6: Generating RTP part 3...")
            script_4c = find_script('param', '4-rtp_generator_part3.py')
            cmd = [sys.executable, script_4c, "--pdb", str(parser_file.name), "--rtp", str(rtp_pickle.name), "--json", str(json_file.name), "--output", str(rtp_modified.name)]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"    OK: Created {rtp_modified.name}")
            
            # 7. Acetylation replacement
            modified_pdb = glycan_dir / f"{basename}_modified.pdb"
            if modified_pdb.exists():
                print(f"  Step 7: Acetylation replacement...")
                script_5 = find_script('param', '5-acetylation_replacement.py')
                cmd = [sys.executable, script_5, str(pdb_file.name), str(modified_pdb.name)]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"    OK: Acetylation complete")
            
            # 8. Clean RTP
            if rtp_modified.exists():
                print(f"  Step 8: Cleaning RTP...")
                dir_letter = basename[0] if basename else "X"
                script_6 = find_script('param', '6-clean_rtp.py')
                cmd = [sys.executable, script_6, str(rtp_modified.name), str(rtp_unique.name), dir_letter]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"    OK: Created {rtp_unique.name}")
            
            print(f"  ✓ Successfully processed {basename}")
            
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Error processing {basename}: {e}")
            if e.stderr:
                print(f"  stderr: {e.stderr[:200]}")
        except Exception as e:
            print(f"  ✗ Unexpected error: {e}")
        finally:
            # Change back to original directory
            os.chdir(original_cwd)
    
    # Unify RTP/HDB files
    print("\n" + "=" * 60)
    print("Unifying RTP and HDB files...")
    print("=" * 60)
    
    carb_total = f"{json_dir}/carb_unique_total.rtp"
    carb_redundance = f"{json_dir}/carb_redundance_removed.rtp"
    carb_hdb = f"{json_dir}/carb_redundance_removed.hdb"
    
    scripts_unify = [
        ("7-together_part_1.py", ["--input", json_dir, "--output", carb_total]),
        ("7-together_part_2.py", ["--input", carb_total, "--output", carb_redundance]),
        ("7-together_part_3.py", ["--input", carb_redundance, "--output", charmm_rtp]),
        ("7-together_part_4.py", ["--input", carb_redundance]),
        ("7-together_part_5.py", ["--input", carb_redundance, "--output", carb_hdb]),
    ]
    
    for script_name, script_args in scripts_unify:
        try:
            script_path = find_script('param', script_name)
            cmd = [sys.executable, script_path] + script_args
            subprocess.run(cmd, check=True)
            print(f"  ✓ {script_name} completed")
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Error in {script_name}: {e}")
            continue
    
    # Build final glycoprotein
    final_result = None
    if Path(carb_redundance).exists():
        print("\n" + "=" * 60)
        print("Building final glycoprotein...")
        print("=" * 60)
        
        input_protein = input_pdb_abs
        output_protein = f"{output_dir_abs}/PDB_GLYCOPROTEIN/protein_corrected.pdb"
        final_structure = f"{output_dir_abs}/PDB_GLYCOPROTEIN/protein_final_connected.pdb"
        final_valence = f"{output_dir_abs}/PDB_GLYCOPROTEIN/protein_final_valence_corrected.pdb"
        
        try:
            print("  Step 1: Generating glycoprotein...")
            script_glycoprotein = find_script('param', '8-glycoprotein.py')
            cmd = [
                sys.executable, script_glycoprotein,
                "--protein", input_protein,
                "--carbs_dir", json_dir,
                "--output", output_protein,
                "--keep_hydrogens_carb",
                "--keep_hydrogens_prot"
            ]
            subprocess.run(cmd, check=True)
            print(f"    OK: Created {output_protein}")
            
            print("  Step 2: Connecting glycosylation...")
            script_connect = find_script('param', '9-conection_glycosilation_without_TER.py')
            cmd = [
                sys.executable, script_connect,
                "--glycosylated", output_protein,
                "--conect", input_protein,
                "--output", final_structure
            ]
            subprocess.run(cmd, check=True)
            print(f"    OK: Created {final_structure}")
            
            print("  Step 3: Identifying glycosylation...")
            script_identify = find_script('param', 'glycosylation_identifying.py')
            cmd = [sys.executable, script_identify, final_structure, final_valence]
            subprocess.run(cmd, check=True)
            print(f"    OK: Created {final_valence}")
            
            # Generate HDB
            print("  Step 4: Generating HDB...")
            script_hdb = find_script('param', '10-generation_hdb.py')
            cmd = [
                sys.executable, script_hdb,
                carb_redundance,
                charmm_hdb,
                "-o", f"{json_dir}/carb_modified.hdb"
            ]
            subprocess.run(cmd, check=True)
            print(f"    OK: HDB generation complete")
            
            # Generate variants
            print("  Step 5: Generating variants...")
            script_variants = find_script('param', 'glycosylation_variants.py')
            cmd = [
                sys.executable, script_variants,
                "-p", final_valence,
                "-r", carb_redundance,
                "-d", carb_hdb,
                "-o", f"{output_dir_abs}/VALENCE_GLYCAN_VARIANTS"
            ]
            subprocess.run(cmd, check=True)
            print(f"    OK: Variants generated")
            
            # Update CHARMM files with variants
            variant_rtp = f"{output_dir_abs}/VALENCE_GLYCAN_VARIANTS/protein_final_valence_corrected_variants.rtp"
            variant_hdb = f"{output_dir_abs}/VALENCE_GLYCAN_VARIANTS/protein_final_valence_corrected_variants.hdb"
            
            if os.path.exists(variant_rtp) and os.path.exists(variant_hdb):
                print("  Step 6: Updating CHARMM files...")
                with open(charmm_rtp, 'a') as f:
                    with open(variant_rtp, 'r') as v:
                        f.write(v.read())
                with open(charmm_hdb, 'a') as f:
                    with open(variant_hdb, 'r') as v:
                        f.write(v.read())
                
                # Update backups
                subprocess.run(["cp", charmm_rtp, charmm_rtp_backup], check=True)
                subprocess.run(["cp", charmm_hdb, charmm_hdb_backup], check=True)
                print(f"    OK: CHARMM files updated")
            
            final_result = final_valence
        except Exception as e:
            print(f"  ✗ Error in glycoprotein building: {e}")
            final_result = None
    else:
        print(f"\nWARNING: {carb_redundance} not found, skipping final building")
    
    print("\n" + "=" * 60)
    print(f"Parametrization completed successfully!")
    print(f"Output directory: {output_dir_abs}")
    if final_result:
        print(f"Final PDB: {final_result}")
    print("=" * 60)
    
    return final_result


def run_glyco_orient(args=None):
    """Run carbohydrate orientation step"""
    if args is None:
        parser = argparse.ArgumentParser(
            prog='glyco-orient',
            description='Run carbohydrate orientation step',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        
        parser.add_argument('--input-pdb', required=True, 
                           help='Input PDB file from param step')
        parser.add_argument('--param-output-dir', required=True, 
                           help='Output directory from param step')
        parser.add_argument('--output-dir', required=True, 
                           help='Output directory for orientation')
        parser.add_argument('--charmm-dir', 
                           help='CHARMM36 directory (default: from param output)')
        parser.add_argument('--theta-step', type=int, default=10, 
                           help='Theta step for MCMC (default: 10)')
        parser.add_argument('--n-steps', type=int, default=10, 
                           help='Number of steps for MCMC (default: 10)')
        parser.add_argument('--max-cycles', type=int, default=5, 
                           help='Maximum cycles for MCMC (default: 5)')
        parser.add_argument('--radius', type=float, default=300, 
                           help='Radius for orientation (default: 300)')
        parser.add_argument('--use-coulomb', default='no', choices=['yes', 'no'], 
                           help='Use Coulomb potential (default: no)')
        parser.add_argument('--n-workers', type=int, default=1, 
                           help='Number of workers for parallel processing (default: 1)')
        parser.add_argument('--save-individual-glycans', action='store_true', 
                           help='Save individual glycans')
        parser.add_argument('--save-before-after', action='store_true', 
                           help='Save before/after structures')
        parser.add_argument('--verbose', action='store_true', 
                           help='Verbose output')
        
        args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb_abs = get_abs_path(args.input_pdb)
    param_output_dir_abs = get_abs_path(args.param_output_dir)
    output_dir_abs = get_abs_path(args.output_dir)
    charmm_dir_abs = get_abs_path(args.charmm_dir) if args.charmm_dir else None
    
    print("=" * 60)
    print("STEP 3: Carbohydrate Orientation")
    print("=" * 60)
    print(f"Input PDB: {input_pdb_abs}")
    print(f"Param output dir: {param_output_dir_abs}")
    print(f"Output dir: {output_dir_abs}")
    
    # Create directories
    os.makedirs(output_dir_abs, exist_ok=True)
    os.makedirs(f"{output_dir_abs}/JSON_FILES", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED", exist_ok=True)
    os.makedirs(f"{output_dir_abs}/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PDB_CARB_ONLY", exist_ok=True)
    
    # Convert PDB to JSON
    print("\nStep 1: Converting PDB to JSON...")
    script_pdb_json = find_script('orient', '1-pdb_to_json.py')
    pdb_json = f"{output_dir_abs}/JSON_FILES/pdb_to_json.json"
    cmd = [
        sys.executable, script_pdb_json,
        "--input_pdb", input_pdb_abs,
        "--output_json", pdb_json
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"  OK: Created {pdb_json}")
    
    # Add CHARMM36 parameters
    print("\nStep 2: Adding CHARMM36 parameters...")
    script_charmm = find_script('orient', '3-adding_chamm36_parameters.py')
    charmm_dir = charmm_dir_abs or f"{param_output_dir_abs}/charmm36.ff"
    glycan_data_json = f"{output_dir_abs}/JSON_FILES/glycan_data_charmm36.json"
    cmd = [
        sys.executable, script_charmm,
        "--input_json", pdb_json,
        "--charmm_dir", charmm_dir,
        "--output_json", glycan_data_json
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"  OK: Created {glycan_data_json}")
    
    # Optimize glycans
    print("\nStep 3: Optimizing glycans using MCMC...")
    script_mcmc = find_script('orient', '4-optimize_glycans_mcmc.py')
    output_pdb = f"{output_dir_abs}/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/protein_optimized.pdb"
    output_json = f"{output_dir_abs}/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/glycan_optimized.json"
    report_file = f"{output_dir_abs}/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/report.txt"
    
    cmd = [
        sys.executable, script_mcmc,
        "--input_json", glycan_data_json,
        "--output_json", output_json,
        "--output_pdb", output_pdb,
        "--glycans_output_dir", f"{output_dir_abs}/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PDB_CARB_ONLY",
        "--theta_step", str(args.theta_step),
        "--n_steps", str(args.n_steps),
        "--max_cycles", str(args.max_cycles),
        "--radius", str(args.radius),
        "--use_coulomb", args.use_coulomb,
        "--n_workers", str(args.n_workers),
        "--report_file", report_file
    ]
    
    if args.save_individual_glycans:
        cmd.append("--save_individual_glycans")
    if args.save_before_after:
        cmd.append("--save_before_after")
    if args.verbose:
        cmd.append("--verbose")
    
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    print("\n" + "=" * 60)
    print(f"Orientation completed successfully!")
    print(f"Output directory: {output_dir_abs}")
    print(f"Optimized PDB: {output_pdb}")
    print("=" * 60)
    
    return output_pdb


def run_glyco_all(args=None):
    """Run complete pipeline"""
    if args is None:
        parser = argparse.ArgumentParser(
            prog='glyco-all',
            description='Run complete glycosylation pipeline (Steps 1-3)',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        
        # Required arguments
        parser.add_argument('--input-pdb', required=True, help='Input PDB file')
        parser.add_argument('--prep-output-dir', required=True, 
                           help='Output directory for prep step')
        parser.add_argument('--param-output-dir', required=True, 
                           help='Output directory for param step')
        parser.add_argument('--orient-output-dir', required=True, 
                           help='Output directory for orient step')
        
        # Optional input files
        parser.add_argument('--input-tsv', help='Input TSV file (glycosylation table)')
        parser.add_argument('--input-glycosylator-tsv', help='Pre-processed glycosylator TSV')
        
        # Prep parameters
        parser.add_argument('--protein-residue-start', type=int, default=1, 
                           help='Protein residue start number (default: 1)')
        parser.add_argument('--rotate-atoms', default="OD1,CG,ND2,HD22,HD21,HB2,HB3", 
                           help='Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)')
        parser.add_argument('--fixed-atom', default="CB", 
                           help='Fixed atom (default: CB)')
        parser.add_argument('--center-atom', default="CA", 
                           help='Center atom (default: CA)')
        parser.add_argument('--radius-prep', type=float, default=30.0, 
                           help='Radius for prep orientation (default: 30.0)')
        parser.add_argument('--rotation-step', type=int, default=1, 
                           help='Rotation step in degrees (default: 1)')
        
        # Param parameters
        parser.add_argument('--skip-charmm-download', action='store_true',
                           help='Skip CHARMM force field download')
        
        # Orient parameters
        parser.add_argument('--charmm-dir', help='CHARMM36 directory (optional)')
        parser.add_argument('--theta-step', type=int, default=10, 
                           help='Theta step for MCMC (default: 10)')
        parser.add_argument('--n-steps', type=int, default=10, 
                           help='Number of steps for MCMC (default: 10)')
        parser.add_argument('--max-cycles', type=int, default=5, 
                           help='Maximum cycles for MCMC (default: 5)')
        parser.add_argument('--radius-orient', type=float, default=300, 
                           help='Radius for orientation (default: 300)')
        parser.add_argument('--use-coulomb', default='no', choices=['yes', 'no'], 
                           help='Use Coulomb potential (default: no)')
        parser.add_argument('--n-workers', type=int, default=1, 
                           help='Number of workers (default: 1)')
        parser.add_argument('--save-individual-glycans', action='store_true', 
                           help='Save individual glycans')
        parser.add_argument('--save-before-after', action='store_true', 
                           help='Save before/after structures')
        parser.add_argument('--verbose', action='store_true', 
                           help='Verbose output')
        
        args = parser.parse_args()
    
    print("\n" + "=" * 60)
    print("RUNNING COMPLETE GLYCOSYLATION PIPELINE")
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60 + "\n")
    
    # Create a simple namespace for prep
    class PrepArgs:
        pass
    
    prep_args = PrepArgs()
    prep_args.input_pdb = args.input_pdb
    prep_args.input_tsv = getattr(args, 'input_tsv', None)
    prep_args.input_glycosylator_tsv = getattr(args, 'input_glycosylator_tsv', None)
    prep_args.output_dir = args.prep_output_dir
    prep_args.protein_residue_start = args.protein_residue_start
    prep_args.rotate_atoms = args.rotate_atoms
    prep_args.fixed_atom = args.fixed_atom
    prep_args.center_atom = args.center_atom
    prep_args.radius = args.radius_prep
    prep_args.rotation_step = args.rotation_step
    
    # Step 1: Preparation
    final_pdb_prep = run_glyco_prep(prep_args)
    
    # Step 2: Parametrization
    class ParamArgs:
        pass
    
    param_args = ParamArgs()
    param_args.prep_output_dir = args.prep_output_dir
    param_args.input_pdb = final_pdb_prep
    param_args.output_dir = args.param_output_dir
    param_args.skip_charmm_download = getattr(args, 'skip_charmm_download', False)
    
    final_pdb_param = run_glyco_param(param_args)
    
    if final_pdb_param is None:
        print("ERROR: Parametrization failed, stopping pipeline")
        return {"prep": final_pdb_prep, "param": None, "orient": None}
    
    # Step 3: Orientation
    class OrientArgs:
        pass
    
    orient_args = OrientArgs()
    orient_args.input_pdb = final_pdb_param
    orient_args.param_output_dir = args.param_output_dir
    orient_args.output_dir = args.orient_output_dir
    orient_args.charmm_dir = getattr(args, 'charmm_dir', None)
    orient_args.theta_step = args.theta_step
    orient_args.n_steps = args.n_steps
    orient_args.max_cycles = args.max_cycles
    orient_args.radius = args.radius_orient
    orient_args.use_coulomb = args.use_coulomb
    orient_args.n_workers = args.n_workers
    orient_args.save_individual_glycans = args.save_individual_glycans
    orient_args.save_before_after = args.save_before_after
    orient_args.verbose = args.verbose
    
    final_pdb_orient = run_glyco_orient(orient_args)
    
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print(f"\nResults:")
    print(f"  Preparation:    {final_pdb_prep}")
    print(f"  Parametrization: {final_pdb_param}")
    print(f"  Orientation:     {final_pdb_orient}")
    
    return {
        "prep": final_pdb_prep,
        "param": final_pdb_param,
        "orient": final_pdb_orient
    }


def main():
    """Main entry point for CLI."""
    if len(sys.argv) < 2:
        print("Usage: glyco-{prep,param,orient,all} [arguments]")
        print("\nAvailable commands:")
        print("  glyco-prep    - Run glycosylation preparation (Step 1)")
        print("  glyco-param   - Run parametrization (Step 2)")
        print("  glyco-orient  - Run carbohydrate orientation (Step 3)")
        print("  glyco-all     - Run complete pipeline (Steps 1-3)")
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
    elif command == "all":
        run_glyco_all()
    else:
        print(f"Unknown command: {command}", file=sys.stderr)
        print("Available commands: prep, param, orient, all", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
