#!/usr/bin/env python3
"""
Glycosylation Pipeline CLI - Unified interface for all steps
Usage: auto_glyco [options]
"""

import sys
import os
import argparse
import subprocess
import shutil
from pathlib import Path
from datetime import datetime

__all__ = ['main', 'auto_glyco', 'run_auto_glyco']

# Get the directory where this script is located
SCRIPT_DIR = Path(__file__).parent.absolute()
SCRIPTS_DIR = SCRIPT_DIR / "python_scripts"


def find_script(script_name):
    """Find script in python_scripts subdirectories"""
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


def auto_glyco():
    """Run complete glycosylation pipeline (equivalent to run_all.sh)"""
    parser = argparse.ArgumentParser(
        description="Automated Glycosylation Pipeline - Complete workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  auto_glyco --pdb protein.pdb --tsv glycans.tsv
  auto_glyco --pdb protein.pdb --tsv glycans.tsv --output_dir ./results --verbose
  auto_glyco --pdb protein.pdb --tsv glycans.tsv --theta_step 5 --max_cycles 10
        """
    )
    
    # Required arguments
    parser.add_argument("--pdb", required=True, help="Input PDB file path")
    parser.add_argument("--tsv", required=True, help="Input TSV file path")
    
    # Optional arguments
    parser.add_argument("--output_dir", default="./output", 
                       help="Output directory (default: ./output)")
    parser.add_argument("--url_charmm36", 
                       default="https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz",
                       help="URL for CHARMM36 force field download")
    
    # Step 3 parameters
    parser.add_argument("--theta_step", type=int, default=10,
                       help="Rotation step in degrees (default: 10)")
    parser.add_argument("--n_steps", type=int, default=10,
                       help="Number of rotation steps (default: 10)")
    parser.add_argument("--max_cycles", type=int, default=5,
                       help="Maximum MCMC cycles (default: 5)")
    parser.add_argument("--radius", type=float, default=300,
                       help="Radius for clash detection (default: 300)")
    parser.add_argument("--use_coulomb", choices=["yes", "no"], default="no",
                       help="Use Coulomb interactions (default: no)")
    parser.add_argument("--n_workers", type=int, default=12,
                       help="Number of parallel workers (default: 12)")
    parser.add_argument("--report_file", 
                       help="Report file path (default: <output_dir>/STEP3/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/report.txt)")
    
    # Flags
    parser.add_argument("--save_individual_glycans", action="store_true",
                       help="Save individual glycan PDB files")
    parser.add_argument("--save_before_after", action="store_true",
                       help="Save before/after comparison files")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    
    # Protein residue start (for Step 1)
    parser.add_argument("--protein_residue_start", type=int, default=10,
                       help="Protein residue start number (default: 10)")
    
    args = parser.parse_args()
    
    # Validate input files
    input_pdb = Path(args.pdb)
    input_tsv = Path(args.tsv)
    
    if not input_pdb.exists():
        print(f"ERROR: Input PDB not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    if not input_tsv.exists():
        print(f"ERROR: Input TSV not found: {input_tsv}", file=sys.stderr)
        sys.exit(1)
    
    # Convert to absolute paths
    input_pdb = input_pdb.absolute()
    input_tsv = input_tsv.absolute()
    output_base = Path(args.output_dir).absolute()
    
    # Set report file
    if args.report_file:
        report_file = Path(args.report_file).absolute()
    else:
        report_file = output_base / "STEP3" / "PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED" / "report.txt"
    
    # Create output directories (matching run_all.sh structure)
    step1_output = output_base / "STEP1"
    tsv_dir = step1_output / "TSV"
    pdb_dir = step1_output / "PDB_PROTEIN_GLYCOSYLATED"
    glycan_dir = step1_output / "EXTRACTED_CARBOHYDRATES"
    to_top_dir = step1_output / "TO_TOP"
    
    step2_output = output_base / "STEP2"
    topo_dir = step2_output
    json_dir = topo_dir / "JSON"
    pdb_glyco_dir = topo_dir / "PDB_GLYCOPROTEIN"
    variants_dir = topo_dir / "VALENCE_GLYCAN_VARIANTS"
    
    step3_output = output_base / "STEP3"
    json_files_dir = step3_output / "JSON_FILES"
    optimized_dir = step3_output / "PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED"
    glycans_only_dir = optimized_dir / "PDB_CARB_ONLY"
    
    # Create all directories
    for dir_path in [tsv_dir, pdb_dir, glycan_dir, to_top_dir, json_dir, 
                     pdb_glyco_dir, variants_dir, json_files_dir, optimized_dir, 
                     glycans_only_dir]:
        dir_path.mkdir(parents=True, exist_ok=True)
    
    # Get base filenames
    basename = input_tsv.stem
    basename_pdb = input_pdb.stem
    
    print("========================================")
    print("Automated Glycosylation Pipeline")
    print("========================================")
    print(f"Input PDB: {input_pdb}")
    print(f"Input TSV: {input_tsv}")
    print(f"Output directory: {output_base}")
    print(f"CHARMM36 URL: {args.url_charmm36}")
    print("========================================")
    
    # ========================================================================
    # STEP 1: GLYCOSYLATION PREPARATION
    # ========================================================================
    print("\n========================================")
    print("STEP 1: GLYCOSYLATION PREPARATION")
    print("========================================\n")
    
    # Step 1.0: Correct input table
    script_correct = find_script("0-correcting_caselino_table_for_variants.py")
    output_tsv = tsv_dir / f"{basename}_corrected.tsv"
    
    print("Step 1.0: Correcting input table...")
    subprocess.run([
        "python3", script_correct,
        "--input", str(input_tsv),
        "--output", str(output_tsv)
    ], check=True)
    
    # Step 1.1: Convert to IUPAC notation
    script_iupac = find_script("1-iupac_converted.py")
    output_glycosylator = tsv_dir / f"{basename}_glycosylator.tsv"
    
    print("\nStep 1.1: Converting to IUPAC notation...")
    subprocess.run([
        "python3", script_iupac,
        "--input_tsv", str(output_tsv),
        "--output_tsv", str(output_glycosylator),
        "--output_dir", str(glycan_dir)
    ], check=True)
    
    # Step 1.2a: Preparation of asparagine orientations
    script_asn = find_script("asn_orientation.py")
    pdb_asn_output = pdb_dir / f"{basename_pdb}_asn_orientation.pdb"
    
    print("\nStep 1.2a: Optimizing asparagine orientations...")
    subprocess.run([
        "python3", script_asn, str(input_pdb),
        "--rotate-atoms", "OD1,CG,ND2,HD22,HD21,HB2,HB3",
        "--fixed-atom", "CB",
        "--center-atom", "CA",
        "--radius", "30.0",
        "--rotation-step", "1",
        "-o", str(pdb_asn_output)
    ], check=True)
    
    # Step 1.2b: Run glycosylation script
    script_glycosylation = find_script("2-glycosylation_script.py")
    pdb_glycosylated = pdb_dir / f"{basename_pdb}_glycosylated.pdb"
    
    print("\nStep 1.2b: Running glycosylation script...")
    subprocess.run([
        "python3", script_glycosylation,
        "--input_tsv_glycosylator", str(output_glycosylator),
        "--input_pdb_protein", str(pdb_asn_output),
        "--protein_residue_start", str(args.protein_residue_start),
        "--output", str(pdb_glycosylated)
    ], check=True)
    
    # Step 1.3: Correct chain labels and residue numbers
    script_correction = find_script("3-correction_chain_labels_residue_numbers.py")
    pdb_renamed = pdb_dir / f"{basename_pdb}_glycosylated_renumbered.pdb"
    
    print("\nStep 1.3: Correcting chain labels and residue numbers...")
    subprocess.run([
        "python3", script_correction,
        str(pdb_glycosylated),
        str(pdb_renamed)
    ], check=True)
    
    # Step 1.4: Extract coordinates of glycans from PDB
    script_extract = find_script("4-extract_coordinates_of_glycans_from_structure.py")
    pdb_noh = pdb_dir / f"{basename_pdb}_glycosylated_renumbered_without_H.pdb"
    
    print("\nStep 1.4: Extracting glycans coordinates...")
    subprocess.run([
        "python3", script_extract,
        "--input_pdb", str(pdb_renamed),
        "--output_noH", str(pdb_noh),
        "--output_dir", str(to_top_dir)
    ], check=True)
    
    print("\nSTEP 1 completed successfully!")
    
    # ========================================================================
    # STEP 2: PARAMETRIZATION SCRIPTS
    # ========================================================================
    print("\n========================================")
    print("STEP 2: PARAMETRIZATION SCRIPTS")
    print("========================================\n")
    
    # Change to topo directory for CHARMM operations
    original_cwd = os.getcwd()
    os.chdir(topo_dir)
    
    # Setup CHARMM force field
    charmm_dir = topo_dir / "charmm36.ff"
    charmm_rtp = charmm_dir / "carb.rtp"
    charmm_hdb = charmm_dir / "carb.hdb"
    charmm_rtp_backup = charmm_dir / "carb.rtp.backup"
    charmm_hdb_backup = charmm_dir / "carb.hdb.backup"
    
    # Download or restore CHARMM force field
    if not charmm_hdb_backup.exists() or not charmm_rtp_backup.exists():
        print("Backup files not found. Downloading CHARMM force field...")
        
        if charmm_dir.exists():
            shutil.rmtree(charmm_dir)
        
        # Download and extract
        subprocess.run(["wget", args.url_charmm36, "-O", "charmm36.ff.tgz"], check=True)
        subprocess.run(["tar", "-xzf", "charmm36.ff.tgz"], check=True)
        
        # Find extracted directory
        extracted_dirs = list(topo_dir.glob("charmm36-*")) + list(topo_dir.glob("charmm36*"))
        if extracted_dirs:
            extracted_dirs[0].rename(charmm_dir)
        else:
            print("ERROR: Could not find extracted CHARMM directory", file=sys.stderr)
            sys.exit(1)
        
        # Remove tarball
        (topo_dir / "charmm36.ff.tgz").unlink()
        
        # Create backups
        shutil.copy(charmm_hdb, charmm_hdb_backup)
        shutil.copy(charmm_rtp, charmm_rtp_backup)
        print(f"CHARMM force field downloaded and backups created at: {charmm_dir}")
    else:
        print("Backup files found. Restoring from backups...")
        charmm_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(charmm_hdb_backup, charmm_hdb)
        shutil.copy(charmm_rtp_backup, charmm_rtp)
        print("Force field restored from backups.")
    
    # Restore residuetypes.dat if backup exists (optional, continue if fails)
    residuetypes_backup = Path.home() / "programs/GROMACS/share/gromacs/top/residuetypes.dat.backup"
    if residuetypes_backup.exists():
        try:
            shutil.copy(residuetypes_backup, 
                       Path.home() / "programs/GROMACS/share/gromacs/top/residuetypes.dat")
        except:
            pass  # Non-critical
    
    # Step 2.1: Generate JSON for each glycan
    script_json_gen = find_script("0-JSON_generator.py")
    print("\nGenerating JSON for each glycan...")
    subprocess.run([
        "python3", script_json_gen,
        "--base_dir", str(step1_output),
        "--output_dir", str(json_dir)
    ], check=True)
    
    # Process each glycan directory
    script_parser = find_script("1-parser_pdb.py")
    script_carb_rtp = find_script("2-parser_carb_rtp.py")
    script_compare = find_script("3-comparison_pdb_rtp.py")
    script_rtp1 = find_script("4-rtp_generator_part1.py")
    script_rtp2 = find_script("4-rtp_generator_part2.py")
    script_rtp3 = find_script("4-rtp_generator_part3.py")
    script_acetyl = find_script("5-acetylation_replacement.py")
    script_clean = find_script("6-clean_rtp.py")
    
    print("\nStarting processing for each glycan directory...")
    
    for glycan_dir_path in json_dir.iterdir():
        if not glycan_dir_path.is_dir():
            continue
        
        dir_basename = glycan_dir_path.name
        pdb_file = glycan_dir_path / f"{dir_basename}.pdb"
        json_file = glycan_dir_path / f"{dir_basename}.json"
        parser_file = glycan_dir_path / f"{dir_basename}_parser.pkl"
        rtp_pickle = glycan_dir_path / "carb_residues.pkl"
        rtp_modified = glycan_dir_path / "carb_modified.rtp"
        rtp_unique = glycan_dir_path / "carb_unique.rtp"
        
        if not pdb_file.exists():
            print(f"WARNING: PDB file not found: {pdb_file}, skipping...")
            continue
        
        print(f"\n----------------------------------------")
        print(f"Processing: {dir_basename}")
        
        os.chdir(glycan_dir_path)
        
        try:
            subprocess.run(["python3", script_parser, str(pdb_file), "-o", str(parser_file)], check=True)
            subprocess.run(["python3", script_carb_rtp, str(charmm_rtp), "-o", str(rtp_pickle)], check=True)
            subprocess.run(["python3", script_compare, "--pdb", str(parser_file), "--rtp", str(rtp_pickle)], check=True)
            subprocess.run(["python3", script_rtp1, "-p", str(parser_file), "-r", str(rtp_pickle)], check=True)
            subprocess.run(["python3", script_rtp2, "--pdb", str(parser_file), "--rtp", str(rtp_pickle), "--json", str(json_file)], check=True)
            subprocess.run(["python3", script_rtp3, "--pdb", str(parser_file), "--rtp", str(rtp_pickle), "--json", str(json_file), "--output", str(rtp_modified)], check=True)
            subprocess.run(["python3", script_acetyl, str(pdb_file), str(glycan_dir_path / f"{dir_basename}_modified.pdb")], check=True)
            
            dir_letter = dir_basename[0]
            subprocess.run(["python3", script_clean, str(rtp_modified), str(rtp_unique), dir_letter], check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR processing {dir_basename}, skipping...")
            os.chdir(topo_dir)
            continue
        
        os.chdir(topo_dir)
        print(f"Finished {dir_basename}")
    
    # Unification of RTP / HDB files
    script_together1 = find_script("7-together_part_1.py")
    script_together2 = find_script("7-together_part_2.py")
    script_together3 = find_script("7-together_part_3.py")
    script_together4 = find_script("7-together_part_4.py")
    script_together5 = find_script("7-together_part_5.py")
    
    print("\nUnifying RTP files...")
    
    unique_rtp_files = list(json_dir.glob("*/carb_unique.rtp"))
    if unique_rtp_files:
        subprocess.run(["python3", script_together1, "--input", str(json_dir), "--output", str(json_dir / "carb_unique_total.rtp")], check=True)
        subprocess.run(["python3", script_together2, "--input", str(json_dir / "carb_unique_total.rtp"), "--output", str(json_dir / "carb_redundance_removed.rtp")], check=True)
        subprocess.run(["python3", script_together3, "--input", str(json_dir / "carb_redundance_removed.rtp"), "--output", str(charmm_rtp)], check=True)
        subprocess.run(["python3", script_together4, "--input", str(json_dir / "carb_redundance_removed.rtp")], check=True)
        subprocess.run(["python3", script_together5, "--input", str(json_dir / "carb_redundance_removed.rtp"), "--output", str(json_dir / "carb_redundance_removed.hdb")], check=True)
    else:
        print("WARNING: No carb_unique.rtp files found. Skipping unification steps.")
    
    # Final Glycoprotein Construction
    script_glycoprot = find_script("8-glycoprotein.py")
    script_connect = find_script("9-conection_glycosilation_without_TER.py")
    script_identify = find_script("glycosylation_identifying.py")
    
    print("\nBuilding final glycoprotein...")
    
    input_protein = pdb_renamed
    if not input_protein.exists():
        protein_candidates = list(pdb_dir.glob("*_glycosylated_renumbered.pdb"))
        if protein_candidates:
            input_protein = protein_candidates[0]
        else:
            print(f"ERROR: Could not find *_glycosylated_renumbered.pdb in {pdb_dir}", file=sys.stderr)
            sys.exit(1)
    
    print(f"Found input protein: {input_protein}")
    
    output_protein = pdb_glyco_dir / f"{basename_pdb}_glycosylated_corrected.pdb"
    final_structure = pdb_glyco_dir / f"{basename_pdb}_glycosylated_final_connected.pdb"
    final_structure_2 = pdb_glyco_dir / f"{basename_pdb}_glycosylated_final_valence_corrected.pdb"
    final_structure_2_noh = pdb_glyco_dir / f"{basename_pdb}_glycosylated_final_valence_corrected_noh.pdb"
    
    glycan_dirs = list(json_dir.glob("*/"))
    if glycan_dirs:
        subprocess.run([
            "python3", script_glycoprot,
            "--protein", str(input_protein),
            "--carbs_dir", str(json_dir),
            "--output", str(output_protein),
            "--keep_hydrogens_carb",
            "--keep_hydrogens_prot"
        ], check=True)
        
        subprocess.run([
            "python3", script_connect,
            "--glycosylated", str(output_protein),
            "--conect", str(input_protein),
            "--output", str(final_structure)
        ], check=True)
        
        subprocess.run([
            "python3", script_identify,
            str(final_structure),
            str(final_structure_2)
        ], check=True)
        
        if not final_structure_2.exists() and final_structure_2_noh.exists():
            final_structure_2 = final_structure_2_noh
        
        print(f"Final structure file: {final_structure_2}")
    else:
        print("WARNING: No glycan directories found, skipping glycoprotein construction steps.")
    
    # HDB Update
    script_hdb = find_script("10-generation_hdb.py")
    
    print("\nUpdating HDB files...")
    
    if (json_dir / "carb_redundance_removed.rtp").exists():
        subprocess.run([
            "python3", script_hdb,
            str(json_dir / "carb_redundance_removed.rtp"),
            str(charmm_hdb),
            "-o", str(json_dir / "carb_modified.hdb")
        ], check=True)
    else:
        print("WARNING: carb_redundance_removed.rtp not found, skipping HDB generation")
    
    # Variants Generation
    script_variants = find_script("glycosylation_variants.py")
    
    print("\nGenerating glycan variants...")
    
    if final_structure_2.exists() and (json_dir / "carb_redundance_removed.rtp").exists() and (json_dir / "carb_redundance_removed.hdb").exists():
        try:
            subprocess.run([
                "python3", script_variants,
                "-p", str(final_structure_2),
                "-r", str(json_dir / "carb_redundance_removed.rtp"),
                "-d", str(json_dir / "carb_redundance_removed.hdb"),
                "-o", str(variants_dir)
            ], check=True)
        except subprocess.CalledProcessError:
            print("WARNING: glycosylation_variants.py failed, but continuing...")
    else:
        print("WARNING: Required files for variant generation not found")
    
    # Include variants in CHARMM files
    print("\nIncluding generated variants in CHARMM force field files...")
    
    variant_rtp = next(variants_dir.glob("*_variants.rtp"), None) or next(variants_dir.glob("*.rtp"), None)
    variant_hdb = next(variants_dir.glob("*_variants.hdb"), None) or next(variants_dir.glob("*.hdb"), None)
    
    if variant_rtp and variant_hdb:
        print(f"Found variant RTP: {variant_rtp}")
        print(f"Found variant HDB: {variant_hdb}")
        
        with open(charmm_rtp, 'a') as f:
            f.write(variant_rtp.read_text())
        
        with open(charmm_hdb, 'a') as f:
            f.write(variant_hdb.read_text())
        
        shutil.copy(charmm_hdb, charmm_hdb_backup)
        shutil.copy(charmm_rtp, charmm_rtp_backup)
        print("Variants successfully included in CHARMM force field files.")
    else:
        print("WARNING: Variant files not found, skipping inclusion.")
    
    print("\nSTEP 2 completed successfully!")
    
    # ========================================================================
    # STEP 3: CARBOHYDRATE ORIENTATION
    # ========================================================================
    print("\n========================================")
    print("STEP 3: CARBOHYDRATE ORIENTATION")
    print("========================================\n")
    
    script_pdb_to_json = find_script("1-pdb_to_json.py")
    script_add_charmm = find_script("3-adding_chamm36_parameters.py")
    script_optimize = find_script("4-optimize_glycans_mcmc.py")
    
    # Find the glycoprotein variants directory
    glycoprotein_step2_dir = variants_dir if variants_dir.exists() else json_dir
    if not glycoprotein_step2_dir.exists():
        glycoprotein_step2_dir = json_dir
    
    # Find the protein PDB file
    protein_step2 = next(glycoprotein_step2_dir.glob("*_variants.pdb"), None)
    if not protein_step2:
        protein_step2 = next(pdb_glyco_dir.glob("*_final_valence_corrected.pdb"), None)
    if not protein_step2:
        protein_step2 = next(pdb_glyco_dir.glob("*.pdb"), None)
        if protein_step2 and "without_H" in str(protein_step2):
            protein_step2 = next(pdb_glyco_dir.glob("*.pdb"), None)
    
    if not protein_step2 or not protein_step2.exists():
        print(f"ERROR: Could not find protein PDB file in {step2_output}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found input protein for Step 3: {protein_step2}")
    
    # Extract basename for output files
    step3_basename = protein_step2.stem
    for suffix in ["_variants", "_final_valence_corrected", "_glycosylated_final_valence_corrected"]:
        step3_basename = step3_basename.replace(suffix, "")
    
    print(f"Using basename: {step3_basename}")
    
    # Define output files for Step 3
    pdb_to_json_output = json_files_dir / "pdb_to_json.json"
    charmm36_json_output = json_files_dir / "glycan_data_charmm36.json"
    optimized_json_output = optimized_dir / "glycan_optimized.json"
    optimized_pdb_output = optimized_dir / f"{step3_basename}_optimized.pdb"
    
    # Step 3.1: Convert PDB to JSON
    print("\nStep 3.1: Converting PDB to JSON...")
    subprocess.run([
        "python3", script_pdb_to_json,
        "--input_pdb", str(protein_step2),
        "--output_json", str(pdb_to_json_output)
    ], check=True)
    
    # Step 3.2: Add CHARMM36 parameters
    print("\nStep 3.2: Adding CHARMM36 parameters...")
    subprocess.run([
        "python3", script_add_charmm,
        "--input_json", str(pdb_to_json_output),
        "--charmm_dir", str(charmm_dir),
        "--output_json", str(charmm36_json_output)
    ], check=True)
    
    # Step 3.3: Optimize glycans using MCMC
    print("\nStep 3.3: Optimizing glycans using MCMC...")
    
    cmd = [
        "python3", script_optimize,
        "--input_json", str(charmm36_json_output),
        "--output_json", str(optimized_json_output),
        "--output_pdb", str(optimized_pdb_output),
        "--glycans_output_dir", str(glycans_only_dir),
        "--theta_step", str(args.theta_step),
        "--n_steps", str(args.n_steps),
        "--max_cycles", str(args.max_cycles),
        "--radius", str(args.radius),
        "--use_coulomb", args.use_coulomb,
        "--n_workers", str(args.n_workers),
        "--report_file", str(report_file)
    ]
    
    if args.save_individual_glycans:
        cmd.append("--save_individual_glycans")
    if args.save_before_after:
        cmd.append("--save_before_after")
    if args.verbose:
        cmd.append("--verbose")
    
    subprocess.run(cmd, check=True)
    
    # Change back to original directory
    os.chdir(original_cwd)
    
    print("\nSTEP 3 completed successfully!")
    
    # ========================================================================
    # COMPLETION
    # ========================================================================
    print("\n========================================")
    print("All steps completed successfully!")
    print("========================================")
    print(f"Results saved in: {output_base}")
    print("")
    print("STEP 1 (Glycosylation Preparation):")
    print(f"  - TSV files: {tsv_dir}")
    print(f"  - PDB files: {pdb_dir}")
    print(f"  - Glycans: {glycan_dir}")
    print("")
    print("STEP 2 (Parametrization):")
    print(f"  - JSON files: {json_dir}")
    print(f"  - Glycoprotein PDB: {pdb_glyco_dir}")
    print(f"  - Variants: {variants_dir}")
    print(f"  - CHARMM36: {charmm_dir}")
    print("")
    print("STEP 3 (Carbohydrate Orientation):")
    print(f"  - JSON files: {json_files_dir}")
    print(f"  - Optimized structures: {optimized_dir}")
    print(f"  - Individual glycans: {glycans_only_dir}")
    print("")
    print("Main output files:")
    print(f"  - Optimized PDB: {optimized_pdb_output}")
    print(f"  - Report: {report_file}")
    print("========================================")


def main():
    """Main entry point for CLI."""
    auto_glyco()


if __name__ == "__main__":
    main()
