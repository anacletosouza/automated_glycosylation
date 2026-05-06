#!/usr/bin/env python3
"""
Command-line interface for automated glycosylation pipeline.

Author: Anacleto Silva de Souza
License: MIT
"""

import argparse
import subprocess
import sys
import os
import shutil
import urllib.request
import tarfile
from pathlib import Path
from typing import List, Optional, Dict, Any

def get_package_dir() -> Path:
    """Get the package installation directory."""
    return Path(__file__).parent

def get_script_path(script_name: str, subdir: str) -> str:
    """
    Get the full path to a script in the package.
    
    Args:
        script_name: Name of the script file
        subdir: Subdirectory within scripts (e.g., '1_glycosylation_preparation')
    """
    package_dir = get_package_dir()
    script_path = package_dir / "scripts" / subdir / script_name
    
    if not script_path.exists():
        # Try with .py extension
        script_path = package_dir / "scripts" / subdir / f"{script_name}.py"
    
    if not script_path.exists():
        raise FileNotFoundError(f"Script not found: {script_name} in {subdir}")
    
    return str(script_path)

def run_python_script(script_name: str, args: List[str], subdir: str) -> subprocess.CompletedProcess:
    """Run a Python script with the given arguments."""
    script_path = get_script_path(script_name, subdir)
    cmd = [sys.executable, script_path] + args
    return subprocess.run(cmd, capture_output=False)

def create_directory(path: Path) -> None:
    """Create directory if it doesn't exist."""
    path.mkdir(parents=True, exist_ok=True)

def download_charmm(output_dir: Path, charmm_url: str = None) -> Path:
    """
    Download CHARMM36 force field.
    
    Args:
        output_dir: Directory to download to
        charmm_url: Custom URL for CHARMM download
    
    Returns:
        Path to CHARMM directory
    """
    if charmm_url is None:
        charmm_url = "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz"
    
    charmm_dir = output_dir / "charmm36.ff"
    charmm_rtp = charmm_dir / "carb.rtp"
    charmm_hdb = charmm_dir / "carb.hdb"
    charmm_rtp_backup = charmm_dir / "carb.rtp.backup"
    charmm_hdb_backup = charmm_dir / "carb.hdb.backup"
    
    # Check if already exists and we have backups
    if charmm_dir.exists() and charmm_rtp_backup.exists() and charmm_hdb_backup.exists():
        print("CHARMM directory with backups already exists. Using existing.")
        return charmm_dir
    
    print("Downloading CHARMM36 force field...")
    
    # Remove existing directory if it exists without backups
    if charmm_dir.exists():
        shutil.rmtree(charmm_dir)
    
    # Download and extract
    tgz_file = output_dir / "charmm36.ff.tgz"
    urllib.request.urlretrieve(charmm_url, tgz_file)
    
    with tarfile.open(tgz_file, "r:gz") as tar:
        tar.extractall(output_dir)
    
    # Rename extracted directory (the name may vary)
    for item in output_dir.iterdir():
        if item.is_dir() and "charmm36" in item.name.lower() and item != charmm_dir:
            item.rename(charmm_dir)
            break
    
    tgz_file.unlink()  # Remove tar file
    
    # Create backups
    if charmm_rtp.exists():
        shutil.copy2(charmm_rtp, charmm_rtp_backup)
        shutil.copy2(charmm_hdb, charmm_hdb_backup)
    
    print("CHARMM36 download complete")
    return charmm_dir

def run_glyco_prep():
    """Run glycosylation preparation pipeline (Step 1)."""
    parser = argparse.ArgumentParser(
        description="Step 1: Glycosylation preparation - Adds glycans to protein structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  glyco-prep -i input.pdb -o output_dir
  
  # With custom parameters
  glyco-prep -i input.pdb -o output_dir --glycan_sites_tsv glycosylation_sites.tsv
        """
    )
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True,
                        help="Input PDB file (protein structure)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for results")
    
    # Optional arguments for glycosylation sites
    parser.add_argument("--glycan_sites_tsv", default=None,
                        help="TSV file with glycosylation sites (optional, auto-detect if not provided)")
    
    # Optional arguments for asparagine orientation
    parser.add_argument("--rotate-atoms", 
                        default="OD1,CG,ND2,HD22,HD21,HB2,HB3",
                        help="Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)")
    parser.add_argument("--fixed-atom",
                        default="CB",
                        help="Fixed atom for rotation (default: CB)")
    parser.add_argument("--center-atom",
                        default="CA",
                        help="Center atom for rotation (default: CA)")
    parser.add_argument("--radius",
                        type=float,
                        default=30.0,
                        help="Radius for neighbor detection (default: 30.0)")
    parser.add_argument("--rotation-step",
                        type=int,
                        default=1,
                        help="Rotation step in degrees (default: 1)")
    
    # Optional arguments for protein residue start
    parser.add_argument("--protein-residue-start",
                        type=int,
                        default=1,
                        help="Starting residue number for protein (default: 1)")
    
    # Other options
    parser.add_argument("--keep-temp",
                        action="store_true",
                        help="Keep temporary files")
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb = Path(args.input).resolve()
    output_path = Path(args.output_dir).resolve()
    
    if not input_pdb.exists():
        print(f"Error: Input PDB file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directories
    pdb_glycosylated = output_path / "PDB_PROTEIN_GLYCOSYLATED"
    tsv_dir = output_path / "TSV"
    extracted_dir = output_path / "EXTRACTED_CARBOHYDRATES"
    to_top_dir = output_path / "TO_TOP"
    
    for d in [pdb_glycosylated, tsv_dir, extracted_dir, to_top_dir]:
        create_directory(d)
    
    # Step 0: Correct Caselino table (if TSV provided)
    if args.asn_tsv and Path(args.asn_tsv).exists():
        asn_tsv_path = Path(args.asn_tsv).resolve()
        corrected_tsv = tsv_dir / "glycosylation_sites_corrected.tsv"
        
        print("Step 0: Correcting glycosylation sites table...")
        result = run_python_script("0-correcting_caselino_table_for_variants.py",
                                   ["--input", str(asn_tsv_path), "--output", str(corrected_tsv)],
                                   "1_glycosylation_preparation")
        if result.returncode != 0:
            print("Error in table correction", file=sys.stderr)
            sys.exit(result.returncode)
        
        # Step 1: Convert to IUPAC notation
        glycosylator_tsv = tsv_dir / "glycosylation_sites_glycosylator.tsv"
        
        print("Step 1: Converting to IUPAC notation...")
        result = run_python_script("1-iupac_converted.py",
                                   ["--input_tsv", str(corrected_tsv), 
                                    "--output_tsv", str(glycosylator_tsv),
                                    "--output_dir", str(extracted_dir)],
                                   "1_glycosylation_preparation")
        if result.returncode != 0:
            print("Error in IUPAC conversion", file=sys.stderr)
            sys.exit(result.returncode)
        
        asn_input_tsv = str(glycosylator_tsv)
    else:
        asn_input_tsv = None
    
    # Step 2a: Optimize asparagine orientations
    asn_output = pdb_glycosylated / "protein_asn_orientation.pdb"
    
    print("Step 2a: Optimizing asparagine orientations...")
    asn_args = [
        str(input_pdb),
        "--rotate-atoms", args.rotate_atoms,
        "--fixed-atom", args.fixed_atom,
        "--center-atom", args.center_atom,
        "--radius", str(args.radius),
        "--rotation-step", str(args.rotation_step),
        "-o", str(asn_output)
    ]
    
    result = run_python_script("asn_orientation.py", asn_args, "1_glycosylation_preparation")
    if result.returncode != 0:
        print("Error in asparagine orientation", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Step 2b: Run glycosylation script
    glycosylated_output = pdb_glycosylated / "protein_glycosylated.pdb"
    
    print("Step 2b: Running glycosylation script...")
    glyco_args = [
        "--input_pdb_protein", str(asn_output),
        "--protein_residue_start", str(args.protein_residue_start),
        "--output", str(glycosylated_output)
    ]
    
    if asn_input_tsv:
        glyco_args.extend(["--input_tsv_glycosylator", asn_input_tsv])
    
    result = run_python_script("2-glycosylation_script.py", glyco_args, "1_glycosylation_preparation")
    if result.returncode != 0:
        print("Error in glycosylation script", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Step 3: Correct chain labels and residue numbers
    renamed_pdb = pdb_glycosylated / "protein_glycosylated_renumbered.pdb"
    
    print("Step 3: Correcting chain labels and residue numbers...")
    correct_args = [str(glycosylated_output), str(renamed_pdb)]
    
    result = run_python_script("3-correction_chain_labels_residue_numbers.py", 
                               correct_args, "1_glycosylation_preparation")
    if result.returncode != 0:
        print("Error in chain correction", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Step 4: Extract coordinates of glycans
    output_noh_pdb = pdb_glycosylated / "protein_glycosylated_renumbered_without_H.pdb"
    
    print("Step 4: Extracting glycans coordinates...")
    extract_args = [
        "--input_pdb", str(renamed_pdb),
        "--output_noH", str(output_noh_pdb),
        "--output_dir", str(to_top_dir)
    ]
    
    result = run_python_script("4-extract_coordinates_of_glycans_from_structure.py",
                               extract_args, "1_glycosylation_preparation")
    if result.returncode != 0:
        print("Error in glycan extraction", file=sys.stderr)
        sys.exit(result.returncode)
    
    print(f"\nGlycosylation preparation completed successfully!")
    print(f"Results saved to: {args.output_dir}")
    print(f"Glycosylated protein: {renamed_pdb}")

def run_glyco_param():
    """Run parametrization pipeline (Step 2)."""
    parser = argparse.ArgumentParser(
        description="Step 2: Parametrization - Generate topology and force field parameters",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("-i", "--input-pdb", required=True,
                        help="Input PDB file from step 1 (glycosylated_protein_renumbered.pdb)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for topology files")
    
    # Optional CHARMM download
    parser.add_argument("--download-charmm",
                        action="store_true",
                        help="Download CHARMM36 force field (default: False, use existing)")
    parser.add_argument("--charmm-url",
                        default="https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz",
                        help="URL for CHARMM force field download")
    
    # Optional: path to existing CHARMM directory
    parser.add_argument("--charmm-dir",
                        default=None,
                        help="Path to existing CHARMM36 force field directory")
    
    # Optional: skip download if backup exists
    parser.add_argument("--force-download",
                        action="store_true",
                        help="Force download even if backup exists")
    
    # Optional: number of CPUs for parallel processing
    parser.add_argument("--n-cpus",
                        type=int,
                        default=1,
                        help="Number of CPUs for parallel processing (default: 1)")
    
    # Optional: keep intermediate files
    parser.add_argument("--keep-intermediate",
                        action="store_true",
                        help="Keep intermediate files")
    
    # Optional: glycosylation sites directory (TO_TOP from step 1)
    parser.add_argument("--glycans-dir",
                        default=None,
                        help="Directory containing extracted glycans (TO_TOP from step 1)")
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb = Path(args.input_pdb).resolve()
    output_path = Path(args.output_dir).resolve()
    
    if not input_pdb.exists():
        print(f"Error: Input PDB file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directories
    json_dir = output_path / "JSON"
    pdb_dir = output_path / "PDB_GLYCOPROTEIN"
    valence_dir = output_path / "VALENCE_GLYCAN_VARIANTS"
    
    for d in [json_dir, pdb_dir, valence_dir]:
        create_directory(d)
    
    # Determine CHARMM directory
    if args.charmm_dir:
        charmm_dir = Path(args.charmm_dir).resolve()
    else:
        charmm_dir = output_path / "charmm36.ff"
    
    # Download CHARMM if requested or if doesn't exist
    if args.download_charmm or args.force_download or not charmm_dir.exists():
        charmm_dir = download_charmm(output_path, args.charmm_url)
    
    charmm_rtp = charmm_dir / "carb.rtp"
    charmm_hdb = charmm_dir / "carb.hdb"
    charmm_rtp_backup = charmm_dir / "carb.rtp.backup"
    charmm_hdb_backup = charmm_dir / "carb.hdb.backup"
    
    # Restore from backups if they exist
    if charmm_rtp_backup.exists() and charmm_hdb_backup.exists():
        shutil.copy2(charmm_rtp_backup, charmm_rtp)
        shutil.copy2(charmm_hdb_backup, charmm_hdb)
        print("Restored CHARMM files from backups")
    
    if not charmm_rtp.exists() or not charmm_hdb.exists():
        print(f"Error: CHARMM files not found in {charmm_dir}", file=sys.stderr)
        print("Please provide --charmm-dir or use --download-charmm", file=sys.stderr)
        sys.exit(1)
    
    # Determine glycans directory
    if args.glycans_dir:
        glycans_dir = Path(args.glycans_dir).resolve()
    else:
        # Try to find TO_TOP from step 1
        step1_dir = input_pdb.parent.parent
        glycans_dir = step1_dir / "TO_TOP"
    
    if not glycans_dir.exists():
        print(f"Warning: Glycans directory not found at {glycans_dir}", file=sys.stderr)
        print("Please provide --glycans-dir with path to TO_TOP directory from step 1", file=sys.stderr)
    
    # Step: Generate JSONs for each glycan
    print("Generating JSON for each glycan...")
    
    json_gen_args = [
        "--base_dir", str(glycans_dir.parent) if glycans_dir.exists() else str(output_path),
        "--output_dir", str(json_dir)
    ]
    
    result = run_python_script("0-JSON_generator.py", json_gen_args, "2_parametrization_scripts")
    if result.returncode != 0:
        print("Error in JSON generation", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Process each glycan directory
    print("Processing each glycan...")
    
    for glycan_dir in json_dir.iterdir():
        if not glycan_dir.is_dir():
            continue
        
        basename = glycan_dir.name
        print(f"Processing: {basename}")
        
        pdb_file = glycan_dir / f"{basename}.pdb"
        json_file = glycan_dir / f"{basename}.json"
        parser_file = glycan_dir / f"{basename}_parser.pkl"
        rtp_pickle = glycan_dir / "carb_residues.pkl"
        rtp_modified = glycan_dir / "carb_modified.rtp"
        rtp_unique = glycan_dir / "carb_unique.rtp"
        
        # Run parser_pdb
        if pdb_file.exists():
            result = run_python_script("1-parser_pdb.py", [str(pdb_file), "-o", str(parser_file)], 
                                       "2_parametrization_scripts")
            if result.returncode != 0:
                continue
        
        # Run parser_carb_rtp
        result = run_python_script("2-parser_carb_rtp.py", [str(charmm_rtp), "-o", str(rtp_pickle)],
                                   "2_parametrization_scripts")
        if result.returncode != 0:
            continue
        
        # Run comparison_pdb_rtp
        if parser_file.exists():
            result = run_python_script("3-comparison_pdb_rtp.py",
                                       ["--pdb", str(parser_file), "--rtp", str(rtp_pickle)],
                                       "2_parametrization_scripts")
        
        # Run rtp_generator_part1
        if parser_file.exists():
            result = run_python_script("4-rtp_generator_part1.py",
                                       ["-p", str(parser_file), "-r", str(rtp_pickle)],
                                       "2_parametrization_scripts")
        
        # Run rtp_generator_part2
        if parser_file.exists() and json_file.exists():
            result = run_python_script("4-rtp_generator_part2.py",
                                       ["--pdb", str(parser_file), "--rtp", str(rtp_pickle), 
                                        "--json", str(json_file)],
                                       "2_parametrization_scripts")
        
        # Run rtp_generator_part3
        if parser_file.exists() and json_file.exists():
            result = run_python_script("4-rtp_generator_part3.py",
                                       ["--pdb", str(parser_file), "--rtp", str(rtp_pickle), 
                                        "--json", str(json_file), "--output", str(rtp_modified)],
                                       "2_parametrization_scripts")
        
        # Run acetylation_replacement
        if pdb_file.exists():
            modified_pdb = glycan_dir / f"{basename}_modified.pdb"
            result = run_python_script("5-acetylation_replacement.py",
                                       [str(pdb_file), str(modified_pdb)],
                                       "2_parametrization_scripts")
        
        # Run clean_rtp
        if rtp_modified.exists():
            dir_letter = basename[0] if basename else "a"
            result = run_python_script("6-clean_rtp.py",
                                       [str(rtp_modified), str(rtp_unique), dir_letter],
                                       "2_parametrization_scripts")
    
    # Unification steps
    print("Unifying RTP/HDB files...")
    
    result = run_python_script("7-together_part_1.py",
                               ["--input", str(json_dir), "--output", str(json_dir / "carb_unique_total.rtp")],
                               "2_parametrization_scripts")
    
    result = run_python_script("7-together_part_2.py",
                               ["--input", str(json_dir / "carb_unique_total.rtp"), 
                                "--output", str(json_dir / "carb_redundance_removed.rtp")],
                               "2_parametrization_scripts")
    
    result = run_python_script("7-together_part_3.py",
                               ["--input", str(json_dir / "carb_redundance_removed.rtp"), 
                                "--output", str(charmm_rtp)],
                               "2_parametrization_scripts")
    
    result = run_python_script("7-together_part_4.py",
                               ["--input", str(json_dir / "carb_redundance_removed.rtp")],
                               "2_parametrization_scripts")
    
    result = run_python_script("7-together_part_5.py",
                               ["--input", str(json_dir / "carb_redundance_removed.rtp"), 
                                "--output", str(json_dir / "carb_redundance_removed.hdb")],
                               "2_parametrization_scripts")
    
    # Final glycoprotein
    output_protein = pdb_dir / "glycoprotein_corrected.pdb"
    
    result = run_python_script("8-glycoprotein.py",
                               ["--protein", str(input_pdb), "--carbs_dir", str(json_dir),
                                "--output", str(output_protein), "--keep_hydrogens_carb", "--keep_hydrogens_prot"],
                               "2_parametrization_scripts")
    
    # Connection steps
    final_structure = pdb_dir / "glycoprotein_final_connected.pdb"
    
    result = run_python_script("9-conection_glycosilation_without_TER.py",
                               ["--glycosylated", str(output_protein), "--conect", str(input_pdb),
                                "--output", str(final_structure)],
                               "2_parametrization_scripts")
    
    final_structure_2 = pdb_dir / "glycoprotein_final_valence_corrected.pdb"
    
    result = run_python_script("glycosylation_identifying.py",
                               [str(final_structure), str(final_structure_2)],
                               "2_parametrization_scripts")
    
    # Generation HDB
    result = run_python_script("10-generation_hdb.py",
                               [str(json_dir / "carb_redundance_removed.rtp"), str(charmm_hdb),
                                "-o", str(json_dir / "carb_modified.hdb")],
                               "2_parametrization_scripts")
    
    # Generate variants
    result = run_python_script("glycosylation_variants.py",
                               ["-p", str(final_structure_2), 
                                "-r", str(json_dir / "carb_redundance_removed.rtp"),
                                "-d", str(json_dir / "carb_redundance_removed.hdb"),
                                "-o", str(valence_dir)],
                               "2_parametrization_scripts")
    
    # Include variants in CHARMM files
    variant_rtp = valence_dir / f"{final_structure_2.stem}_variants.rtp"
    variant_hdb = valence_dir / f"{final_structure_2.stem}_variants.hdb"
    
    if variant_rtp.exists() and variant_hdb.exists():
        print("Appending variants to CHARMM files...")
        with open(charmm_rtp, 'a') as f:
            f.write(variant_rtp.read_text())
        with open(charmm_hdb, 'a') as f:
            f.write(variant_hdb.read_text())
        
        # Update backups
        shutil.copy2(charmm_rtp, charmm_rtp_backup)
        shutil.copy2(charmm_hdb, charmm_hdb_backup)
    
    print(f"\nParametrization completed successfully!")
    print(f"Results saved to: {args.output_dir}")

def run_glyco_orient():
    """Run carbohydrate orientation pipeline (Step 3)."""
    parser = argparse.ArgumentParser(
        description="Step 3: Carbohydrate orientation - Optimize glycan orientations using MCMC",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("-i", "--input-pdb", required=True,
                        help="Input PDB file from step 2 (final valence corrected PDB)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for optimized structures")
    
    # Optional CHARMM directory
    parser.add_argument("--charmm-dir",
                        default=None,
                        help="CHARMM36 force field directory (default: auto-detect)")
    
    # Optional MCMC parameters
    parser.add_argument("--theta-step",
                        type=int,
                        default=10,
                        help="Theta step for MCMC (default: 10)")
    parser.add_argument("--n-steps",
                        type=int,
                        default=10,
                        help="Number of MCMC steps (default: 10)")
    parser.add_argument("--max-cycles",
                        type=int,
                        default=5,
                        help="Maximum cycles for MCMC (default: 5)")
    parser.add_argument("--radius",
                        type=float,
                        default=300.0,
                        help="Radius for interactions (default: 300.0)")
    
    # Optional Coulomb option
    parser.add_argument("--use-coulomb",
                        choices=['yes', 'no'],
                        default='no',
                        help="Use Coulomb interactions (default: no)")
    
    # Optional number of workers
    parser.add_argument("--n-workers",
                        type=int,
                        default=1,
                        help="Number of workers for parallel processing (default: 1)")
    
    # Optional flags
    parser.add_argument("--save-individual-glycans",
                        action="store_true",
                        help="Save individual glycan PDB files")
    parser.add_argument("--save-before-after",
                        action="store_true",
                        help="Save before/after comparison files")
    parser.add_argument("--verbose",
                        action="store_true",
                        help="Verbose output")
    
    # Optional report file
    parser.add_argument("--report-file",
                        default=None,
                        help="Report file path (default: auto-generated)")
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb = Path(args.input_pdb).resolve()
    output_path = Path(args.output_dir).resolve()
    
    if not input_pdb.exists():
        print(f"Error: Input PDB file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directories
    json_dir = output_path / "JSON_FILES"
    pdb_optimized_dir = output_path / "PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED"
    pdb_carb_only = pdb_optimized_dir / "PDB_CARB_ONLY"
    
    for d in [json_dir, pdb_optimized_dir, pdb_carb_only]:
        create_directory(d)
    
    # Determine CHARMM directory
    if args.charmm_dir:
        charmm_dir = Path(args.charmm_dir).resolve()
    else:
        # Try to find from step 2 output
        step2_dir = input_pdb.parent.parent.parent
        charmm_dir = step2_dir / "charmm36.ff"
        if not charmm_dir.exists():
            charmm_dir = step2_dir / "2-GLYCOPROTEIN_TOPOLOGY" / "charmm36.ff"
    
    if not charmm_dir.exists():
        print(f"Warning: CHARMM directory not found at {charmm_dir}", file=sys.stderr)
        print("Please provide --charmm-dir argument", file=sys.stderr)
        sys.exit(1)
    
    # Step 1: Convert PDB to JSON
    print("Step 1: Converting PDB to JSON...")
    pdb_json = json_dir / "pdb_to_json.json"
    
    result = run_python_script("1-pdb_to_json.py",
                               ["--input_pdb", str(input_pdb), "--output_json", str(pdb_json)],
                               "3_carbohydrate_orientation")
    if result.returncode != 0:
        print("Error in PDB to JSON conversion", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Step 2: Add CHARMM36 parameters
    print("Step 2: Adding CHARMM36 parameters...")
    glycan_json = json_dir / "glycan_data_charmm36.json"
    
    result = run_python_script("3-adding_chamm36_parameters.py",
                               ["--input_json", str(pdb_json), "--charmm_dir", str(charmm_dir),
                                "--output_json", str(glycan_json)],
                               "3_carbohydrate_orientation")
    if result.returncode != 0:
        print("Error adding CHARMM parameters", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Step 3: Optimize glycans
    print("Step 3: Optimizing glycans using MCMC...")
    optimized_json = pdb_optimized_dir / "glycan_optimized.json"
    output_pdb = pdb_optimized_dir / "glycoprotein_optimized.pdb"
    report_file = args.report_file if args.report_file else pdb_optimized_dir / "report.txt"
    
    mcmc_args = [
        "--input_json", str(glycan_json),
        "--output_json", str(optimized_json),
        "--output_pdb", str(output_pdb),
        "--glycans_output_dir", str(pdb_carb_only),
        "--theta_step", str(args.theta_step),
        "--n_steps", str(args.n_steps),
        "--max_cycles", str(args.max_cycles),
        "--radius", str(args.radius),
        "--use_coulomb", args.use_coulomb,
        "--n_workers", str(args.n_workers),
        "--report_file", str(report_file)
    ]
    
    if args.save_individual_glycans:
        mcmc_args.append("--save_individual_glycans")
    if args.save_before_after:
        mcmc_args.append("--save_before_after")
    if args.verbose:
        mcmc_args.append("--verbose")
    
    result = run_python_script("4-optimize_glycans_mcmc.py", mcmc_args, "3_carbohydrate_orientation")
    if result.returncode != 0:
        print("Error in MCMC optimization", file=sys.stderr)
        sys.exit(result.returncode)
    
    print(f"\nCarbohydrate orientation completed successfully!")
    print(f"Results saved to: {args.output_dir}")
    print(f"Optimized structure: {output_pdb}")

def run_all_pipeline():
    """Run the complete pipeline (steps 1, 2, and 3)."""
    parser = argparse.ArgumentParser(
        description="Run complete automated glycosylation pipeline (Steps 1-3)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with required inputs
  glyco-all --input protein.pdb --output-dir results
  
  # With all options
  glyco-all --input protein.pdb --output-dir results --glycan_sites_tsv sites.tsv \\
            --download-charmm --n-cpus 4 --n-workers 4 --verbose
        """
    )
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True,
                        help="Input PDB file (protein structure)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for all results")
    
    # Optional arguments
    parser.add_argument("--glycan_sites_tsv", default=None,
                        help="TSV file with glycosylation sites")
    parser.add_argument("--download-charmm", action="store_true",
                        help="Download CHARMM36 force field")
    parser.add_argument("--charmm-url", default=None,
                        help="Custom URL for CHARMM download")
    parser.add_argument("--n-cpus", type=int, default=1,
                        help="Number of CPUs for parallel processing (default: 1)")
    parser.add_argument("--n-workers", type=int, default=1,
                        help="Number of workers for MCMC (default: 1)")
    parser.add_argument("--keep-temp", action="store_true",
                        help="Keep temporary files")
    parser.add_argument("--verbose", action="store_true",
                        help="Verbose output")
    
    # Asparagine orientation parameters
    parser.add_argument("--rotate-atoms", default="OD1,CG,ND2,HD22,HD21,HB2,HB3",
                        help="Atoms to rotate")
    parser.add_argument("--fixed-atom", default="CB",
                        help="Fixed atom for rotation")
    parser.add_argument("--center-atom", default="CA",
                        help="Center atom for rotation")
    parser.add_argument("--radius", type=float, default=30.0,
                        help="Radius for neighbor detection")
    parser.add_argument("--rotation-step", type=int, default=1,
                        help="Rotation step in degrees")
    
    # MCMC parameters
    parser.add_argument("--theta-step", type=int, default=10,
                        help="Theta step for MCMC")
    parser.add_argument("--n-steps", type=int, default=10,
                        help="Number of MCMC steps")
    parser.add_argument("--max-cycles", type=int, default=5,
                        help="Maximum cycles for MCMC")
    parser.add_argument("--mcmc-radius", type=float, default=300.0,
                        help="Radius for MCMC interactions")
    parser.add_argument("--use-coulomb", choices=['yes', 'no'], default='no',
                        help="Use Coulomb interactions")
    
    # Save options
    parser.add_argument("--save-individual-glycans", action="store_true",
                        help="Save individual glycan PDB files")
    parser.add_argument("--save-before-after", action="store_true",
                        help="Save before/after comparison files")
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb = Path(args.input).resolve()
    output_path = Path(args.output_dir).resolve()
    
    if not input_pdb.exists():
        print(f"Error: Input PDB file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    # Create main output directory
    step1_dir = output_path / "1-GLYCOPROTEIN_PREPARATION"
    step2_dir = output_path / "2-GLYCOPROTEIN_TOPOLOGY"
    step3_dir = output_path / "3-MINIMIZATION_CARBOHYDRATE"
    
    for d in [step1_dir, step2_dir, step3_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Step 1
    print("\n" + "="*60)
    print("STEP 1: GLYCOSYLATION PREPARATION")
    print("="*60)
    
    # Build step 1 command
    step1_cmd = [sys.executable, "-c", f"""
import sys
sys.argv = ['glyco-prep', '-i', '{input_pdb}', '-o', '{step1_dir}']
if '{args.asn_tsv}' and '{args.asn_tsv}' != 'None':
    sys.argv.extend(['--glycan_sites_tsv', '{args.asn_tsv}'])
sys.argv.extend(['--rotate-atoms', '{args.rotate_atoms}'])
sys.argv.extend(['--fixed-atom', '{args.fixed_atom}'])
sys.argv.extend(['--center-atom', '{args.center_atom}'])
sys.argv.extend(['--radius', '{args.radius}'])
sys.argv.extend(['--rotation-step', '{args.rotation_step}'])
if {args.keep_temp}:
    sys.argv.append('--keep-temp')

from automated_glycosylation.cli import run_glyco_prep
run_glyco_prep()
"""]
    
    result = subprocess.run(step1_cmd)
    if result.returncode != 0:
        print("Error in Step 1", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Get the output PDB from step 1
    step1_pdb = step1_dir / "PDB_PROTEIN_GLYCOSYLATED" / "protein_glycosylated_renumbered.pdb"
    
    if not step1_pdb.exists():
        print(f"Error: Step 1 output not found: {step1_pdb}", file=sys.stderr)
        sys.exit(1)
    
    # Step 2
    print("\n" + "="*60)
    print("STEP 2: PARAMETRIZATION")
    print("="*60)
    
    step2_cmd = [sys.executable, "-c", f"""
import sys
sys.argv = ['glyco-param', '-i', '{step1_pdb}', '-o', '{step2_dir}', '--n-cpus', '{args.n_cpus}']
if {args.download_charmm}:
    sys.argv.append('--download-charmm')
if '{args.charmm_url}' and '{args.charmm_url}' != 'None':
    sys.argv.extend(['--charmm-url', '{args.charmm_url}'])
if {args.keep_temp}:
    sys.argv.append('--keep-intermediate')

from automated_glycosylation.cli import run_glyco_param
run_glyco_param()
"""]
    
    result = subprocess.run(step2_cmd)
    if result.returncode != 0:
        print("Error in Step 2", file=sys.stderr)
        sys.exit(result.returncode)
    
    # Get the output PDB from step 2
    step2_pdb = step2_dir / "VALENCE_GLYCAN_VARIANTS" / "glycoprotein_final_valence_corrected_variants.pdb"
    
    if not step2_pdb.exists():
        # Try alternative location
        step2_pdb = step2_dir / "PDB_GLYCOPROTEIN" / "glycoprotein_final_valence_corrected.pdb"
    
    if not step2_pdb.exists():
        print(f"Error: Step 2 output not found", file=sys.stderr)
        sys.exit(1)
    
    # Step 3
    print("\n" + "="*60)
    print("STEP 3: CARBOHYDRATE ORIENTATION")
    print("="*60)
    
    step3_cmd = [sys.executable, "-c", f"""
import sys
sys.argv = ['glyco-orient', '-i', '{step2_pdb}', '-o', '{step3_dir}',
           '--theta-step', '{args.theta_step}',
           '--n-steps', '{args.n_steps}',
           '--max-cycles', '{args.max_cycles}',
           '--radius', '{args.mcmc_radius}',
           '--use-coulomb', '{args.use_coulomb}',
           '--n-workers', '{args.n_workers}']
if {args.save_before_after}:
    sys.argv.append('--save-before-after')
if {args.save_individual_glycans}:
    sys.argv.append('--save-individual-glycans')
if {args.verbose}:
    sys.argv.append('--verbose')

from automated_glycosylation.cli import run_glyco_orient
run_glyco_orient()
"""]
    
    result = subprocess.run(step3_cmd)
    if result.returncode != 0:
        print("Error in Step 3", file=sys.stderr)
        sys.exit(result.returncode)
    
    print("\n" + "="*60)
    print("COMPLETE PIPELINE FINISHED SUCCESSFULLY!")
    print(f"All results saved to: {args.output_dir}")
    print("="*60)

def main():
    """Main entry point for CLI."""
    if len(sys.argv) < 2:
        print("Usage: glyco-{prep,param,orient,all} [arguments]")
        print("\nAvailable commands:")
        print("  glyco-prep    - Run glycosylation preparation (Step 1)")
        print("  glyco-param   - Run parametrization (Step 2)")
        print("  glyco-orient  - Run carbohydrate orientation (Step 3)")
        print("  glyco-all     - Run complete pipeline (Steps 1-3)")
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
        print(f"Unknown command: {command}", file=sys.stderr)
        print("Available commands: prep, param, orient, all", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
