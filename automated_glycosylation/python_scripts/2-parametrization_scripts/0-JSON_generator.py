#!/usr/bin/env python3
"""
Script to process glycan data from multiple sources and create organized directories
with JSON files for carbohydrate parameterization.

Author: Anacleto
Date: 2026-01-08
Modified: 2026-05-07 - Generic file names
"""

import os
import json
import pandas as pd
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re
import argparse
import glob

# Define the mapping from PDB residue names to CHARMM names
CARB_PDB_CHARMM_DICT = {
    "NDG": "BGLCNA",
    "FCA": "AFUC", 
    "NAG": "BGLCNA",
    "BMA": "BMAN",
    "MAN": "AMAN",
    "A2G": "AGALNA",
    "GAL": "BLGAL",
    "SIA": "ANE5AC"
    # Add more mappings as needed
}

def find_file(base_dir: Path, pattern: str) -> Optional[Path]:
    """Find a file matching pattern in directory."""
    matches = list(base_dir.glob(pattern))
    if matches:
        return matches[0]
    return None

def extract_residue_names(sequence_poly: str) -> List[str]:
    if pd.isna(sequence_poly) or not sequence_poly:
        return []
    
    residues = []
    seen = set()
    parts = sequence_poly.split('_')
    for part in parts:
        residue_name = re.sub(r'\d+$', '', part)
        if residue_name and residue_name not in seen:
            residues.append(residue_name)
            seen.add(residue_name)
    return residues

def create_residue_list(residue_names: List[str]) -> List[Dict[str, str]]:
    residues = []
    for res_name in residue_names:
        charmm_name = CARB_PDB_CHARMM_DICT.get(res_name, "*****")
        residues.append({
            "pdb_name": res_name,
            "charmm_name": charmm_name
        })
    return residues

def process_glycans(base_dir: Path, output_base_dir: Path):
    """Main function to process all glycan data."""
    
    # Find input files dynamically (generic names)
    topol_carb_path = find_file(base_dir / "TO_TOP", "topol_carb*.tsv")
    if not topol_carb_path:
        topol_carb_path = base_dir / "TO_TOP" / "topol_carb.tsv"
    
    # Find the glycosylator TSV file (generic pattern)
    tsv_dir = base_dir / "TSV"
    glycosylator_path = find_file(tsv_dir, "*_glycosylator.tsv")
    if not glycosylator_path:
        # Try alternative patterns
        glycosylator_path = find_file(tsv_dir, "*.tsv")
        if glycosylator_path:
            # Filter out corrected files
            if "_corrected" in str(glycosylator_path):
                glycosylator_path = find_file(tsv_dir, "*_glycosylator.tsv")
    
    pdb_dir = base_dir / "TO_TOP" / "PDB"
    
    # Log file for warnings
    warning_log_path = output_base_dir / "json_warnings.txt"
    
    # Create output directory if it doesn't exist
    output_base_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if required files exist
    if not topol_carb_path.exists():
        print(f"Error: topol_carb file not found in {base_dir / 'TO_TOP'}")
        print(f"Looking for pattern: topol_carb*.tsv")
        return
    
    if not glycosylator_path or not glycosylator_path.exists():
        print(f"Error: Glycosylator TSV file not found in {tsv_dir}")
        print(f"Looking for pattern: *_glycosylator.tsv")
        return
    
    if not pdb_dir.exists():
        print(f"Error: PDB directory not found: {pdb_dir}")
        return
    
    print(f"Using topol_carb file: {topol_carb_path}")
    print(f"Using glycosylator file: {glycosylator_path}")
    print(f"Using PDB directory: {pdb_dir}")
    
    # Read data files
    print("\nReading data files...")
    try:
        topol_carb_df = pd.read_csv(topol_carb_path, sep='\t')
        glycosylator_df = pd.read_csv(glycosylator_path, sep='\t')
    except Exception as e:
        print(f"Error reading input files: {e}")
        return
    
    if len(topol_carb_df) != len(glycosylator_df):
        print(f"Warning: Dataframes have different lengths!")
        print(f"topol_carb.tsv: {len(topol_carb_df)} rows")
        print(f"glycosylator file: {len(glycosylator_df)} rows")
    
    warnings = []
    
    print(f"\nProcessing {len(topol_carb_df)} glycans...")
    for idx, row in topol_carb_df.iterrows():
        glycan_binding = row['glycan_binding']
        pdb_filename = f"{glycan_binding}.pdb"
        pdb_path = pdb_dir / pdb_filename
        
        if not pdb_path.exists():
            print(f"Warning: PDB file not found for {glycan_binding}: {pdb_filename}")
            warnings.append(f"Missing PDB: {pdb_filename}")
            continue
        
        if idx < len(glycosylator_df):
            iupac_sequence = glycosylator_df.loc[idx, 'iupac_glycosylator']
            if pd.isna(iupac_sequence):
                iupac_sequence = ""
        else:
            iupac_sequence = ""
            warnings.append(f"No IUPAC sequence for {glycan_binding} (index out of range)")
        
        sequence_poly = row['sequence_poly']
        residue_names = extract_residue_names(sequence_poly)
        
        for res_name in residue_names:
            if res_name not in CARB_PDB_CHARMM_DICT:
                warnings.append(f"Missing CHARMM mapping for {res_name} in {glycan_binding}")
        
        output_dir = output_base_dir / glycan_binding
        output_dir.mkdir(exist_ok=True)
        shutil.copy2(pdb_path, output_dir / pdb_filename)
        
        residues = create_residue_list(residue_names)
        json_data = {
            "sequence": str(iupac_sequence) if not pd.isna(iupac_sequence) else "",
            "residues": residues
        }
        
        with open(output_dir / f"{glycan_binding}.json", 'w') as json_file:
            json.dump(json_data, json_file, indent=2)
        
        print(f"Created: {glycan_binding}/")
    
    if warnings:
        print(f"\nFound {len(warnings)} warnings. Writing to {warning_log_path}")
        with open(warning_log_path, 'w') as warn_file:
            warn_file.write("JSON Warnings Log\n")
            warn_file.write("=" * 50 + "\n\n")
            for i, warning in enumerate(warnings, 1):
                warn_file.write(f"{i}. {warning}\n")
        
        unique_warnings = set(warnings)
        print(f"\nUnique warnings:")
        for warning in unique_warnings:
            print(f"  - {warning}")
    
    print("\n" + "=" * 50)
    print("PROCESSING SUMMARY")
    print("=" * 50)
    
    created_dirs = [d for d in output_base_dir.iterdir() if d.is_dir()]
    print(f"Created directories: {len(created_dirs)}")
    json_files = list(output_base_dir.rglob("*.json"))
    print(f"Created JSON files: {len(json_files)}")
    pdb_files = list(output_base_dir.rglob("*.pdb"))
    print(f"Copied PDB files: {len(pdb_files)}")
    
    if warnings:
        print(f"Warnings logged: {len(warnings)} (see {warning_log_path})")
    
    print(f"\nOutput directory: {output_base_dir}")
    print("Processing complete!")

def main():
    parser = argparse.ArgumentParser(description="Process glycan data and generate JSON files.")
    parser.add_argument("--base_dir", type=str, required=True,
                        help="Base directory containing input files (output from previous pipeline)")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory where JSON files and logs will be saved")
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    output_base_dir = Path(args.output_dir)
    
    # Check if base directory exists
    if not base_dir.exists():
        print(f"Error: Base directory not found: {base_dir}")
        return
    
    # Create output directory
    output_base_dir.mkdir(parents=True, exist_ok=True)
    
    print("Glycan Processing Script")
    print("=" * 50)
    print(f"Base directory: {base_dir}")
    print(f"Output directory: {output_base_dir}")
    print("=" * 50)
    
    try:
        process_glycans(base_dir, output_base_dir)
    except Exception as e:
        print(f"\nError during processing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
