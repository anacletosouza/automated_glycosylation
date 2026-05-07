#!/usr/bin/env python3
"""
Script to process glycan data from multiple sources and create organized directories
with JSON files for carbohydrate parameterization.

Author: Anacleto
Date: 2026-01-08
"""

import os
import json
import pandas as pd
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re
import argparse

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
    
    # Input files
    topol_carb_path = base_dir / "TO_TOP" / "topol_carb.tsv"
    caselino_table_path = base_dir / "TSV" / "caselino_2020_tables_glycosylator.tsv"
    pdb_dir = base_dir / "TO_TOP" / "PDB"
    
    # Log file for warnings
    warning_log_path = output_base_dir / "json_warnings.txt"
    
    # Create output directory if it doesn't exist
    output_base_dir.mkdir(parents=True, exist_ok=True)
    
    # Read data files
    print("Reading data files...")
    try:
        topol_carb_df = pd.read_csv(topol_carb_path, sep='\t')
        caselino_df = pd.read_csv(caselino_table_path, sep='\t')
    except Exception as e:
        print(f"Error reading input files: {e}")
        return
    
    if len(topol_carb_df) != len(caselino_df):
        print(f"Warning: Dataframes have different lengths!")
        print(f"topol_carb.tsv: {len(topol_carb_df)} rows")
        print(f"caselino_2020_tables_glycosylator.tsv: {len(caselino_df)} rows")
    
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
        
        if idx < len(caselino_df):
            iupac_sequence = caselino_df.loc[idx, 'iupac_glycosylator']
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
    parser.add_argument("--base_dir", type=str, default="/grain/anacleto/project/AA_simulations/with_membranes/glycosylation_method/Delta",
                        help="Base directory containing input files and PDB directory")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Directory where JSON files and logs will be saved (default: under base_dir/python_scripts/...)")
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    
    if args.output_dir:
        output_base_dir = Path(args.output_dir)
    else:
        output_base_dir = base_dir / "python_scripts" / "2-parametrization_scripts" / "JSON"
    
    required_files = [
        base_dir / "TO_TOP" / "topol_carb.tsv",
        base_dir / "TSV" / "caselino_2020_tables_glycosylator.tsv",
        base_dir / "TO_TOP" / "PDB"
    ]
    for req_file in required_files:
        if not req_file.exists():
            print(f"Error: Required file/directory not found: {req_file}")
            return
    
    print("Glycan Processing Script")
    print("=" * 50)
    try:
        process_glycans(base_dir, output_base_dir)
    except Exception as e:
        print(f"\nError during processing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()

