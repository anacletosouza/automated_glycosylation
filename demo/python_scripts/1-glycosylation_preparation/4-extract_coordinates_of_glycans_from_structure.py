"""
Script to process glycosylated PDB files:
1. Remove hydrogen atoms from the PDB file
2. Extract glycan coordinates and save them to separate PDB files
3. Create a dictionary of glycans
4. Generate a TSV file with glycan-protein linkage information

Usage:
------
python glycan_processing.py --input_pdb input.pdb --output_noH output_noH.pdb --output_dir output_directory

Arguments:
----------
--input_pdb      : Path to input PDB file (glycosylated)
--output_noH     : Path to output PDB file without hydrogens
--output_dir     : Directory to save glycan PDB files and linkage TSV
"""

import os
import numpy as np
from pathlib import Path
import shutil
from collections import defaultdict
import argparse

def remove_hydrogens(input_pdb, output_pdb):
    """
    Remove hydrogen atoms from a PDB file.
    
    Parameters:
    -----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file (without hydrogens)
    """
    print(f"Removing hydrogens from {input_pdb}...")

    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                if element == 'H' or atom_name.startswith('H'):
                    continue
                f_out.write(line)
            elif line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                if element == 'H' or atom_name.startswith('H'):
                    continue
                f_out.write(line)
            else:
                f_out.write(line)

    print(f"Saved PDB without hydrogens to {output_pdb}")

def parse_pdb(pdb_file):
    """
    Parse PDB file and separate ATOM (protein) and HETATM (glycan) records.
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    
    Returns:
    --------
    tuple : (protein_lines, glycan_lines)
        Lists of lines for protein and glycan atoms
    """
    protein_lines = []
    glycan_lines = []

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                protein_lines.append(line)
            elif line.startswith('HETATM'):
                glycan_lines.append(line)
            elif line.startswith('TER'):
                protein_lines.append(line)

    return protein_lines, glycan_lines

def extract_glycans(glycan_lines):
    """
    Extract individual glycans from HETATM records.
    
    Parameters:
    -----------
    glycan_lines : list
        List of HETATM lines from PDB file
    
    Returns:
    --------
    dict : Dictionary with glycan IDs as keys and lists of lines as values
    """
    chains = defaultdict(list)
    for line in glycan_lines:
        chain_id = line[21:22].strip()
        chains[chain_id].append(line)

    glycans = {}
    glycan_counter = 1
    glycan_residues = {'NDG', 'FCA', 'NAG', 'BMA', 'MAN', 'GAL', 'SIA', 'A2G'}

    for chain_id, chain_lines in chains.items():
        print(f"  Processing chain {chain_id} with {len(chain_lines)} glycan lines")

        residue_groups = defaultdict(list)
        for line in chain_lines:
            res_num = int(line[22:26].strip())
            residue_groups[res_num].append(line)

        sorted_residues = sorted(residue_groups.keys())

        start_points = []
        for i, res_num in enumerate(sorted_residues):
            first_line = residue_groups[res_num][0]
            res_name = first_line[17:20].strip()
            if res_name in ['NDG', 'A2G']:
                start_points.append(res_num)

        print(f"    Found {len(start_points)} potential glycan start points: {start_points}")

        for start_idx, start_res in enumerate(start_points):
            end_res = start_res
            if start_idx < len(start_points) - 1:
                next_start = start_points[start_idx + 1]
                end_res = next_start - 1
            else:
                end_res = sorted_residues[-1]

            glycan_residue_numbers = [r for r in sorted_residues if start_res <= r <= end_res]

            glycan_lines_list = []
            for res_num in glycan_residue_numbers:
                glycan_lines_list.extend(residue_groups[res_num])

            glycan_id = f"{chain_id}_{glycan_counter}"
            glycan_counter += 1

            glycans[glycan_id] = {
                'lines': glycan_lines_list,
                'chain': chain_id,
                'residue_numbers': glycan_residue_numbers,
                'start_residue': min(glycan_residue_numbers),
                'end_residue': max(glycan_residue_numbers),
            }

            print(f"    Glycan {glycan_id}: residues {glycan_residue_numbers[0]}-{glycan_residue_numbers[-1]}, "
                  f"atoms: {len(glycan_lines_list)}")

            sorted_residues = [r for r in sorted_residues if r > end_res]
            if not sorted_residues:
                break

    return glycans

def save_glycan_pdbs(glycans, output_dir):
    """
    Save each glycan to a separate PDB file and extract accurate residue info.
    
    Parameters:
    -----------
    glycans : dict
        Dictionary of glycans from extract_glycans()
    output_dir : str
        Directory to save PDB files
    
    Returns:
    --------
    Path : Path to PDB directory
    """
    pdb_dir = Path(output_dir) / "PDB"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving glycan PDB files to {pdb_dir}...")

    filenames_used = {}

    for glycan_id, glycan_data in glycans.items():
        base_filename = f"{glycan_id}.pdb"
        output_path = pdb_dir / base_filename

        with open(output_path, 'w') as f:
            for line in glycan_data['lines']:
                f.write(line)
            f.write("TER\n")

        residue_info = {}
        with open(output_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    if res_num not in residue_info:
                        residue_info[res_num] = res_name
        
        sorted_residue_numbers = sorted(residue_info.keys())
        
        sequence_parts = []
        residue_counter = 1
        for res_num in sorted_residue_numbers:
            res_name = residue_info[res_num]
            sequence_parts.append(f"{res_name}{residue_counter}")
            residue_counter += 1
        
        sequence_poly = "_".join(sequence_parts)
        simple_sequence_parts = [residue_info[res_num] for res_num in sorted_residue_numbers]
        simple_sequence_poly = "_".join(simple_sequence_parts)
        actual_residue_count = len(sorted_residue_numbers)

        glycan_data['residue_sequence'] = sequence_poly
        glycan_data['simple_residue_sequence'] = simple_sequence_poly
        glycan_data['residues'] = simple_sequence_parts
        glycan_data['unit_number'] = actual_residue_count
        glycan_data['residue_info'] = residue_info

        print(f"  Saved {glycan_id}: {base_filename}")
        print(f"    Full sequence: {sequence_poly}")
        print(f"    Simple sequence: {simple_sequence_poly}")
        print(f"    Actual residue count: {actual_residue_count}")

    return pdb_dir

def calculate_atom_distance(coord1, coord2):
    """
    Calculate Euclidean distance between two 3D coordinates.
    
    Parameters:
    -----------
    coord1, coord2 : tuple or list
        (x, y, z) coordinates
    
    Returns:
    --------
    float : Euclidean distance
    """
    return np.sqrt((coord1[0] - coord2[0])**2 +
                   (coord1[1] - coord2[1])**2 +
                   (coord1[2] - coord2[2])**2)

def extract_coordinates_from_line(line):
    """
    Extract coordinates from a PDB line.
    
    Parameters:
    -----------
    line : str
        PDB line
    
    Returns:
    --------
    tuple : (x, y, z) coordinates
    """
    try:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        return (x, y, z)
    except (ValueError, IndexError):
        return None

def find_glycan_linkages(protein_lines, glycans):
    """
    Find linkages between glycans and protein residues.
    
    Parameters:
    -----------
    protein_lines : list
        List of ATOM lines from protein
    glycans : dict
        Dictionary of glycans from extract_glycans()
    
    Returns:
    --------
    list : List of dictionaries with linkage information
    """
    linkages = []
    protein_atoms = {}

    for line in protein_lines:
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            res_num = int(line[22:26].strip())
            atom_key = f"{chain_id}_{res_num}_{res_name}_{atom_name}"

            coords = extract_coordinates_from_line(line)
            if coords:
                protein_atoms[atom_key] = {
                    'coordinates': coords,
                    'res_name': res_name,
                    'chain': chain_id,
                    'res_num': res_num,
                    'atom_name': atom_name,
                    'line': line
                }

    for glycan_id, glycan_data in glycans.items():
        glycan_chain = glycan_data['chain']
        glycan_lines = glycan_data['lines']
        c1_coords = None
        c1_line = None

        for line in glycan_lines:
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            res_num = int(line[22:26].strip())
            if res_num == glycan_data['start_residue'] and atom_name == 'C1':
                c1_coords = extract_coordinates_from_line(line)
                c1_line = line
                break

        if not c1_coords:
            first_residue_type = glycan_data['residues'][0]
            for line in glycan_lines:
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                if res_name == first_residue_type and atom_name == 'C1':
                    c1_coords = extract_coordinates_from_line(line)
                    c1_line = line
                    break

        if not c1_coords:
            print(f"Warning: No C1 atom found for glycan {glycan_id}")
            continue

        linking_type = "N-linked" if glycan_data['residues'][0] == 'NDG' else "O-linked"
        target_atoms = []

        if linking_type == "N-linked":
            for atom_key, atom_data in protein_atoms.items():
                if atom_data['res_name'] == 'ASN' and atom_data['atom_name'] == 'ND2':
                    target_atoms.append(atom_data)
        else:
            for atom_key, atom_data in protein_atoms.items():
                if (atom_data['res_name'] == 'SER' and atom_data['atom_name'] == 'OG') or \
                   (atom_data['res_name'] == 'THR' and atom_data['atom_name'] == 'OG1'):
                    target_atoms.append(atom_data)

        min_distance = float('inf')
        closest_atom = None

        for atom_data in target_atoms:
            distance = calculate_atom_distance(c1_coords, atom_data['coordinates'])
            if distance < min_distance:
                min_distance = distance
                closest_atom = atom_data

        if closest_atom and min_distance < 5.0:
            glycan_res_num = int(c1_line[22:26].strip())
            glycan_atom_name = c1_line[12:16].strip()
            glycan_res_name = c1_line[17:20].strip()

            linkage_info = {
                'site_protein_residue': f"{closest_atom['res_name']}{closest_atom['res_num']}",
                'glycan_binding': glycan_id,
                'sequence_poly': glycan_data['residue_sequence'],
                'simple_sequence_poly': glycan_data['simple_residue_sequence'],
                'unit_number': glycan_data['unit_number'],
                'protein_residue_number': closest_atom['res_num'],
                'protein_chain': closest_atom['chain'],
                'protein_atom': closest_atom['atom_name'],
                'glycan_atom': f"{glycan_res_name}:{glycan_atom_name}",
                'distance': round(min_distance, 3),
                'linking_type': linking_type
            }
            linkages.append(linkage_info)
            print(f"  Found linkage: {linkage_info['site_protein_residue']} "
                  f"-> {glycan_id} ({linking_type}, distance: {min_distance:.3f} Å)")
        else:
            print(f"Warning: No suitable linkage found for glycan {glycan_id} "
                  f"(closest distance: {min_distance:.3f} Å)")

    return linkages

def save_linkages_tsv(linkages, output_dir):
    """
    Save linkage information to a TSV file.
    
    Parameters:
    -----------
    linkages : list
        List of linkage dictionaries from find_glycan_linkages()
    output_dir : str
        Directory to save TSV file
    """
    output_path = Path(output_dir) / "topol_carb.tsv"
    headers = [
        'site_protein_residue',
        'glycan_binding',
        'sequence_poly',
        'unit_number',
        'protein_residue_number',
        'protein_chain',
        'protein_atom',
        'glycan_atom',
        'distance',
        'linking_type'
    ]

    print(f"Saving linkage information to {output_path}...")

    with open(output_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for linkage in linkages:
            row = [
                linkage['site_protein_residue'],
                linkage['glycan_binding'],
                linkage['sequence_poly'],
                str(linkage['unit_number']),
                str(linkage['protein_residue_number']),
                linkage['protein_chain'],
                linkage['protein_atom'],
                linkage['glycan_atom'],
                str(linkage['distance']),
                linkage['linking_type']
            ]
            f.write('\t'.join(row) + '\n')

    print(f"  Saved {len(linkages)} linkages to TSV file")

def main():
    """Main function to execute all steps."""
    parser = argparse.ArgumentParser(description="Process glycosylated PDB files and extract glycan info")
    parser.add_argument("--input_pdb", required=True, help="Path to input glycosylated PDB file")
    parser.add_argument("--output_noH", required=True, help="Path to output PDB file without hydrogens")
    parser.add_argument("--output_dir", required=True, help="Directory to save glycan PDB files and TSV")

    args = parser.parse_args()

    input_pdb = args.input_pdb
    output_pdb_noH = args.output_noH
    topol_carb_dir = args.output_dir

    print("=" * 70)
    print("GLYCAN PROCESSING SCRIPT")
    print("=" * 70)

    print("\n1. REMOVING HYDROGEN ATOMS")
    print("-" * 40)
    remove_hydrogens(input_pdb, output_pdb_noH)

    print("\n2. PARSING PDB FILE")
    print("-" * 40)
    protein_lines, glycan_lines = parse_pdb(input_pdb)
    print(f"  Found {len(protein_lines)} protein lines")
    print(f"  Found {len(glycan_lines)} glycan lines")

    print("\n3. EXTRACTING GLYCANS")
    print("-" * 40)
    glycans = extract_glycans(glycan_lines)
    print(f"  Found {len(glycans)} individual glycans")

    print("\n4. SAVING GLYCAN PDB FILES AND EXTRACTING RESIDUE INFORMATION")
    print("-" * 40)
    pdb_dir = save_glycan_pdbs(glycans, topol_carb_dir)

    print("\n5. FINAL GLYCAN INFORMATION")
    print("-" * 40)
    for i, (glycan_id, glycan_data) in enumerate(glycans.items(), 1):
        print(f"  Glycan {i}: {glycan_id}")
        print(f"    Chain: {glycan_data['chain']}")
        print(f"    Full sequence: {glycan_data['residue_sequence']}")
        print(f"    Simple sequence: {glycan_data['simple_residue_sequence']}")
        print(f"    Unit number (residue count): {glycan_data['unit_number']}")
        print(f"    Number of atoms: {len(glycan_data['lines'])}")
        print(f"    Residue range: {glycan_data['start_residue']} to {glycan_data['end_residue']}")

    print("\n6. FINDING GLYCAN-PROTEIN LINKAGES")
    print("-" * 40)
    linkages = find_glycan_linkages(protein_lines, glycans)

    print("\n7. SAVING LINKAGE INFORMATION")
    print("-" * 40)
    save_linkages_tsv(linkages, topol_carb_dir)

    print("\n" + "=" * 70)
    print("PROCESSING COMPLETE!")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  - Input PDB: {input_pdb}")
    print(f"  - PDB without hydrogens: {output_pdb_noH}")
    print(f"  - Number of glycans extracted: {len(glycans)}")
    print(f"  - Glycan PDB files saved to: {pdb_dir}")
    print(f"  - Number of linkages found: {len(linkages)}")
    print(f"  - Linkage TSV file: {Path(topol_carb_dir) / 'topol_carb.tsv'}")

    return glycans, linkages

if __name__ == "__main__":
    glycans, linkages = main()

