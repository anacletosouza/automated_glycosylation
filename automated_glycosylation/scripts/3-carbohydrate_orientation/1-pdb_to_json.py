"""
Script to process glycosylated PDB files and extract complete information into JSON:
1. Remove hydrogen atoms from the PDB file
2. Extract complete protein and glycan information
3. Find glycan-protein linkages
4. Save all data to a JSON file that can reconstruct the original PDB

Usage:
------
python glycan_processing.py --input_pdb input.pdb --output_json output_data.json

Arguments:
----------
--input_pdb      : Path to input PDB file (glycosylated)
--output_json    : Path to output JSON file with all data
"""

import os
import numpy as np
from pathlib import Path
import json
import argparse
import re

def get_unique_hetatm_residues(pdb_file):
    """
    Get unique residue names from HETATM records in a PDB file.
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    
    Returns:
    --------
    list : Sorted list of unique residue names
    """
    residues = set()
    with open(pdb_file) as f:
        for line in f:
            if line.upper().startswith("HETATM"):
                # Extract residue name (columns 17-20)
                residue_name = line[17:20].strip()
                residues.add(residue_name)
    return sorted(residues)

def safe_int_parse(value, default=None):
    """
    Safely parse integer values, handling non-numeric strings.
    
    Parameters:
    -----------
    value : str
        String to parse as integer
    default : any
        Default value to return if parsing fails
    
    Returns:
    --------
    int or default : Parsed integer or default value
    """
    if not value:
        return default
    
    # Remove any non-digit characters
    cleaned = re.sub(r'[^0-9]', '', value)
    if cleaned:
        try:
            return int(cleaned)
        except ValueError:
            return default
    return default

def safe_float_parse(value, default=None):
    """
    Safely parse float values, handling non-numeric strings.
    
    Parameters:
    -----------
    value : str
        String to parse as float
    default : any
        Default value to return if parsing fails
    
    Returns:
    --------
    float or default : Parsed float or default value
    """
    if not value:
        return default
    
    # Handle scientific notation and regular floats
    try:
        return float(value)
    except ValueError:
        # Try to extract any numeric part
        match = re.search(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', value)
        if match:
            try:
                return float(match.group())
            except ValueError:
                return default
        return default

def parse_atom_line(line):
    """
    Parse a PDB ATOM/HETATM line into a dictionary with all information.
    
    Parameters:
    -----------
    line : str
        PDB line starting with ATOM or HETATM
    
    Returns:
    --------
    dict : Dictionary with all atom information
    """
    record = line[0:6].strip()
    
    # Extract fields with safe parsing
    atom_number = safe_int_parse(line[6:11].strip())
    atom_name = line[12:16].strip()
    alt_loc = line[16:17].strip()
    residue_name = line[17:20].strip()
    chain_id = line[21:22].strip()
    
    # Handle residue number - might contain chain ID in some non-standard formats
    residue_number_str = line[22:26].strip()
    residue_number = safe_int_parse(residue_number_str)
    
    # If residue_number couldn't be parsed, try to extract from a different format
    if residue_number is None and residue_number_str:
        # Some PDBs have format like "A 17" or "A17"
        # Try to extract numbers from the string
        for part in residue_number_str.split():
            temp = safe_int_parse(part)
            if temp is not None:
                residue_number = temp
                break
    
    icode = line[26:27].strip()
    x = safe_float_parse(line[30:38].strip(), 0.0)
    y = safe_float_parse(line[38:46].strip(), 0.0)
    z = safe_float_parse(line[46:54].strip(), 0.0)
    occupancy = safe_float_parse(line[54:60].strip(), 1.00)
    temp_factor = safe_float_parse(line[60:66].strip(), 0.00)
    element = line[76:78].strip()
    charge = line[78:80].strip()
    
    atom_data = {
        'record': record,
        'atom_number': atom_number,
        'atom_name': atom_name,
        'alt_loc': alt_loc,
        'residue_name': residue_name,
        'chain_id': chain_id,
        'residue_number': residue_number,
        'icode': icode,
        'x': x,
        'y': y,
        'z': z,
        'occupancy': occupancy,
        'temp_factor': temp_factor,
        'element': element,
        'charge': charge,
        'original_line': line.rstrip()
    }
    
    return atom_data

def is_glycan_residue(residue_name):
    """
    Determine if a residue is a glycan residue based on common patterns.
    
    Parameters:
    -----------
    residue_name : str
        Name of the residue
    
    Returns:
    --------
    bool : True if it's a glycan residue
    """
    # Remove spaces
    residue = residue_name.strip()
    
    # Check if it starts with common carbohydrate patterns
    if len(residue) >= 2:
        prefix = residue[:2]
        # Common patterns for glycan residues
        if prefix in ['ND', 'NA', 'FC', 'BM', 'MA', 'A2', 'GA', 'SI', 'GL', 'XY']:
            return True
    
    # Check common full names
    common_glycans = {'FUC', 'MAN', 'GAL', 'GLC', 'XYS', 'SIA'}
    if residue in common_glycans:
        return True
    
    return False

def parse_pdb_complete(pdb_file):
    """
    Parse complete PDB file, separating protein atoms, glycan atoms, and other records.
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    
    Returns:
    --------
    dict : Dictionary with all PDB information
    """
    # Get unique HETATM residues from the PDB file
    unique_hetatm_residues = get_unique_hetatm_residues(pdb_file)
    
    # Automatically detect glycan residues based on patterns
    glycan_residues = set()
    for residue in unique_hetatm_residues:
        if is_glycan_residue(residue):
            glycan_residues.add(residue)
    
    print(f"  Unique HETATM residues found: {sorted(unique_hetatm_residues)}")
    print(f"  Automatically detected glycan residues: {sorted(glycan_residues)}")
    
    protein_atoms = []
    glycan_atoms = []
    other_lines = []
    current_model = 1
    
    with open(pdb_file, 'r') as f:
        line_num = 0
        for line in f:
            line_num += 1
            stripped_line = line.rstrip()
            
            if line.startswith('MODEL'):
                current_model = int(line[10:14].strip()) if line[10:14].strip() else 1
            
            elif line.startswith('ATOM'):
                try:
                    atom_data = parse_atom_line(line)
                    atom_data['model'] = current_model
                    protein_atoms.append(atom_data)
                except Exception as e:
                    print(f"  Warning: Error parsing ATOM line {line_num}: {e}")
                    print(f"    Line: {stripped_line}")
                    other_lines.append(stripped_line)
            
            elif line.startswith('HETATM'):
                try:
                    atom_data = parse_atom_line(line)
                    atom_data['model'] = current_model
                    
                    # Use automatically determined glycan residues
                    if atom_data['residue_name'] in glycan_residues:
                        glycan_atoms.append(atom_data)
                    else:
                        other_lines.append(stripped_line)
                except Exception as e:
                    print(f"  Warning: Error parsing HETATM line {line_num}: {e}")
                    print(f"    Line: {stripped_line}")
                    other_lines.append(stripped_line)
            
            elif line.startswith('TER'):
                other_lines.append(stripped_line)
            
            elif line.startswith('ENDMDL'):
                other_lines.append(stripped_line)
            
            elif line.startswith('END'):
                other_lines.append(stripped_line)
            
            elif line.startswith('CONECT'):
                other_lines.append(stripped_line)
            
            elif line.strip():
                other_lines.append(stripped_line)
    
    return {
        'protein_atoms': protein_atoms,
        'glycan_atoms': glycan_atoms,
        'other_lines': other_lines,
        'detected_glycan_residues': sorted(glycan_residues),
        'all_hetatm_residues': sorted(unique_hetatm_residues)
    }

def extract_glycans_from_atoms(glycan_atoms):
    """
    Extract individual glycans preserving the exact order
    they appear in the PDB file.
    """

    residue_map = {}
    residue_order = []  # <<< ordem real do PDB

    for atom in glycan_atoms:
        key = (atom['chain_id'], atom['residue_number'])

        if key not in residue_map:
            residue_map[key] = {
                'residue_name': atom['residue_name'],
                'atoms': [],
                'residue_number': atom['residue_number'],
                'chain': atom['chain_id']
            }
            residue_order.append(key)  # <<< salva ordem de aparecimento

        residue_map[key]['atoms'].append(atom)

    glycans = {}
    current_glycan = []
    current_glycan_id = 1  # numeração GLOBAL

    for key in residue_order:
        chain_id, res_num = key
        residue_name = residue_map[key]['residue_name']

        start_new = False

        # regra original preservada
        if residue_name.startswith(('ND', 'A2')):
            start_new = True
        elif not current_glycan:
            start_new = True

        if start_new and current_glycan:
            first_chain = current_glycan[0][0] or "_"
            glycan_id = f"{first_chain}_{current_glycan_id}"

            glycans[glycan_id] = create_glycan_data(
                glycan_id, current_glycan, residue_map
            )

            current_glycan_id += 1
            current_glycan = []

        current_glycan.append(key)

    # último glicano
    if current_glycan:
        first_chain = current_glycan[0][0] or "_"
        glycan_id = f"{first_chain}_{current_glycan_id}"

        glycans[glycan_id] = create_glycan_data(
            glycan_id, current_glycan, residue_map
        )

    return glycans


def create_glycan_data(glycan_id, residue_keys, residue_map):
    atoms = []
    residue_numbers = []
    residue_names = []
    chains = set()

    for key in residue_keys:
        info = residue_map[key]
        atoms.extend(info['atoms'])
        residue_numbers.append(info['residue_number'])
        residue_names.append(info['residue_name'])
        chains.add(info['chain'])

    return {
        'id': glycan_id,
        'chain': next(iter(chains)) if chains else '',
        'atoms': atoms,
        'residue_numbers': residue_numbers,
        'start_residue': min(residue_numbers),
        'end_residue': max(residue_numbers),
        'simple_sequence': "_".join(residue_names),
        'full_sequence': "_".join(
            f"{r}{i+1}" for i, r in enumerate(residue_names)
        ),
        'unit_number': len(residue_numbers),
        'residue_info': dict(zip(residue_numbers, residue_names)),
        'residues': residue_names
    }

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

def find_glycan_linkages(protein_atoms, glycans):
    """
    Find linkages between glycans and protein residues with complete coordinate information.
    
    Parameters:
    -----------
    protein_atoms : list
        List of protein atom dictionaries
    glycans : dict
        Dictionary of glycans
    
    Returns:
    --------
    list : List of dictionaries with complete linkage information including coordinates
    """
    linkages = []
    
    # Pre-process protein atoms for faster lookup
    # Create dictionaries for different atom types
    asn_nd2_atoms = []
    ser_og_atoms = []
    thr_og1_atoms = []
    
    for atom in protein_atoms:
        if atom['residue_name'] == 'ASN' and atom['atom_name'] == 'ND2':
            asn_nd2_atoms.append(atom)
        elif atom['residue_name'] == 'SER' and atom['atom_name'] == 'OG':
            ser_og_atoms.append(atom)
        elif atom['residue_name'] == 'THR' and atom['atom_name'] == 'OG1':
            thr_og1_atoms.append(atom)
    
    print(f"  Found {len(asn_nd2_atoms)} ASN-ND2 atoms (N-linked sites)")
    print(f"  Found {len(ser_og_atoms)} SER-OG atoms (O-linked sites)")
    print(f"  Found {len(thr_og1_atoms)} THR-OG1 atoms (O-linked sites)")
    
    processed_glycans = 0
    linked_glycans = set()
    
    for glycan_id, glycan_data in glycans.items():
        processed_glycans += 1
        print(f"  Processing glycan {glycan_id} ({processed_glycans}/{len(glycans)})")
        
        # Find C1 atom in the glycan
        c1_atom = None
        for atom in glycan_data['atoms']:
            # Look for C1 atom - common in carbohydrates
            if 'C1' in atom['atom_name']:
                c1_atom = atom
                break
        
        # If C1 not found, try alternative naming
        if not c1_atom:
            for atom in glycan_data['atoms']:
                if atom['atom_name'] == 'C1' or atom['atom_name'] == 'C1A' or atom['atom_name'] == 'C1B':
                    c1_atom = atom
                    break
        
        # If still not found, use first atom of first residue
        if not c1_atom and glycan_data['atoms']:
            c1_atom = glycan_data['atoms'][0]
            print(f"    Warning: No C1 atom found for glycan {glycan_id}, using first atom instead")
        
        if not c1_atom:
            print(f"    Warning: No atoms found for glycan {glycan_id}")
            continue
        
        c1_coords = (c1_atom['x'], c1_atom['y'], c1_atom['z'])
        
        # Determine expected linkage type based on first residue name
        first_res_name = glycan_data['residue_info'][glycan_data['start_residue']]
        
        # IMPORTANT: A21 and A22 are N-linked glycans (modified N-linked sites)
        # ND* residues are always N-linked
        if first_res_name.startswith('ND') or first_res_name.startswith('A2'):
            expected_type = 'N-linked'
            search_atoms = asn_nd2_atoms
        else:
            # For other residues, could be O-linked or other types
            expected_type = 'O-linked'
            search_atoms = ser_og_atoms + thr_og1_atoms
        
        print(f"    First residue: {first_res_name}, expected type: {expected_type}")
        
        # Find the closest appropriate protein atom
        min_distance = float('inf')
        closest_protein_atom = None
        actual_type = expected_type
        
        for atom in search_atoms:
            atom_coords = (atom['x'], atom['y'], atom['z'])
            distance = calculate_atom_distance(c1_coords, atom_coords)
            if distance < min_distance:
                min_distance = distance
                closest_protein_atom = atom
        
        # Check if we found a linkage
        if closest_protein_atom and min_distance < 2.5:  # 2.5 Å is a good threshold for covalent bonds
            # Verify the linkage type matches
            if expected_type == 'N-linked' and closest_protein_atom['residue_name'] == 'ASN':
                actual_type = 'N-linked'
            elif expected_type == 'O-linked' and closest_protein_atom['residue_name'] in ['SER', 'THR']:
                actual_type = 'O-linked'
            else:
                print(f"    Warning: Type mismatch for glycan {glycan_id}")
                actual_type = 'Unknown'
            
            protein_atom_data = {
                'atom_number': closest_protein_atom['atom_number'],
                'atom_name': closest_protein_atom['atom_name'],
                'residue_name': closest_protein_atom['residue_name'],
                'chain_id': closest_protein_atom['chain_id'],
                'residue_number': closest_protein_atom['residue_number'],
                'x': closest_protein_atom['x'],
                'y': closest_protein_atom['y'],
                'z': closest_protein_atom['z'],
                'element': closest_protein_atom['element'],
                'occupancy': closest_protein_atom['occupancy'],
                'temp_factor': closest_protein_atom['temp_factor'],
                'original_line': closest_protein_atom['original_line']
            }
            
            glycan_atom_data = {
                'atom_number': c1_atom['atom_number'],
                'atom_name': c1_atom['atom_name'],
                'residue_name': c1_atom['residue_name'],
                'chain_id': c1_atom['chain_id'],
                'residue_number': c1_atom['residue_number'],
                'x': c1_atom['x'],
                'y': c1_atom['y'],
                'z': c1_atom['z'],
                'element': c1_atom['element'],
                'occupancy': c1_atom['occupancy'],
                'temp_factor': c1_atom['temp_factor'],
                'original_line': c1_atom['original_line']
            }
            
            linkage_info = {
                'glycan_binding': glycan_id,
                'unit_number': glycan_data['unit_number'],
                'protein_residue_number': closest_protein_atom['residue_number'],
                'protein_chain': closest_protein_atom['chain_id'],
                'protein_atom': closest_protein_atom['atom_name'],
                'glycan_atom': f"{c1_atom['residue_name']}:{c1_atom['atom_name']}",
                'distance': round(min_distance, 3),
                'linking_type': actual_type,
                'protein_atom_complete': protein_atom_data,
                'glycan_atom_complete': glycan_atom_data,
                'site_protein_residue': f"{closest_protein_atom['residue_name']}{closest_protein_atom['residue_number']}",
                'sequence_poly': glycan_data['full_sequence'],
                'simple_sequence_poly': glycan_data['simple_sequence'],
                'glycan_first_residue': first_res_name
            }
            
            linkages.append(linkage_info)
            linked_glycans.add(glycan_id)
            print(f"    Found linkage: {linkage_info['site_protein_residue']} -> {glycan_id} "
                  f"({actual_type}, distance: {min_distance:.3f} Å)")
        else:
            print(f"    No suitable linkage found for glycan {glycan_id} "
                  f"(closest distance: {min_distance:.3f} Å)")
    
    # If we have more glycans than linkages, check for missing linkages
    if len(linkages) < len(glycans):
        print(f"\n  WARNING: Found {len(linkages)} linkages for {len(glycans)} glycans")
        print(f"  Looking for missing linkages...")
        
        # Find glycans without linkages
        missing_glycans = [gid for gid in glycans.keys() if gid not in linked_glycans]
        
        for glycan_id in missing_glycans:
            glycan_data = glycans[glycan_id]
            print(f"  Processing missing glycan: {glycan_id}")
            
            # Try to find C1 atom
            c1_atom = None
            for atom in glycan_data['atoms']:
                if 'C1' in atom['atom_name']:
                    c1_atom = atom
                    break
            
            if not c1_atom and glycan_data['atoms']:
                c1_atom = glycan_data['atoms'][0]
            
            if not c1_atom:
                print(f"    No C1 atom found for missing glycan {glycan_id}")
                continue
            
            c1_coords = (c1_atom['x'], c1_atom['y'], c1_atom['z'])
            
            # Try broader search with all potential linking atoms
            all_potential_atoms = asn_nd2_atoms + ser_og_atoms + thr_og1_atoms
            
            min_distance = float('inf')
            closest_protein_atom = None
            
            for atom in all_potential_atoms:
                atom_coords = (atom['x'], atom['y'], atom['z'])
                distance = calculate_atom_distance(c1_coords, atom_coords)
                if distance < min_distance:
                    min_distance = distance
                    closest_protein_atom = atom
            
            # Use a slightly larger threshold for missing glycans
            if closest_protein_atom and min_distance < 3.0:
                # Determine type
                if closest_protein_atom['residue_name'] == 'ASN' and closest_protein_atom['atom_name'] == 'ND2':
                    linking_type = 'N-linked'
                elif closest_protein_atom['residue_name'] == 'SER' and closest_protein_atom['atom_name'] == 'OG':
                    linking_type = 'O-linked'
                elif closest_protein_atom['residue_name'] == 'THR' and closest_protein_atom['atom_name'] == 'OG1':
                    linking_type = 'O-linked'
                else:
                    linking_type = 'Unknown'
                
                protein_atom_data = {
                    'atom_number': closest_protein_atom['atom_number'],
                    'atom_name': closest_protein_atom['atom_name'],
                    'residue_name': closest_protein_atom['residue_name'],
                    'chain_id': closest_protein_atom['chain_id'],
                    'residue_number': closest_protein_atom['residue_number'],
                    'x': closest_protein_atom['x'],
                    'y': closest_protein_atom['y'],
                    'z': closest_protein_atom['z'],
                    'element': closest_protein_atom['element'],
                    'occupancy': closest_protein_atom['occupancy'],
                    'temp_factor': closest_protein_atom['temp_factor'],
                    'original_line': closest_protein_atom['original_line']
                }
                
                glycan_atom_data = {
                    'atom_number': c1_atom['atom_number'],
                    'atom_name': c1_atom['atom_name'],
                    'residue_name': c1_atom['residue_name'],
                    'chain_id': c1_atom['chain_id'],
                    'residue_number': c1_atom['residue_number'],
                    'x': c1_atom['x'],
                    'y': c1_atom['y'],
                    'z': c1_atom['z'],
                    'element': c1_atom['element'],
                    'occupancy': c1_atom['occupancy'],
                    'temp_factor': c1_atom['temp_factor'],
                    'original_line': c1_atom['original_line']
                }
                
                linkage_info = {
                    'glycan_binding': glycan_id,
                    'unit_number': glycan_data['unit_number'],
                    'protein_residue_number': closest_protein_atom['residue_number'],
                    'protein_chain': closest_protein_atom['chain_id'],
                    'protein_atom': closest_protein_atom['atom_name'],
                    'glycan_atom': f"{c1_atom['residue_name']}:{c1_atom['atom_name']}",
                    'distance': round(min_distance, 3),
                    'linking_type': linking_type,
                    'protein_atom_complete': protein_atom_data,
                    'glycan_atom_complete': glycan_atom_data,
                    'site_protein_residue': f"{closest_protein_atom['residue_name']}{closest_protein_atom['residue_number']}",
                    'sequence_poly': glycan_data['full_sequence'],
                    'simple_sequence_poly': glycan_data['simple_sequence'],
                    'glycan_first_residue': glycan_data['residue_info'][glycan_data['start_residue']],
                    'note': 'Found in secondary search'
                }
                
                linkages.append(linkage_info)
                print(f"    Found missing linkage: {linkage_info['site_protein_residue']} -> {glycan_id} "
                      f"({linking_type}, distance: {min_distance:.3f} Å)")
    
    print(f"\n  Total linkages found: {len(linkages)}")
    print(f"  Total glycans: {len(glycans)}")
    
    return linkages

def create_complete_data_structure(pdb_data, glycans, linkages):
    """
    Create complete data structure with all information.
    
    Parameters:
    -----------
    pdb_data : dict
        Parsed PDB data
    glycans : dict
        Glycan data
    linkages : list
        Linkage information
    
    Returns:
    --------
    dict : Complete data structure
    """
    # Prepare protein data
    protein_data = []
    for atom in pdb_data['protein_atoms']:
        protein_data.append({
            'record': atom['record'],
            'atom_number': atom['atom_number'],
            'atom_name': atom['atom_name'],
            'alt_loc': atom['alt_loc'],
            'residue_name': atom['residue_name'],
            'chain_id': atom['chain_id'],
            'residue_number': atom['residue_number'],
            'icode': atom['icode'],
            'x': atom['x'],
            'y': atom['y'],
            'z': atom['z'],
            'occupancy': atom['occupancy'],
            'temp_factor': atom['temp_factor'],
            'element': atom['element'],
            'charge': atom['charge'],
            'model': atom['model']
        })
    
    # Prepare glycan data
    glycan_data = {}
    for glycan_id, glycan_info in glycans.items():
        atoms_list = []
        for atom in glycan_info['atoms']:
            atoms_list.append({
                'record': atom['record'],
                'atom_number': atom['atom_number'],
                'atom_name': atom['atom_name'],
                'alt_loc': atom['alt_loc'],
                'residue_name': atom['residue_name'],
                'chain_id': atom['chain_id'],
                'residue_number': atom['residue_number'],
                'icode': atom['icode'],
                'x': atom['x'],
                'y': atom['y'],
                'z': atom['z'],
                'occupancy': atom['occupancy'],
                'temp_factor': atom['temp_factor'],
                'element': atom['element'],
                'charge': atom['charge'],
                'model': atom['model']
            })
        
        glycan_data[glycan_id] = {
            'id': glycan_id,
            'chain': glycan_info['chain'],
            'atoms': atoms_list,
            'residue_numbers': glycan_info['residue_numbers'],
            'start_residue': glycan_info['start_residue'],
            'end_residue': glycan_info['end_residue'],
            'simple_sequence': glycan_info['simple_sequence'],
            'full_sequence': glycan_info['full_sequence'],
            'unit_number': glycan_info['unit_number'],
            'residue_info': glycan_info['residue_info'],
            'residues': glycan_info['residues']
        }
    
    # Prepare other lines
    other_lines = pdb_data['other_lines']
    
    # Create complete structure
    complete_data = {
        'metadata': {
            'total_protein_atoms': len(protein_data),
            'total_glycans': len(glycan_data),
            'total_linkages': len(linkages),
            'glycan_ids': list(glycan_data.keys()),
            'detected_glycan_residues': pdb_data.get('detected_glycan_residues', []),
            'all_hetatm_residues': pdb_data.get('all_hetatm_residues', [])
        },
        'protein': protein_data,
        'glycans': glycan_data,
        'linkages': linkages,
        'other_lines': other_lines
    }
    
    return complete_data

def save_to_json(data, output_json):
    """
    Save complete data structure to JSON file.
    
    Parameters:
    -----------
    data : dict
        Complete data structure
    output_json : str
        Path to output JSON file
    """
    print(f"Saving complete data to {output_json}...")
    
    # Custom JSON encoder to handle numpy types
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                              np.int16, np.int32, np.int64, np.uint8,
                              np.uint16, np.uint32, np.uint64)):
                return int(obj)
            elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return super(NumpyEncoder, self).default(obj)
    
    with open(output_json, 'w') as f:
        json.dump(data, f, indent=2, cls=NumpyEncoder)
    
    print(f"  Saved {len(data['protein'])} protein atoms")
    print(f"  Saved {len(data['glycans'])} glycans")
    print(f"  Saved {len(data['linkages'])} linkages")
    print(f"  Detected glycan residues: {data['metadata']['detected_glycan_residues']}")

def reconstruct_pdb_from_json(json_data, output_pdb):
    """
    Reconstruct PDB file from JSON data (for verification).
    
    Parameters:
    -----------
    json_data : dict
        JSON data structure
    output_pdb : str
        Path to output PDB file
    """
    with open(output_pdb, 'w') as f:
        # Write protein atoms
        for atom in json_data['protein']:
            line = f"{atom['record']:6s}{atom['atom_number']:5d} {atom['atom_name']:4s}" \
                   f"{atom['alt_loc']:1s}{atom['residue_name']:3s} {atom['chain_id']:1s}" \
                   f"{atom['residue_number']:4d}{atom['icode']:1s}   " \
                   f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}" \
                   f"{atom['occupancy']:6.2f}{atom['temp_factor']:6.2f}" \
                   f"          {atom['element']:2s}{atom['charge']:2s}\n"
            f.write(line)
        
        # Write TER records from other_lines
        for line in json_data['other_lines']:
            if line.startswith('TER'):
                f.write(line + '\n')
        
        # Write glycan atoms
        for glycan_id, glycan_data in json_data['glycans'].items():
            for atom in glycan_data['atoms']:
                line = f"{atom['record']:6s}{atom['atom_number']:5d} {atom['atom_name']:4s}" \
                       f"{atom['alt_loc']:1s}{atom['residue_name']:3s} {atom['chain_id']:1s}" \
                       f"{atom['residue_number']:4d}{atom['icode']:1s}   " \
                       f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}" \
                       f"{atom['occupancy']:6.2f}{atom['temp_factor']:6.2f}" \
                       f"          {atom['element']:2s}{atom['charge']:2s}\n"
                f.write(line)
        
        # Write remaining other_lines
        for line in json_data['other_lines']:
            if not line.startswith('TER'):
                f.write(line + '\n')
        
        # Write END
        f.write("END\n")

def main():
    """Main function to execute all steps."""
    parser = argparse.ArgumentParser(description="Process glycosylated PDB files and extract complete information to JSON")
    parser.add_argument("--input_pdb", required=True, help="Path to input glycosylated PDB file")
    parser.add_argument("--output_json", required=True, help="Path to output JSON file")
    parser.add_argument("--reconstruct_pdb", help="Optional: Path to reconstruct PDB file for verification")
    
    args = parser.parse_args()
    
    input_pdb = args.input_pdb
    output_json = args.output_json
    reconstruct_pdb = args.reconstruct_pdb
    
    print("=" * 70)
    print("GLYCAN PROCESSING SCRIPT - JSON OUTPUT")
    print("=" * 70)
    
    print("\n1. PARSING COMPLETE PDB FILE")
    print("-" * 40)
    pdb_data = parse_pdb_complete(input_pdb)
    print(f"  Found {len(pdb_data['protein_atoms'])} protein atoms")
    print(f"  Found {len(pdb_data['glycan_atoms'])} glycan atoms")
    print(f"  Found {len(pdb_data['other_lines'])} other lines")
    
    print("\n2. EXTRACTING GLYCANS")
    print("-" * 40)
    glycans = extract_glycans_from_atoms(pdb_data['glycan_atoms'])
    print(f"  Found {len(glycans)} individual glycans")
    
    for i, (glycan_id, glycan_data) in enumerate(glycans.items(), 1):
        print(f"  Glycan {i}: {glycan_id}")
        print(f"    Chain: {glycan_data['chain']}")
        print(f"    Sequence: {glycan_data['simple_sequence']}")
        print(f"    Residue count: {glycan_data['unit_number']}")
        print(f"    Atoms: {len(glycan_data['atoms'])}")
        print(f"    Residue range: {glycan_data['start_residue']} to {glycan_data['end_residue']}")
        print(f"    Residues: {glycan_data['residues']}")
    
    print("\n3. FINDING GLYCAN-PROTEIN LINKAGES")
    print("-" * 40)
    linkages = find_glycan_linkages(pdb_data['protein_atoms'], glycans)
    
    print("\n4. CREATING COMPLETE DATA STRUCTURE")
    print("-" * 40)
    complete_data = create_complete_data_structure(pdb_data, glycans, linkages)
    
    print("\n5. SAVING TO JSON FILE")
    print("-" * 40)
    save_to_json(complete_data, output_json)
    
    if reconstruct_pdb:
        print("\n6. RECONSTRUCTING PDB FILE FOR VERIFICATION")
        print("-" * 40)
        reconstruct_pdb_from_json(complete_data, reconstruct_pdb)
        print(f"  Reconstructed PDB saved to {reconstruct_pdb}")
    
    print("\n" + "=" * 70)
    print("PROCESSING COMPLETE!")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  - Input PDB: {input_pdb}")
    print(f"  - Output JSON: {output_json}")
    print(f"  - Total protein atoms: {len(complete_data['protein'])}")
    print(f"  - Total glycans: {len(complete_data['glycans'])}")
    print(f"  - Total linkages: {len(complete_data['linkages'])}")
    print(f"  - Detected glycan residues: {complete_data['metadata']['detected_glycan_residues']}")
    
    if reconstruct_pdb:
        print(f"  - Reconstructed PDB: {reconstruct_pdb}")
    
    return complete_data

if __name__ == "__main__":
    data = main()
