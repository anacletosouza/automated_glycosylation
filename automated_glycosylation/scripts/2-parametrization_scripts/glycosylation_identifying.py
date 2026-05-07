#!/usr/bin/env python3
"""
Script to fix glycosidic linkages in PDB files from glycosylator.
Removes extra hydrogens and identifies proper glycosidic bonds.
"""

import os
import sys
import math
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional

# ============================================================================
# GEOMETRY FUNCTIONS
# ============================================================================

def calculate_atom_distance(coord1: Tuple[float, float, float], 
                          coord2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two 3D coordinates."""
    return math.sqrt(
        (coord1[0] - coord2[0])**2 +
        (coord1[1] - coord2[1])**2 +
        (coord1[2] - coord2[2])**2
    )

def extract_coordinates_from_line(line: str) -> Optional[Tuple[float, float, float]]:
    """Extract coordinates from a PDB line."""
    try:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        return (x, y, z)
    except (ValueError, IndexError):
        return None

# ============================================================================
# PDB PROCESSING FUNCTIONS
# ============================================================================

def remove_hydrogens(input_pdb: str, output_pdb: str):
    """
    Remove hydrogen atoms from a PDB file.
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

def parse_pdb(pdb_file: str):
    """
    Parse PDB file and separate ATOM (protein) and HETATM (glycan) records.
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

def extract_glycans(glycan_lines: List[str]) -> Dict[str, Dict]:
    """
    Extract individual glycans from HETATM records.
    """
    chains = defaultdict(list)
    for line in glycan_lines:
        chain_id = line[21:22].strip()
        chains[chain_id].append(line)
    
    glycans = {}
    glycan_counter = 1
    
    # Common carbohydrate residues
    glycan_residues = {'NDG', 'FCA', 'NAG', 'BMA', 'MAN', 'GAL', 'SIA', 'A2G'}
    
    for chain_id, chain_lines in chains.items():
        print(f"  Processing chain {chain_id} with {len(chain_lines)} glycan lines")
        
        residue_groups = defaultdict(list)
        for line in chain_lines:
            res_num = int(line[22:26].strip())
            residue_groups[res_num].append(line)
        
        sorted_residues = sorted(residue_groups.keys())
        
        # Find start points (NDG or A2G residues)
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
            
            # Extract residue names in order
            residue_info = {}
            for line in glycan_lines_list:
                res_name = line[17:20].strip()
                res_num = int(line[22:26].strip())
                if res_num not in residue_info:
                    residue_info[res_num] = res_name
            
            sorted_residue_numbers = sorted(residue_info.keys())
            
            # Create sequence strings
            sequence_parts = []
            simple_sequence_parts = []
            residue_counter = 1
            for res_num in sorted_residue_numbers:
                res_name = residue_info[res_num]
                sequence_parts.append(f"{res_name}{residue_counter}")
                simple_sequence_parts.append(res_name)
                residue_counter += 1
            
            sequence_poly = "_".join(sequence_parts)
            simple_sequence_poly = "_".join(simple_sequence_parts)
            
            glycans[glycan_id] = {
                'lines': glycan_lines_list,
                'chain': chain_id,
                'residue_numbers': glycan_residue_numbers,
                'start_residue': min(glycan_residue_numbers),
                'end_residue': max(glycan_residue_numbers),
                'residue_info': residue_info,
                'residue_sequence': sequence_poly,
                'simple_residue_sequence': simple_sequence_poly,
                'residues': simple_sequence_parts,
                'unit_number': len(glycan_residue_numbers),
            }
            
            print(f"    Glycan {glycan_id}: residues {glycan_residue_numbers[0]}-{glycan_residue_numbers[-1]}, "
                  f"atoms: {len(glycan_lines_list)}, sequence: {simple_sequence_poly}")
            
            # Remove processed residues
            sorted_residues = [r for r in sorted_residues if r > end_res]
            if not sorted_residues:
                break
    
    print(f"  Total glycans extracted: {len(glycans)}")
    return glycans

def find_protein_glycan_linkages(protein_lines: List[str], glycans: Dict[str, Dict]) -> List[Dict]:
    """
    Find linkages between glycans and protein residues.
    Filters by reasonable distance (1.0-2.0 Å for covalent bonds).
    """
    print("\nFinding protein-glycan linkages...")
    
    linkages = []
    protein_atoms = {}
    
    # Parse protein atoms
    for line in protein_lines:
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            res_num = int(line[22:26].strip())
            
            coords = extract_coordinates_from_line(line)
            if coords:
                protein_atoms[(chain_id, res_num, res_name, atom_name)] = {
                    'coordinates': coords,
                    'res_name': res_name,
                    'chain': chain_id,
                    'res_num': res_num,
                    'atom_name': atom_name,
                    'line': line
                }
    
    # Find linkages for each glycan
    for glycan_id, glycan_data in glycans.items():
        glycan_chain = glycan_data['chain']
        glycan_lines = glycan_data['lines']
        
        # Find C1 atom in first residue
        c1_coords = None
        c1_line = None
        
        # First try: look in start residue
        for line in glycan_lines:
            res_num = int(line[22:26].strip())
            if res_num == glycan_data['start_residue']:
                atom_name = line[12:16].strip()
                if atom_name == 'C1' or 'C1' in atom_name:
                    c1_coords = extract_coordinates_from_line(line)
                    c1_line = line
                    break
        
        # Second try: look for any C1 in first residue type
        if not c1_coords:
            first_residue_type = glycan_data['residues'][0]
            for line in glycan_lines:
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                if res_name == first_residue_type and (atom_name == 'C1' or 'C1' in atom_name):
                    c1_coords = extract_coordinates_from_line(line)
                    c1_line = line
                    break
        
        if not c1_coords:
            print(f"  Warning: No C1 atom found for glycan {glycan_id}")
            continue
        
        # Determine linking type
        linking_type = "N-linked" if glycan_data['residues'][0] == 'NDG' else "O-linked"
        
        # Find target protein atoms
        target_atoms = []
        for (chain, res_num, res_name, atom_name), atom_data in protein_atoms.items():
            if chain != glycan_chain:
                continue
            
            if linking_type == "N-linked":
                if res_name == 'ASN' and atom_name == 'ND2':
                    target_atoms.append(atom_data)
            else:  # O-linked
                if (res_name == 'SER' and atom_name == 'OG') or \
                   (res_name == 'THR' and atom_name == 'OG1'):
                    target_atoms.append(atom_data)
        
        # Find closest protein atom
        min_distance = float('inf')
        closest_atom = None
        
        for atom_data in target_atoms:
            distance = calculate_atom_distance(c1_coords, atom_data['coordinates'])
            if distance < min_distance:
                min_distance = distance
                closest_atom = atom_data
        
        # Filter by reasonable distance for covalent bond
        if closest_atom and 1.0 <= min_distance <= 2.0:  # Typical covalent bond range
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
            print(f"  Found {linking_type} linkage: {linkage_info['site_protein_residue']} "
                  f"-> {glycan_id} (distance: {min_distance:.3f} Å)")
        else:
            if closest_atom:
                print(f"  Skipped {linking_type} linkage for {glycan_id}: "
                      f"distance {min_distance:.3f} Å is outside covalent bond range")
            else:
                print(f"  Warning: No suitable protein atom found for glycan {glycan_id}")
    
    print(f"  Total protein-glycan linkages found: {len(linkages)}")
    return linkages

def find_glycosidic_linkages_within_glycans(glycans: Dict[str, Dict]) -> List[Dict]:
    """
    Find glycosidic linkages between consecutive residues within glycans.
    IMPROVED VERSION: Better detection of glycosidic bonds between carboidratos.
    """
    print("\nFinding glycosidic linkages within glycans...")
    
    all_linkages = []
    
    for glycan_id, glycan_data in glycans.items():
        if glycan_data['unit_number'] < 2:
            continue
        
        print(f"  Analyzing glycan {glycan_id} ({glycan_data['unit_number']} residues)")
        
        # Create a dictionary to map residue numbers to atom data
        residue_atoms = defaultdict(dict)
        
        for line in glycan_data['lines']:
            res_num = int(line[22:26].strip())
            atom_name = line[12:16].strip()
            coords = extract_coordinates_from_line(line)
            element = line[76:78].strip()
            
            if coords:
                residue_atoms[res_num][atom_name] = {
                    'name': atom_name,
                    'coordinates': coords,
                    'element': element,
                    'line': line
                }
        
        # Get residue numbers in order
        residue_numbers = sorted(glycan_data['residue_info'].keys())
        
        # Check linkages between consecutive residues
        for i in range(len(residue_numbers) - 1):
            res1_num = residue_numbers[i]
            res2_num = residue_numbers[i + 1]
            
            res1_name = glycan_data['residue_info'][res1_num]
            res2_name = glycan_data['residue_info'][res2_num]
            
            # Look for C1 (anomeric carbon) in residue 1
            # Try different possible names for C1
            c1_candidates = []
            for atom_name in residue_atoms[res1_num]:
                if atom_name == 'C1' or 'C1' in atom_name:
                    c1_candidates.append(residue_atoms[res1_num][atom_name])
            
            if not c1_candidates:
                # If no C1 found, try to find any carbon that could be anomeric
                for atom_name, atom in residue_atoms[res1_num].items():
                    if atom['element'] == 'C' and atom_name.startswith('C'):
                        c1_candidates.append(atom)
            
            if not c1_candidates:
                print(f"    Warning: No suitable carbon atom found in residue {res1_name}{res1_num}")
                continue
            
            # Look for potential glycosidic oxygens in residue 2
            # Common glycosidic oxygen positions
            potential_oxygens = []
            for atom_name, atom in residue_atoms[res2_num].items():
                if atom['element'] == 'O':
                    # Check for common glycosidic oxygen names
                    if atom_name in ['O2', 'O3', 'O4', 'O6', 'O1']:
                        potential_oxygens.append(atom)
            
            if not potential_oxygens:
                # If no named oxygens, try any oxygen
                for atom_name, atom in residue_atoms[res2_num].items():
                    if atom['element'] == 'O':
                        potential_oxygens.append(atom)
            
            if not potential_oxygens:
                print(f"    Warning: No oxygen atom found in residue {res2_name}{res2_num}")
                continue
            
            # Now check all combinations of C1 candidates and oxygen candidates
            best_linkage = None
            best_distance = float('inf')
            
            for c1_atom in c1_candidates:
                c1_coords = c1_atom['coordinates']
                
                for o_atom in potential_oxygens:
                    o_coords = o_atom['coordinates']
                    dist = calculate_atom_distance(c1_coords, o_coords)
                    
                    # Check if distance is reasonable for glycosidic bond
                    # Typical C-O glycosidic bond: 1.42 ± 0.15 Å
                    if 1.2 <= dist <= 1.6 and dist < best_distance:
                        # Also check if there's a carbon attached to this oxygen in residue 2
                        # (to confirm it's a glycosidic oxygen, not a hydroxyl oxygen)
                        o_has_carbon_neighbor = False
                        for atom_name2, atom2 in residue_atoms[res2_num].items():
                            if atom2['element'] == 'C':
                                dist_co = calculate_atom_distance(o_coords, atom2['coordinates'])
                                if 1.2 <= dist_co <= 1.6:  # C-O bond distance
                                    o_has_carbon_neighbor = True
                                    break
                        
                        if o_has_carbon_neighbor:
                            best_distance = dist
                            best_linkage = {
                                'c1_atom': c1_atom,
                                'o_atom': o_atom,
                                'distance': dist
                            }
            
            if best_linkage:
                c1_atom = best_linkage['c1_atom']
                o_atom = best_linkage['o_atom']
                distance = best_linkage['distance']
                
                # Determine linkage type based on oxygen name
                linkage_type = None
                if o_atom['name'] == 'O2':
                    linkage_type = '1-2'
                elif o_atom['name'] == 'O3':
                    linkage_type = '1-3'
                elif o_atom['name'] == 'O4':
                    linkage_type = '1-4'
                elif o_atom['name'] == 'O6':
                    linkage_type = '1-6'
                elif o_atom['name'] == 'O1':
                    linkage_type = '2-1'  # For sialic acid
                else:
                    # Unknown oxygen position
                    linkage_type = f"1-? (O:{o_atom['name']})"
                
                linkage_info = {
                    'glycan': glycan_id,
                    'residue1': f"{res1_name}{res1_num}",
                    'residue2': f"{res2_name}{res2_num}",
                    'atom1': c1_atom['name'],
                    'atom2': o_atom['name'],
                    'linkage_type': linkage_type,
                    'distance': distance,
                    'is_valid': True
                }
                
                all_linkages.append(linkage_info)
                print(f"    Found {linkage_type} linkage: {res1_name}{res1_num}:{c1_atom['name']} "
                      f"- {res2_name}{res2_num}:{o_atom['name']} (dist: {distance:.3f} Å)")
        
        # Also check for sialic acid linkages (which can be 2-6 or 2-3)
        # SIA residues often link to GAL residues
        for i in range(len(residue_numbers)):
            res_num = residue_numbers[i]
            res_name = glycan_data['residue_info'][res_num]
            
            if res_name == 'SIA':
                # Look for C2 in SIA residue
                c2_atom = None
                for atom_name, atom in residue_atoms[res_num].items():
                    if atom_name == 'C2' or 'C2' in atom_name:
                        c2_atom = atom
                        break
                
                if c2_atom:
                    # Look for potential linkages to next residue
                    if i < len(residue_numbers) - 1:
                        next_res_num = residue_numbers[i + 1]
                        next_res_name = glycan_data['residue_info'][next_res_num]
                        
                        # Check for oxygens in next residue
                        for atom_name, atom in residue_atoms[next_res_num].items():
                            if atom['element'] == 'O' and atom_name in ['O6', 'O3', 'O4']:
                                dist = calculate_atom_distance(c2_atom['coordinates'], atom['coordinates'])
                                if 1.2 <= dist <= 1.6:
                                    linkage_type = f"2-{atom_name[1:]}"  # 2-6, 2-3, etc.
                                    
                                    linkage_info = {
                                        'glycan': glycan_id,
                                        'residue1': f"{res_name}{res_num}",
                                        'residue2': f"{next_res_name}{next_res_num}",
                                        'atom1': c2_atom['name'],
                                        'atom2': atom['name'],
                                        'linkage_type': linkage_type,
                                        'distance': dist,
                                        'is_valid': True
                                    }
                                    
                                    all_linkages.append(linkage_info)
                                    print(f"    Found {linkage_type} linkage (sialic acid): "
                                          f"{res_name}{res_num}:{c2_atom['name']} - "
                                          f"{next_res_name}{next_res_num}:{atom['name']} (dist: {dist:.3f} Å)")
    
    print(f"  Total glycosidic linkages found: {len(all_linkages)}")
    return all_linkages

def identify_hydrogens_on_glycosidic_oxygens(orig_glycans: Dict[str, Dict], 
                                           glycosidic_linkages: List[Dict]) -> Tuple[List[str], List[str]]:
    """
    Identify hydrogen atoms that need to be removed from glycosidic oxygens.
    Also identify overlapping intra-glycan oxygens that need to be removed.
    Returns tuple: (hydrogens_to_remove, oxygens_to_remove)
    """
    print("\nIdentifying hydrogens on glycosidic oxygens and overlapping oxygens...")
    
    hydrogens_to_remove = []
    oxygens_to_remove = []
    
    # Build lookup for glycan atoms by residue and atom name
    for glycan_id, glycan_data in orig_glycans.items():
        # Create atom lookup with serial numbers and coordinates
        atom_lookup = {}
        coord_lookup = {}
        
        for line in glycan_data['lines']:
            try:
                serial = line[6:11].strip()
                atom_name = line[12:16].strip()
                res_num = int(line[22:26].strip())
                element = line[76:78].strip()
                coords = extract_coordinates_from_line(line)
                
                if not coords:
                    continue
                    
                key = (res_num, atom_name)
                atom_lookup[key] = {
                    'serial': serial,
                    'element': element,
                    'line': line,
                    'coordinates': coords
                }
                coord_lookup[serial] = coords
            except (ValueError, IndexError):
                continue
        
        # Check linkages in this glycan
        for linkage in glycosidic_linkages:
            if linkage['glycan'] != glycan_id:
                continue
            
            # Extract residue number from residue string
            import re
            res2_str = linkage['residue2']
            
            # Try to extract the number part
            match = re.search(r'(\d+)$', res2_str)
            if match:
                res2_num = int(match.group(1))
            else:
                # If no number found, try to get from the original residue numbers
                res2_name = linkage['residue2'][:3]  # First 3 chars are residue name
                for res_num, res_name in glycan_data['residue_info'].items():
                    if res_name == res2_name:
                        res2_num = res_num
                        break
                else:
                    continue
            
            oxygen_name = linkage['atom2']
            oxygen_key = (res2_num, oxygen_name)
            
            if oxygen_key not in atom_lookup:
                continue
            
            # Get oxygen serial and coordinates
            oxygen_serial = atom_lookup[oxygen_key]['serial']
            oxygen_coords = atom_lookup[oxygen_key]['coordinates']
            
            print(f"  Processing glycosidic oxygen: {linkage['residue2']}:{oxygen_name} (serial: {oxygen_serial})")
            
            # Find hydrogens close to this oxygen (within typical O-H bond distance ~1.0 Å)
            # Use a larger radius (1.3 Å) to be safe
            search_radius = 1.3  # Å
            
            for (res_num_h, atom_name_h), atom_data in atom_lookup.items():
                if atom_data['element'] != 'H':
                    continue
                
                if res_num_h != res2_num:
                    continue
                
                h_coords = atom_data['coordinates']
                distance = calculate_atom_distance(oxygen_coords, h_coords)
                
                if distance <= search_radius:
                    h_serial = atom_data['serial']
                    if h_serial not in hydrogens_to_remove:
                        hydrogens_to_remove.append(h_serial)
                        print(f"    Found hydrogen within {search_radius} Å: H{h_serial}:{atom_name_h} "
                              f"(dist: {distance:.3f} Å)")
    
    # Now identify overlapping intra-glycan oxygens
    print("\nIdentifying overlapping intra-glycan oxygens...")
    
    # Collect all oxygen coordinates from all glycans
    all_oxygen_data = []
    
    for glycan_id, glycan_data in orig_glycans.items():
        for line in glycan_data['lines']:
            try:
                serial = line[6:11].strip()
                atom_name = line[12:16].strip()
                res_num = int(line[22:26].strip())
                element = line[76:78].strip()
                coords = extract_coordinates_from_line(line)
                
                if element == 'O' and coords:
                    all_oxygen_data.append({
                        'serial': serial,
                        'atom_name': atom_name,
                        'res_num': res_num,
                        'glycan_id': glycan_id,
                        'coordinates': coords,
                        'line': line
                    })
            except (ValueError, IndexError):
                continue
    
    print(f"  Total oxygens found: {len(all_oxygen_data)}")
    
    # Find overlapping oxygens (distance < 0.5 Å)
    overlap_distance = 0.5  # Å
    
    # Use a grid to speed up proximity search
    grid_size = 1.0  # Å
    grid = defaultdict(list)
    
    for oxygen in all_oxygen_data:
        x, y, z = oxygen['coordinates']
        grid_x = int(x / grid_size)
        grid_y = int(y / grid_size)
        grid_z = int(z / grid_size)
        grid[(grid_x, grid_y, grid_z)].append(oxygen)
    
    # Check each grid cell and neighboring cells
    processed_pairs = set()
    
    for (gx, gy, gz), oxygens_in_cell in grid.items():
        # Check pairs within this cell
        for i in range(len(oxygens_in_cell)):
            o1 = oxygens_in_cell[i]
            
            for j in range(i + 1, len(oxygens_in_cell)):
                o2 = oxygens_in_cell[j]
                
                pair_key = tuple(sorted((o1['serial'], o2['serial'])))
                if pair_key in processed_pairs:
                    continue
                
                distance = calculate_atom_distance(o1['coordinates'], o2['coordinates'])
                
                if distance < overlap_distance:
                    processed_pairs.add(pair_key)
                    
                    print(f"    Found overlapping oxygens: O{o1['serial']} ({o1['glycan_id']}:{o1['res_num']}:{o1['atom_name']}) "
                          f"and O{o2['serial']} ({o2['glycan_id']}:{o2['res_num']}:{o2['atom_name']}) "
                          f"(dist: {distance:.3f} Å)")
                    
                    # Decide which one to keep and which to remove
                    # Keep the one that is part of a glycosidic linkage
                    o1_is_glycosidic = False
                    o2_is_glycosidic = False
                    
                    for linkage in glycosidic_linkages:
                        if linkage['glycan'] == o1['glycan_id']:
                            # Extract residue number from linkage
                            match = re.search(r'(\d+)$', linkage['residue2'])
                            if match and int(match.group(1)) == o1['res_num'] and linkage['atom2'] == o1['atom_name']:
                                o1_is_glycosidic = True
                        
                        if linkage['glycan'] == o2['glycan_id']:
                            match = re.search(r'(\d+)$', linkage['residue2'])
                            if match and int(match.group(1)) == o2['res_num'] and linkage['atom2'] == o2['atom_name']:
                                o2_is_glycosidic = True
                    
                    # Remove the oxygen that is NOT part of a glycosidic linkage
                    # If both or neither are glycosidic, remove the one with higher serial number
                    if o1_is_glycosidic and not o2_is_glycosidic:
                        # Keep o1, remove o2
                        if o2['serial'] not in oxygens_to_remove:
                            oxygens_to_remove.append(o2['serial'])
                            print(f"      Removing O{o2['serial']} (not glycosidic)")
                    elif o2_is_glycosidic and not o1_is_glycosidic:
                        # Keep o2, remove o1
                        if o1['serial'] not in oxygens_to_remove:
                            oxygens_to_remove.append(o1['serial'])
                            print(f"      Removing O{o1['serial']} (not glycosidic)")
                    else:
                        # Both or neither are glycosidic, remove the one with higher serial number
                        if o1['serial'] > o2['serial']:
                            if o1['serial'] not in oxygens_to_remove:
                                oxygens_to_remove.append(o1['serial'])
                                print(f"      Removing O{o1['serial']} (higher serial)")
                        else:
                            if o2['serial'] not in oxygens_to_remove:
                                oxygens_to_remove.append(o2['serial'])
                                print(f"      Removing O{o2['serial']} (higher serial)")
        
        # Check neighboring cells for close oxygens
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    
                    neighbor_key = (gx + dx, gy + dy, gz + dz)
                    if neighbor_key in grid:
                        for o1 in oxygens_in_cell:
                            for o2 in grid[neighbor_key]:
                                pair_key = tuple(sorted((o1['serial'], o2['serial'])))
                                if pair_key in processed_pairs:
                                    continue
                                
                                distance = calculate_atom_distance(o1['coordinates'], o2['coordinates'])
                                
                                if distance < overlap_distance:
                                    processed_pairs.add(pair_key)
                                    
                                    print(f"    Found overlapping oxygens (neighbor cells): O{o1['serial']} ({o1['glycan_id']}:{o1['res_num']}:{o1['atom_name']}) "
                                          f"and O{o2['serial']} ({o2['glycan_id']}:{o2['res_num']}:{o2['atom_name']}) "
                                          f"(dist: {distance:.3f} Å)")
                                    
                                    # Similar logic to decide which to remove
                                    o1_is_glycosidic = False
                                    o2_is_glycosidic = False
                                    
                                    for linkage in glycosidic_linkages:
                                        if linkage['glycan'] == o1['glycan_id']:
                                            match = re.search(r'(\d+)$', linkage['residue2'])
                                            if match and int(match.group(1)) == o1['res_num'] and linkage['atom2'] == o1['atom_name']:
                                                o1_is_glycosidic = True
                                        
                                        if linkage['glycan'] == o2['glycan_id']:
                                            match = re.search(r'(\d+)$', linkage['residue2'])
                                            if match and int(match.group(1)) == o2['res_num'] and linkage['atom2'] == o2['atom_name']:
                                                o2_is_glycosidic = True
                                    
                                    if o1_is_glycosidic and not o2_is_glycosidic:
                                        if o2['serial'] not in oxygens_to_remove:
                                            oxygens_to_remove.append(o2['serial'])
                                            print(f"      Removing O{o2['serial']} (not glycosidic)")
                                    elif o2_is_glycosidic and not o1_is_glycosidic:
                                        if o1['serial'] not in oxygens_to_remove:
                                            oxygens_to_remove.append(o1['serial'])
                                            print(f"      Removing O{o1['serial']} (not glycosidic)")
                                    else:
                                        if o1['serial'] > o2['serial']:
                                            if o1['serial'] not in oxygens_to_remove:
                                                oxygens_to_remove.append(o1['serial'])
                                                print(f"      Removing O{o1['serial']} (higher serial)")
                                        else:
                                            if o2['serial'] not in oxygens_to_remove:
                                                oxygens_to_remove.append(o2['serial'])
                                                print(f"      Removing O{o2['serial']} (higher serial)")
    
    print(f"\n  Total hydrogens on glycosidic oxygens identified for removal: {len(hydrogens_to_remove)}")
    print(f"  Total overlapping oxygens identified for removal: {len(oxygens_to_remove)}")
    
    # Also find hydrogens attached to oxygens that will be removed
    if oxygens_to_remove:
        print("\n  Finding hydrogens attached to oxygens that will be removed...")
        
        # Build a map of oxygen coordinates to serials
        oxygen_serial_to_coords = {}
        for oxygen in all_oxygen_data:
            if oxygen['serial'] in oxygens_to_remove:
                oxygen_serial_to_coords[oxygen['serial']] = oxygen['coordinates']
        
        # Find hydrogens close to these oxygens
        for glycan_id, glycan_data in orig_glycans.items():
            for line in glycan_data['lines']:
                try:
                    serial = line[6:11].strip()
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip()
                    coords = extract_coordinates_from_line(line)
                    
                    if element != 'H' or not coords:
                        continue
                    
                    # Check if this hydrogen is close to any oxygen being removed
                    for oxy_serial, oxy_coords in oxygen_serial_to_coords.items():
                        distance = calculate_atom_distance(coords, oxy_coords)
                        if distance <= 1.3:  # Typical O-H bond distance
                            if serial not in hydrogens_to_remove:
                                hydrogens_to_remove.append(serial)
                                print(f"    Found hydrogen H{serial} attached to oxygen O{oxy_serial} "
                                      f"(dist: {distance:.3f} Å)")
                except (ValueError, IndexError):
                    continue
        
        print(f"  Additional hydrogens found (attached to removed oxygens): "
              f"{len(hydrogens_to_remove) - len([h for h in hydrogens_to_remove if h in locals().get('prev_h_count', [])])}")
    
    return hydrogens_to_remove, oxygens_to_remove

def identify_anomeric_hydrogens(orig_glycans: Dict[str, Dict], 
                               protein_linkages: List[Dict]) -> List[str]:
    """
    Identify anomeric hydrogens (H1) that need to be removed from reducing ends.
    """
    print("\nIdentifying anomeric hydrogens...")
    
    hydrogens_to_remove = []
    
    for glycan_id, glycan_data in orig_glycans.items():
        # Check if this glycan has a protein linkage
        has_protein_linkage = any(l['glycan_binding'] == glycan_id for l in protein_linkages)
        
        if has_protein_linkage:
            # Get first residue number and name
            first_res_num = glycan_data['start_residue']
            first_res_name = glycan_data['residues'][0]
            
            # Look for H1 hydrogens in first residue
            for line in glycan_data['lines']:
                try:
                    res_num = int(line[22:26].strip())
                    if res_num == first_res_num:
                        atom_name = line[12:16].strip()
                        element = line[76:78].strip()
                        
                        # Check for various H1 naming patterns
                        is_anomeric_h = False
                        if element == 'H':
                            if atom_name == 'H1' or 'H1' in atom_name:
                                is_anomeric_h = True
                            elif atom_name.startswith('H') and '1' in atom_name:
                                # Check if it's H1 with variations (H1A, H1B, etc.)
                                if atom_name[1:].replace('A', '').replace('B', '').replace('C', '').replace('D', '').startswith('1'):
                                    is_anomeric_h = True
                        
                        if is_anomeric_h:
                            serial = line[6:11].strip()
                            if serial not in hydrogens_to_remove:
                                hydrogens_to_remove.append(serial)
                                print(f"  Marked anomeric hydrogen for removal: {first_res_name}{first_res_num}:{atom_name}")
                except (ValueError, IndexError):
                    continue
    
    print(f"  Total anomeric hydrogens identified for removal: {len(hydrogens_to_remove)}")
    return hydrogens_to_remove

def create_fixed_pdb(input_pdb: str, output_pdb: str, 
                    hydrogens_to_remove: List[str], oxygens_to_remove: List[str]):
    """
    Create a fixed PDB file with specified hydrogens and oxygens removed.
    """
    print(f"\nCreating fixed PDB: {output_pdb}")
    
    hydrogens_set = set(hydrogens_to_remove)
    oxygens_set = set(oxygens_to_remove)
    atoms_removed = 0
    h_removed = 0
    o_removed = 0
    
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    serial = line[6:11].strip()
                    element = line[76:78].strip()
                    
                    if serial in hydrogens_set:
                        atoms_removed += 1
                        h_removed += 1
                        continue
                    elif serial in oxygens_set:
                        atoms_removed += 1
                        o_removed += 1
                        continue
                except (ValueError, IndexError):
                    pass
            f_out.write(line)
    
    print(f"  Removed {atoms_removed} atoms total:")
    print(f"    - {h_removed} hydrogen atoms")
    print(f"    - {o_removed} oxygen atoms (overlapping)")

def save_linkage_report(glycosidic_linkages: List[Dict], 
                       protein_linkages: List[Dict], 
                       output_file: str):
    """Save a detailed report of all linkages."""
    print(f"\nSaving linkage report to: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("GLYCOSIDIC LINKAGE ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("SUMMARY:\n")
        f.write(f"  Total glycans analyzed: {len(set(l['glycan'] for l in glycosidic_linkages))}\n")
        f.write(f"  Total glycosidic linkages (carb-carb): {len(glycosidic_linkages)}\n")
        f.write(f"  Total protein-glycan linkages: {len(protein_linkages)}\n\n")
        
        f.write("GLYCOSIDIC LINKAGES WITHIN GLYCANS:\n")
        f.write("-" * 80 + "\n")
        if glycosidic_linkages:
            # Group by glycan
            by_glycan = defaultdict(list)
            for linkage in glycosidic_linkages:
                by_glycan[linkage['glycan']].append(linkage)
            
            for glycan_id, linkages in sorted(by_glycan.items()):
                f.write(f"\nGlycan: {glycan_id}\n")
                for linkage in linkages:
                    f.write(f"  {linkage['residue1']}:{linkage['atom1']} - "
                           f"{linkage['residue2']}:{linkage['atom2']}\n")
                    f.write(f"    Type: {linkage['linkage_type']}, "
                           f"Distance: {linkage['distance']:.3f} Å\n")
        else:
            f.write("No glycosidic linkages found within glycans.\n\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("PROTEIN-GLYCAN LINKAGES:\n")
        f.write("-" * 80 + "\n")
        if protein_linkages:
            # Group by protein residue
            by_protein = defaultdict(list)
            for linkage in protein_linkages:
                key = linkage['site_protein_residue']
                by_protein[key].append(linkage)
            
            for protein_res, linkages in sorted(by_protein.items()):
                f.write(f"\nProtein site: {protein_res}\n")
                for linkage in linkages:
                    f.write(f"  Glycan: {linkage['glycan_binding']}\n")
                    f.write(f"    Sequence: {linkage['simple_sequence_poly']}\n")
                    f.write(f"    Units: {linkage['unit_number']}, Type: {linkage['linking_type']}\n")
                    f.write(f"    Protein atom: {linkage['protein_atom']}\n")
                    f.write(f"    Glycan atom: {linkage['glycan_atom']}\n")
                    f.write(f"    Distance: {linkage['distance']} Å\n")
        else:
            f.write("No protein-glycan linkages found.\n")
    
    print(f"  Report saved successfully")

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_glycosidic_pdb.py <input_pdb> [output_pdb]")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    
    if len(sys.argv) >= 3:
        output_pdb = sys.argv[2]
    else:
        base, ext = os.path.splitext(input_pdb)
        output_pdb = f"{base}_fixed{ext}"
    
    if not os.path.exists(input_pdb):
        print(f"ERROR: Input PDB file not found: {input_pdb}")
        sys.exit(1)
    
    print("=" * 80)
    print("GLYCOSIDIC LINKAGE FIXER")
    print("=" * 80)
    
    # Step 1: First remove all hydrogens to simplify analysis
    temp_pdb = output_pdb.replace(".pdb", "_noh.pdb")
    remove_hydrogens(input_pdb, temp_pdb)
    
    # Step 2: Parse PDB
    protein_lines, glycan_lines = parse_pdb(temp_pdb)
    print(f"\nParsed PDB:")
    print(f"  Protein atoms: {len(protein_lines)}")
    print(f"  Glycan atoms: {len(glycan_lines)}")
    
    # Step 3: Extract glycans
    glycans = extract_glycans(glycan_lines)
    
    if not glycans:
        print("\nNo glycans found. Exiting.")
        return
    
    # Step 4: Find protein-glycan linkages
    protein_linkages = find_protein_glycan_linkages(protein_lines, glycans)
    
    # Step 5: Find glycosidic linkages within glycans
    glycosidic_linkages = find_glycosidic_linkages_within_glycans(glycans)
    
    # Step 6: Identify hydrogens to remove (need original PDB with hydrogens)
    print("\n" + "=" * 80)
    print("ANALYZING ORIGINAL PDB FOR HYDROGENS AND OVERLAPPING OXYGENS")
    print("=" * 80)
    
    # Parse original PDB to get glycans with hydrogens
    orig_protein_lines, orig_glycan_lines = parse_pdb(input_pdb)
    orig_glycans = extract_glycans(orig_glycan_lines)
    
    # Identify hydrogens on glycosidic oxygens AND overlapping oxygens
    hydrogens_on_oxygens, overlapping_oxygens = identify_hydrogens_on_glycosidic_oxygens(
        orig_glycans, glycosidic_linkages
    )
    
    # Identify anomeric hydrogens
    anomeric_hydrogens = identify_anomeric_hydrogens(orig_glycans, protein_linkages)
    
    # Combine all atoms to remove
    all_hydrogens_to_remove = hydrogens_on_oxygens + anomeric_hydrogens
    
    # Step 7: Create final fixed PDB
    create_fixed_pdb(input_pdb, output_pdb, all_hydrogens_to_remove, overlapping_oxygens)
    
    # Step 8: Save linkage report
    report_file = output_pdb.replace(".pdb", "_linkages.txt")
    save_linkage_report(glycosidic_linkages, protein_linkages, report_file)
    
    # Step 9: Cleanup
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)
    
    # Final summary
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Input PDB: {input_pdb}")
    print(f"Fixed PDB: {output_pdb}")
    print(f"Linkage report: {report_file}")
    print(f"\nStatistics:")
    print(f"  Glycans identified: {len(glycans)}")
    print(f"  Glycosidic linkages (carb-carb): {len(glycosidic_linkages)}")
    print(f"  Protein-glycan linkages: {len(protein_linkages)}")
    print(f"  Hydrogens on glycosidic oxygens removed: {len(hydrogens_on_oxygens)}")
    print(f"  Anomeric hydrogens removed: {len(anomeric_hydrogens)}")
    print(f"  Overlapping oxygens removed: {len(overlapping_oxygens)}")
    print(f"  Total atoms removed: {len(all_hydrogens_to_remove) + len(overlapping_oxygens)}")
    
    # Show some examples of glycosidic linkages found
    if glycosidic_linkages:
        print(f"\nExamples of glycosidic linkages found:")
        for linkage in glycosidic_linkages[:10]:  # Show first 10
            print(f"  {linkage['glycan']}: {linkage['residue1']}:{linkage['atom1']} - "
                  f"{linkage['residue2']}:{linkage['atom2']} ({linkage['linkage_type']}, {linkage['distance']:.3f} Å)")
    
    # Show top protein-glycan linkages
    if protein_linkages:
        print(f"\nTop protein-glycan linkages (sorted by distance):")
        sorted_linkages = sorted(protein_linkages, key=lambda x: x['distance'])[:10]
        for linkage in sorted_linkages:
            print(f"  {linkage['site_protein_residue']} -> {linkage['glycan_binding']}: "
                  f"{linkage['distance']} Å ({linkage['linking_type']})")
    
    print("=" * 80)

if __name__ == "__main__":
    main()
