#!/usr/bin/env python3
"""
Script to process PDB glycoproteins and create carbohydrate residue variants
with charge neutralization and proper HDB generation compatible with GROMACS.
"""

import os
import re
import json
import sys
import argparse
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional

# List of standard amino acids to exclude
STANDARD_AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
               'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
               'TYR', 'VAL'}

# ============================================================================
# PDB PARSER
# ============================================================================

class CarbohydratePDBParser:
    def __init__(self, pdb_content):
        self.pdb_content = pdb_content
        self.residues = []
        self.atoms = []
        self.parse_pdb()
    
    def parse_pdb(self):
        lines = self.pdb_content.strip().split('\n')
        
        for line in lines:
            if line.startswith('HETATM'):
                atom_data = self.parse_atom_line(line)
                if atom_data:
                    self.atoms.append(atom_data)
        
        self.group_atoms_by_residue()
    
    def parse_atom_line(self, line):
        try:
            serial = int(line[6:11].strip())
            name = line[12:16].strip()
            alt_loc = line[16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            return {
                'serial': serial,
                'name': name,
                'alt_loc': alt_loc,
                'res_name': res_name,
                'chain_id': chain_id,
                'res_seq': res_seq,
                'x': x, 'y': y, 'z': z,
                'element': line[76:78].strip() if len(line) > 76 else ''
            }
        except (ValueError, IndexError) as e:
            print(f"Error parsing line: {line}")
            return None
    
    def group_atoms_by_residue(self):
        residues_dict = {}
        
        for atom in self.atoms:
            key = (atom['res_name'], atom['chain_id'], atom['res_seq'])
            if key not in residues_dict:
                residues_dict[key] = {
                    'name_pdb': atom['res_name'],
                    'chain_carb': atom['chain_id'],
                    'res_seq': atom['res_seq'],
                    'atoms': []
                }
            residues_dict[key]['atoms'].append(atom)
        
        self.residues = sorted(residues_dict.values(), key=lambda x: x['res_seq'])

# ============================================================================
# RTP PARSER
# ============================================================================

def parse_rtp_file(file_path):
    residues = {}
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Split into residue sections
    pattern = r'\[ (\w{3}) \]'
    sections = re.split(pattern, content)
    
    for i in range(1, len(sections), 2):
        res_name = sections[i].strip()
        section_content = sections[i + 1]
        
        # Parse comments
        lines = section_content.strip().split('\n')
        comments = []
        atoms_start = bonds_start = impropers_start = -1
        
        for idx, line in enumerate(lines):
            if line.strip().startswith('[ atoms ]'):
                atoms_start = idx
            elif line.strip().startswith('[ bonds ]'):
                bonds_start = idx
            elif line.strip().startswith('[ impropers ]'):
                impropers_start = idx
            elif idx < atoms_start or atoms_start == -1:
                comments.append(line.strip())
        
        # Parse atoms
        atoms = []
        if atoms_start != -1:
            if bonds_start != -1:
                atom_lines = lines[atoms_start + 1:bonds_start]
            elif impropers_start != -1:
                atom_lines = lines[atoms_start + 1:impropers_start]
            else:
                atom_lines = lines[atoms_start + 1:]
            
            for line in atom_lines:
                line = line.strip()
                if line and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) >= 4:
                        atoms.append({
                            'name': parts[0],
                            'type': parts[1],
                            'charge': float(parts[2]),
                            'charge_group': int(parts[3])
                        })
        
        # Parse bonds
        bonds = []
        if bonds_start != -1:
            if impropers_start != -1:
                bond_lines = lines[bonds_start + 1:impropers_start]
            else:
                bond_lines = lines[bonds_start + 1:]
            
            for line in bond_lines:
                line = line.strip()
                if line and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) >= 2:
                        bonds.append((parts[0], parts[1]))
        
        # Parse impropers
        impropers = []
        if impropers_start != -1:
            improper_lines = lines[impropers_start + 1:]
            for line in improper_lines:
                line = line.strip()
                if line and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) >= 4:
                        impropers.append(parts[:4])
        
        # Extract additional info from comments
        name_charmm = ''
        name_iupac = ''
        for line in comments:
            if 'Original CHARMM name:' in line:
                name_charmm = line.split('Original CHARMM name:')[1].split(';')[0].strip()
            if 'IUPAC name:' in line:
                name_iupac = line.split('IUPAC name:')[1].split(';')[0].strip()
        
        residues[res_name] = {
            'name_charmm': name_charmm,
            'name_iupac': name_iupac,
            'comments': comments,
            'atoms': atoms,
            'bonds': bonds,
            'impropers': impropers
        }
    
    return residues

# ============================================================================
# HDB PARSER - CORRECTED FOR GROMACS FORMAT
# ============================================================================

def parse_hdb_file(file_path):
    """Parse HDB file into a dictionary by residue name (PDB codes)"""
    residues = {}
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    i = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Skip comments and empty lines
        if not line or line.startswith(';') or line.startswith('#'):
            i += 1
            continue
        
        parts = line.split()
        
        # Check if this is a residue header line
        if len(parts) >= 2:
            try:
                res_name = parts[0]
                num_donors = int(parts[1])
                
                donors = []
                i += 1
                donors_read = 0
                
                # Read donor lines
                while i < len(lines) and donors_read < num_donors:
                    donor_line = lines[i].strip()
                    
                    # Skip comment lines within donor section
                    if not donor_line or donor_line.startswith(';') or donor_line.startswith('#'):
                        i += 1
                        continue
                    
                    donors.append(donor_line)
                    donors_read += 1
                    i += 1
                
                if donors_read == num_donors:
                    residues[res_name] = {
                        'num_donors': num_donors,
                        'donors': donors
                    }
                else:
                    print(f"Warning: Expected {num_donors} donors for {res_name}, found {donors_read}")
                
            except (ValueError, IndexError):
                i += 1
        else:
            i += 1
    
    return residues

def extract_rtp_mappings(rtp_file):
    correspondences = []
    
    with open(rtp_file, 'r') as f:
        content = f.read()
    
    # Find all [ RES ] sections with CHARMM names
    pattern = r'\[ (\w{3}) \]\s*(;.*?Original CHARMM name:\s*([\w\d]+))?'
    matches = re.findall(pattern, content, re.DOTALL)
    
    for match in matches:
        pdb_code = match[0].strip()
        if match[2]:  # CHARMM name found
            charmm_name = match[2].strip()
            correspondences.append((pdb_code, charmm_name))
        else:
            # Try to find CHARMM name in comments after the residue name
            res_pattern = rf'\[ {pdb_code} \](.*?)(?:\[|\Z)'
            res_match = re.search(res_pattern, content, re.DOTALL)
            if res_match:
                section = res_match.group(1)
                charmm_match = re.search(r'Original CHARMM name:\s*([\w\d]+)', section)
                if charmm_match:
                    correspondences.append((pdb_code, charmm_match.group(1).strip()))
    
    return correspondences

# ============================================================================
# HELPER FUNCTIONS FOR HDB VALIDATION AND FORMATTING
# ============================================================================

def validate_hdb_donor_line(donor_line: str) -> Tuple[bool, str]:
    """
    Validate HDB donor line format for GROMACS compatibility.
    Returns (is_valid, error_message)
    """
    parts = donor_line.split()
    
    if len(parts) < 4:
        return False, f"Too few fields: {donor_line}"
    
    try:
        donor_type = int(parts[0])
        num_atoms = int(parts[1])
    except ValueError:
        return False, f"Invalid type or count: {donor_line}"
    
    # Check if number of atoms matches what's expected
    if len(parts[2:]) != num_atoms:
        return False, f"Atom count mismatch: expected {num_atoms}, got {len(parts[2:])} in {donor_line}"
    
    # Validate based on donor type
    if donor_type == 1:  # OH group (hydrogen on oxygen)
        if num_atoms != 3:
            return False, f"Type 1 (OH) requires 3 atoms, got {num_atoms} in {donor_line}"
    elif donor_type == 2:  # NH group (hydrogen on nitrogen)
        if num_atoms != 3:
            return False, f"Type 2 (NH) requires 3 atoms, got {num_atoms} in {donor_line}"
    elif donor_type == 3:  # CH group (hydrogen on carbon)
        if num_atoms != 4:
            return False, f"Type 3 (CH) requires 4 atoms, got {num_atoms} in {donor_line}"
    elif donor_type == 4:  # NH2 group (hydrogen on NH2)
        if num_atoms != 4:
            return False, f"Type 4 (NH2) requires 4 atoms, got {num_atoms} in {donor_line}"
    else:
        # Unknown type, accept but warn
        print(f"  WARNING: Unknown donor type {donor_type} in line: {donor_line}")
    
    return True, ""

def fix_hdb_donor_line(donor_line: str) -> str:
    """
    Fix common HDB donor line formatting issues.
    Returns corrected line or original if no issues found.
    """
    parts = donor_line.split()
    
    if len(parts) < 4:
        return donor_line
    
    try:
        donor_type = int(parts[0])
        num_atoms = int(parts[1])
    except ValueError:
        return donor_line
    
    # Count actual atoms
    actual_atoms = len(parts[2:])
    
    # If atom count doesn't match, adjust
    if actual_atoms != num_atoms:
        # Rebuild line with correct count
        corrected_parts = [str(donor_type), str(actual_atoms)] + parts[2:]
        return ' '.join(corrected_parts)
    
    return donor_line

def format_hdb_donor_line(donor_line: str) -> str:
    """
    Format HDB donor line with proper spacing for GROMACS.
    Maintains the exact spacing shown in the example.
    """
    parts = donor_line.split()
    
    if len(parts) < 4:
        return donor_line
    
    # Format with proper spacing
    formatted = f"{parts[0]:<7} {parts[1]:<7} {parts[2]:<7}"
    
    # Add remaining atoms with proper spacing
    for i in range(3, len(parts)):
        formatted += f"{parts[i]:<7}"
    
    return formatted.rstrip()

# ============================================================================
# MAIN PROCESSOR
# ============================================================================

class PDBProcessor:
    def __init__(self, pdb_file: str, rtp_file: str, hdb_file: str, output_dir: str = None):
        self.pdb_file = pdb_file
        self.rtp_file = rtp_file
        self.hdb_file = hdb_file
        self.output_dir = output_dir
        
        self.pdb_residues = {}
        self.rtp_residues = {}
        self.hdb_residues = {}
        self.rtp_mappings = []
        
        self.variant_mapping = {}
        self.variants_info = {}
        self.variant_rtp_sections = {}
        self.variant_hdb_sections = {}
        
        self.report_lines = []
        self.atom_patterns = {}  # Store atom patterns for each variant type
    
    def load_pdb(self):
        print(f"Loading PDB: {self.pdb_file}")
        
        with open(self.pdb_file, 'r') as f:
            pdb_content = f.read()
        
        parser = CarbohydratePDBParser(pdb_content)
        
        for residue in parser.residues:
            if residue['name_pdb'] in STANDARD_AA:
                continue
            
            key = (residue['name_pdb'], residue['chain_carb'], residue['res_seq'])
            atom_names = {atom['name'] for atom in residue['atoms']}
            
            self.pdb_residues[key] = {
                'res_name': residue['name_pdb'],
                'chain': residue['chain_carb'],
                'res_seq': residue['res_seq'],
                'atoms': residue['atoms'],
                'atom_names': atom_names
            }
        
        print(f"Total carbohydrate residues loaded: {len(self.pdb_residues)}")
    
    def load_rtp(self):
        print(f"Loading RTP: {self.rtp_file}")
        
        self.rtp_residues = parse_rtp_file(self.rtp_file)
        self.rtp_mappings = extract_rtp_mappings(self.rtp_file)
        
        print(f"Total RTP residues loaded: {len(self.rtp_residues)}")
        print(f"Total RTP mappings found: {len(self.rtp_mappings)}")
        
        # Debug: show first few residues
        for i, (res_name, data) in enumerate(list(self.rtp_residues.items())[:3]):
            print(f"  {res_name}: {len(data['atoms'])} atoms")
    
    def load_hdb(self):
        print(f"Loading HDB: {self.hdb_file}")
        
        self.hdb_residues = parse_hdb_file(self.hdb_file)
        
        print(f"Total HDB residues loaded: {len(self.hdb_residues)}")
        
        # Debug: Show sample HDB data
        print("Sample HDB entries:")
        for i, (res_name, data) in enumerate(list(self.hdb_residues.items())[:3]):
            print(f"  {res_name}: {data['num_donors']} donors")
            for j, donor in enumerate(data['donors'][:2]):
                print(f"    {format_hdb_donor_line(donor)}")
    
    def identify_variants(self):
        print("\nIdentifying residue variants...")
        
        # Group residues by type
        residues_by_type = defaultdict(list)
        for key, residue in self.pdb_residues.items():
            residues_by_type[residue['res_name']].append((key, residue))
        
        # Process each residue type
        for res_name, residues in residues_by_type.items():
            self._process_residue_type(res_name, residues)
    
    def _process_residue_type(self, res_name: str, residues: List[Tuple[Tuple, dict]]):
        # Find RTP data for this residue
        rtp_residue = None
        charmm_name = None
        
        # Check if residue exists in RTP
        if res_name in self.rtp_residues:
            rtp_residue = self.rtp_residues[res_name]
            charmm_name = rtp_residue.get('name_charmm', '')
        else:
            # Try to find by CHARMM name
            for pdb_code, charmm in self.rtp_mappings:
                if pdb_code == res_name and charmm in self.rtp_residues:
                    rtp_residue = self.rtp_residues[charmm]
                    charmm_name = charmm
                    break
        
        if not rtp_residue:
            self.report_lines.append(f"\nResidue {res_name}: Not found in RTP, skipping...")
            print(f"  WARNING: No RTP data found for {res_name}")
            return
        
        rtp_atoms = {atom['name'] for atom in rtp_residue['atoms']}
        
        # Store patterns for this residue type
        if res_name not in self.atom_patterns:
            self.atom_patterns[res_name] = {}
        
        variant_counter = 1
        
        for (orig_res_name, chain, seq), residue in residues:
            pdb_atom_names = residue['atom_names']
            
            # Check if this pattern already exists
            found_variant = None
            for variant_name, pattern_atoms in self.atom_patterns[res_name].items():
                if pdb_atom_names == pattern_atoms:
                    found_variant = variant_name
                    break
            
            # If not found, create new variant
            if found_variant is None:
                # Create variant name (first 2 letters + number)
                variant_name = f"{res_name[:2]}{variant_counter}"
                variant_counter += 1
                
                # Store this atom pattern
                self.atom_patterns[res_name][variant_name] = pdb_atom_names
                found_variant = variant_name
                
                # Identify missing and extra atoms
                missing_atoms = rtp_atoms - pdb_atom_names
                extra_atoms = pdb_atom_names - rtp_atoms
                
                # Store variant info
                self.variants_info[found_variant] = {
                    'original': res_name,
                    'charmm_name': charmm_name,
                    'atoms_present': sorted(pdb_atom_names),
                    'atoms_missing': sorted(missing_atoms),
                    'atoms_extra': sorted(extra_atoms),
                    'count': 0,
                    'rtp_atoms_count': len(rtp_atoms),
                    'pdb_atoms_count': len(pdb_atom_names)
                }
            
            # Update variant count and mapping
            self.variants_info[found_variant]['count'] += 1
            self.variant_mapping[(orig_res_name, chain, seq)] = found_variant
            
            # Add to report
            self.report_lines.append(
                f"Residue {res_name} {chain}:{seq} -> Variant {found_variant}"
            )
        
        # Create RTP and HDB for each variant of this type
        for variant_name, atom_set in self.atom_patterns[res_name].items():
            self._create_variant_rtp(res_name, variant_name, atom_set, rtp_residue)
            self._create_variant_hdb(res_name, variant_name, atom_set, charmm_name)
    
    def _adjust_charges(self, variant_atoms, missing_atoms, rtp_residue):
        """Adjust charges by adding missing charges to C1"""
        if not missing_atoms:
            return variant_atoms
        
        # Calculate total missing charge
        missing_charge = 0.0
        for atom in rtp_residue['atoms']:
            if atom['name'] in missing_atoms:
                missing_charge += atom['charge']
        
        if abs(missing_charge) > 0.001:
            # Find C1 atom
            c1_index = -1
            for i, atom in enumerate(variant_atoms):
                if atom['name'] == 'C1':
                    c1_index = i
                    break
            
            if c1_index >= 0:
                old_charge = variant_atoms[c1_index]['charge']
                variant_atoms[c1_index]['charge'] += missing_charge
                print(f"  Adjusted C1 charge for missing atoms: {old_charge:.4f} -> {variant_atoms[c1_index]['charge']:.4f}")
            else:
                # If no C1, add to first atom
                if variant_atoms:
                    old_charge = variant_atoms[0]['charge']
                    variant_atoms[0]['charge'] += missing_charge
                    print(f"  Adjusted {variant_atoms[0]['name']} charge: {old_charge:.4f} -> {variant_atoms[0]['charge']:.4f}")
        
        return variant_atoms
    
    def _create_variant_rtp(self, original_name: str, variant_name: str, atom_set: Set[str], rtp_residue: dict):
        """Create RTP section for variant"""
        # Filter atoms that are present
        variant_atoms = []
        missing_atoms = []
        
        for rtp_atom in rtp_residue['atoms']:
            if rtp_atom['name'] in atom_set:
                variant_atoms.append(rtp_atom.copy())
            else:
                missing_atoms.append(rtp_atom['name'])
        
        # Adjust charges
        variant_atoms = self._adjust_charges(variant_atoms, missing_atoms, rtp_residue)
        
        # Filter bonds
        variant_bonds = []
        for atom1, atom2 in rtp_residue['bonds']:
            if atom1 in atom_set and atom2 in atom_set:
                variant_bonds.append((atom1, atom2))
        
        # Filter impropers
        variant_impropers = []
        for improper in rtp_residue['impropers']:
            if all(atom in atom_set for atom in improper):
                variant_impropers.append(improper)
        
        # Create comments with all original information
        comments = []
        comments.append(f"; Variant of {original_name}")
        if rtp_residue.get('name_charmm'):
            comments.append(f"; Original CHARMM name: {rtp_residue['name_charmm']}")
        if rtp_residue.get('name_iupac'):
            comments.append(f"; IUPAC name: {rtp_residue['name_iupac']}")
        
        # Add atom information
        comments.append(f"; Atoms in PDB: {', '.join(sorted(atom_set))}")
        if missing_atoms:
            comments.append(f"; Missing atoms: {', '.join(sorted(missing_atoms))}")
        
        # Calculate total charge
        total_charge = sum(atom['charge'] for atom in variant_atoms)
        comments.append(f"; Total charge: {total_charge:.6f}")
        
        # Store variant RTP
        self.variant_rtp_sections[variant_name] = {
            'original': original_name,
            'comments': comments,
            'atoms': variant_atoms,
            'bonds': variant_bonds,
            'impropers': variant_impropers,
            'missing_atoms': missing_atoms,
            'total_charge': total_charge
        }
        
        print(f"  Created RTP for {variant_name}: {len(variant_atoms)} atoms")
    
    def _create_variant_hdb(self, original_name: str, variant_name: str, atom_set: Set[str], charmm_name: str = None):
        """Create HDB section for variant"""
        hdb_data = None
        
        # Try to find HDB data
        if original_name in self.hdb_residues:
            hdb_data = self.hdb_residues[original_name]
        elif charmm_name and charmm_name in self.hdb_residues:
            hdb_data = self.hdb_residues[charmm_name]
        
        if not hdb_data:
            print(f"  WARNING: No HDB data found for {original_name} (CHARMM: {charmm_name})")
            return
        
        variant_donors = []
        invalid_donors = []
        
        for donor_line in hdb_data['donors']:
            # First, fix any formatting issues
            fixed_line = fix_hdb_donor_line(donor_line)
            
            # Validate the line
            is_valid, error_msg = validate_hdb_donor_line(fixed_line)
            
            if not is_valid:
                invalid_donors.append((donor_line, error_msg))
                continue
            
            # Parse the line to check atom presence
            parts = fixed_line.split()
            try:
                num_atoms = int(parts[1])
                atoms = parts[2:2 + num_atoms]
                
                # Check if all atoms are present in this variant
                if all(atom in atom_set for atom in atoms):
                    variant_donors.append(fixed_line)
                else:
                    # Some atoms missing, skip this donor
                    missing_in_variant = [atom for atom in atoms if atom not in atom_set]
                    print(f"    Skipping donor {donor_line}: missing atoms {missing_in_variant}")
            except (ValueError, IndexError):
                invalid_donors.append((donor_line, "Parsing error"))
        
        if invalid_donors:
            print(f"  Found {len(invalid_donors)} invalid HDB donors for {original_name}:")
            for donor, error in invalid_donors[:3]:  # Show only first 3
                print(f"    {donor} - {error}")
        
        if variant_donors:
            self.variant_hdb_sections[variant_name] = {
                'original': original_name,
                'num_donors': len(variant_donors),
                'donors': variant_donors
            }
            print(f"  Created HDB for {variant_name}: {len(variant_donors)} valid donors")
        else:
            print(f"  No valid HDB donors for {variant_name}")
    
    def generate_modified_pdb(self, output_pdb: str):
        print(f"\nGenerating modified PDB: {output_pdb}")
        
        modified_lines = []
        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    chain = line[21:22].strip()
                    try:
                        res_seq = int(line[22:26].strip())
                    except ValueError:
                        res_seq = 0
                    
                    key = (res_name, chain, res_seq)
                    
                    if key in self.variant_mapping:
                        variant_name = self.variant_mapping[key]
                        # Replace residue name
                        new_line = line[:17] + f"{variant_name:>3}" + line[20:]
                        modified_lines.append(new_line)
                    else:
                        modified_lines.append(line)
                else:
                    modified_lines.append(line)
        
        with open(output_pdb, 'w') as f:
            f.writelines(modified_lines)
        
        print(f"Modified PDB saved to: {output_pdb}")
    
    def generate_rtp_file(self, output_rtp: str):
        print(f"Generating RTP with variants: {output_rtp}")
        
        with open(output_rtp, 'w') as f:
            for variant_name in sorted(self.variant_rtp_sections.keys()):
                rtp_section = self.variant_rtp_sections[variant_name]
                
                # Write residue header
                f.write(f"[ {variant_name} ]\n")
                
                # Write comments
                for comment in rtp_section['comments']:
                    f.write(f"{comment}\n")
                
                # Write atoms with proper spacing
                f.write(" [ atoms ]\n")
                for atom in rtp_section['atoms']:
                    f.write(f"  {atom['name']:4} {atom['type']:7} {atom['charge']:7.4f} {atom['charge_group']}\n")
                
                # Write bonds if any
                if rtp_section['bonds']:
                    f.write(" [ bonds ]\n")
                    for atom1, atom2 in rtp_section['bonds']:
                        f.write(f"  {atom1:4} {atom2:4}\n")
                
                # Write impropers if any
                if rtp_section['impropers']:
                    f.write(" [ impropers ]\n")
                    for improper in rtp_section['impropers']:
                        f.write(f"  {' '.join(improper)}\n")
                
                f.write("\n")
        
        print(f"RTP with variants saved to: {output_rtp}")
    
    def generate_hdb_file(self, output_hdb: str):
        print(f"Generating HDB with variants: {output_hdb}")
        
        if not self.variant_hdb_sections:
            print("  WARNING: No HDB sections generated!")
            with open(output_hdb, 'w') as f:
                f.write("; HDB file generated by carbohydrate variant processor\n")
                f.write("; No valid hydrogen donors found for any variant\n")
            return
        
        with open(output_hdb, 'w') as f:
            f.write("; HDB file for carbohydrate variants - GROMACS compatible\n")
            f.write("; Generated from original HDB with missing atoms removed\n")
            f.write("; All donors validated for GROMACS compatibility\n")
            f.write("\n")
            
            for variant_name in sorted(self.variant_hdb_sections.keys()):
                hdb_section = self.variant_hdb_sections[variant_name]
                
                # Write header with proper spacing
                f.write(f"{variant_name:<5} {hdb_section['num_donors']}\n")
                for donor_line in hdb_section['donors']:
                    # Format with proper spacing
                    formatted_line = format_hdb_donor_line(donor_line)
                    f.write(f"{formatted_line}\n")
                f.write("\n")
        
        print(f"HDB with variants saved to: {output_hdb}")
        print(f"  Total HDB variants: {len(self.variant_hdb_sections)}")
        
        # Validation summary
        print("\nHDB Validation Summary:")
        for variant_name, hdb_section in self.variant_hdb_sections.items():
            print(f"  {variant_name}: {hdb_section['num_donors']} donors")
            for donor in hdb_section['donors'][:3]:  # Show first 3
                print(f"    {format_hdb_donor_line(donor)}")
    
    def generate_report(self, output_report: str):
        print(f"Generating report: {output_report}")
        
        with open(output_report, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CARBOHYDRATE VARIANT ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("SUMMARY BY RESIDUE TYPE:\n")
            f.write("-" * 40 + "\n")
            
            # Group variants by original residue
            variants_by_original = defaultdict(list)
            for variant_name, info in self.variants_info.items():
                variants_by_original[info['original']].append(variant_name)
            
            for original_name, variants in sorted(variants_by_original.items()):
                f.write(f"\n{original_name}:\n")
                total_residues = sum(self.variants_info[v]['count'] for v in variants)
                f.write(f"  Total residues: {total_residues}\n")
                
                if variants and self.variants_info[variants[0]]['charmm_name']:
                    f.write(f"  CHARMM name: {self.variants_info[variants[0]]['charmm_name']}\n")
                
                for variant in sorted(variants):
                    info = self.variants_info[variant]
                    f.write(f"  {variant}: {info['count']} residues\n")
                    f.write(f"    Atoms in PDB: {info['pdb_atoms_count']}\n")
                    
                    if info['atoms_missing']:
                        f.write(f"    Missing atoms: {', '.join(info['atoms_missing'])}\n")
                    if info['atoms_extra']:
                        f.write(f"    Extra atoms: {', '.join(info['atoms_extra'])}\n")
            
            f.write("\n\n" + "=" * 80 + "\n")
            f.write("RESIDUE MAPPING DETAILS:\n")
            f.write("=" * 80 + "\n")
            
            for (orig_name, chain, seq), variant in sorted(self.variant_mapping.items()):
                f.write(f"{orig_name} {chain}:{seq:4d} -> {variant}\n")
            
            f.write("\n\n" + "=" * 80 + "\n")
            f.write("STATISTICS:\n")
            f.write("=" * 80 + "\n")
            f.write(f"Total residues processed: {len(self.variant_mapping)}\n")
            f.write(f"Total variants created: {len(self.variants_info)}\n")
            f.write(f"RTP variants generated: {len(self.variant_rtp_sections)}\n")
            f.write(f"HDB variants generated: {len(self.variant_hdb_sections)}\n")
            
            f.write("\n\nCHARGE INFORMATION:\n")
            f.write("-" * 40 + "\n")
            for variant_name, rtp_section in sorted(self.variant_rtp_sections.items()):
                f.write(f"{variant_name}: Total charge = {rtp_section['total_charge']:.6f}\n")
            
            f.write("\n\nHDB VALIDATION:\n")
            f.write("-" * 40 + "\n")
            for variant_name, hdb_section in sorted(self.variant_hdb_sections.items()):
                f.write(f"{variant_name}: {hdb_section['num_donors']} valid hydrogen donors\n")
        
        print(f"Report saved to: {output_report}")
    
    def generate_json(self, output_json: str):
        print(f"Generating JSON: {output_json}")
        
        # Create comprehensive JSON data
        json_data = {
            'pdb_file': self.pdb_file,
            'rtp_file': self.rtp_file,
            'hdb_file': self.hdb_file,
            'variants': self.variants_info,
            'mapping': {
                f"{orig[0]}_{orig[1]}_{orig[2]}": variant 
                for orig, variant in self.variant_mapping.items()
            },
            'variant_rtp': {},
            'variant_hdb': {}
        }
        
        # Add RTP data
        for variant_name, section in self.variant_rtp_sections.items():
            json_data['variant_rtp'][variant_name] = {
                'original': section['original'],
                'atom_count': len(section['atoms']),
                'bond_count': len(section['bonds']),
                'improper_count': len(section['impropers']),
                'missing_atoms': section['missing_atoms'],
                'total_charge': section['total_charge'],
                'atoms': [
                    {
                        'name': atom['name'],
                        'type': atom['type'],
                        'charge': atom['charge'],
                        'charge_group': atom['charge_group']
                    }
                    for atom in section['atoms']
                ]
            }
        
        # Add HDB data
        for variant_name, section in self.variant_hdb_sections.items():
            json_data['variant_hdb'][variant_name] = {
                'original': section['original'],
                'num_donors': section['num_donors'],
                'donors': section['donors']
            }
        
        with open(output_json, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        print(f"JSON saved to: {output_json}")
    
    def run(self):
        print("=" * 80)
        print("CARBOHYDRATE VARIANT PROCESSOR - GROMACS COMPATIBLE")
        print("=" * 80)
        
        self.load_pdb()
        self.load_rtp()
        self.load_hdb()
        
        self.identify_variants()
        
        # Generate output files
        output_prefix = os.path.splitext(os.path.basename(self.pdb_file))[0]
        
        # Create output directory if specified
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
            output_path = os.path.join(self.output_dir, f"{output_prefix}_variants.pdb")
        else:
            output_path = f"{output_prefix}_variants.pdb"
        
        self.generate_modified_pdb(output_path)
        self.generate_rtp_file(output_path.replace('.pdb', '.rtp'))
        self.generate_hdb_file(output_path.replace('.pdb', '.hdb'))
        self.generate_report(output_path.replace('.pdb', '_report.txt'))
        self.generate_json(output_path.replace('.pdb', '_data.json'))
        
        print("\n" + "=" * 80)
        print("PROCESSING COMPLETE!")
        print("=" * 80)

def main():
    parser = argparse.ArgumentParser(
        description="Process PDB glycoproteins and create carbohydrate residue variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -p spike.pdb -r carb.rtp -d carb.hdb
  %(prog)s -p spike.pdb -r carb.rtp -d carb.hdb -o ./output
  %(prog)s --pdb spike_glycosylated_final_connected.pdb --rtp carb_redundance_removed.rtp --hdb carb_redundance_removed.hdb --output ./results
        """
    )
    
    parser.add_argument('-p', '--pdb', required=True,
                      help='Input PDB file')
    parser.add_argument('-r', '--rtp', required=True,
                      help='Input RTP file')
    parser.add_argument('-d', '--hdb', required=True,
                      help='Input HDB file')
    parser.add_argument('-o', '--output', default='.',
                      help='Output directory (default: current directory)')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Verbose output')
    
    args = parser.parse_args()
    
    # Validate input files
    for file_path in [args.pdb, args.rtp, args.hdb]:
        if not os.path.exists(file_path):
            print(f"ERROR: File not found: {file_path}")
            sys.exit(1)
    
    # Create processor and run
    processor = PDBProcessor(
        pdb_file=args.pdb,
        rtp_file=args.rtp,
        hdb_file=args.hdb,
        output_dir=args.output
    )
    
    processor.run()

# Main execution
if __name__ == "__main__":
    main()
