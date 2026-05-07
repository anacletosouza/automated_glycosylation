"""
Script to add CHARMM36 force field parameters to glycan JSON data.
Reads glycan_data.json and adds CHARMM36 parameters (charges, masses, epsilon, sigma).
Saves as glycan_data_charmm36.json.

Usage:
------
python add_charmm_params.py --input_json glycan_data.json --charmm_dir /path/to/charmm36.ff --output_json glycan_data_charmm36.json
"""

import os
import json
import argparse
from pathlib import Path

def parse_atomtypes(atp_file):
    """
    Parse atomtypes.atp file to get mass, charge, epsilon, sigma for each atom type.
    
    Parameters:
    -----------
    atp_file : str
        Path to atomtypes.atp file
    
    Returns:
    --------
    dict : Dictionary mapping atom types to their parameters
    """
    atomtypes = {}
    
    with open(atp_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(';') and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 7:
                    atom_type = parts[0]
                    try:
                        mass = float(parts[1])
                        charge = float(parts[2])
                        ptype = parts[3]
                        sigma = float(parts[4])
                        epsilon = float(parts[5])
                        
                        atomtypes[atom_type] = {
                            'mass': mass,
                            'charge': charge,
                            'ptype': ptype,
                            'sigma': sigma,
                            'epsilon': epsilon
                        }
                    except (ValueError, IndexError):
                        continue
    
    return atomtypes

def parse_ffnonbonded(itp_file):
    """
    Parse ffnonbonded.itp file to get atom type parameters.
    
    Parameters:
    -----------
    itp_file : str
        Path to ffnonbonded.itp file
    
    Returns:
    --------
    dict : Dictionary mapping atom names to atom types and parameters
    """
    atom_params = {}
    current_section = None
    
    with open(itp_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith(';') or line.startswith('#'):
                continue
            
            # Check for section headers
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                continue
            
            if current_section == 'atomtypes':
                parts = line.split()
                if len(parts) >= 6:
                    atom_type = parts[0]
                    try:
                        atom_params[atom_type] = {
                            'bond_type': parts[1] if len(parts) > 1 else '',
                            'mass': float(parts[2]),
                            'charge': float(parts[3]),
                            'ptype': parts[4] if len(parts) > 4 else 'A',
                            'sigma': float(parts[5]) if len(parts) > 5 else 0.1,
                            'epsilon': float(parts[6]) if len(parts) > 6 else 0.0
                        }
                    except (ValueError, IndexError):
                        continue
    
    return atom_params

def parse_rtp_files(charmm_dir):
    """
    Parse all .rtp files to get residue topologies with atom charges.
    
    Parameters:
    -----------
    charmm_dir : str
        Path to CHARMM36 force field directory
    
    Returns:
    --------
    dict : Dictionary mapping residue names to atom information
    """
    rtp_files = ['aminoacids.rtp', 'carb.rtp', 'cgenff.rtp', 'lipid.rtp', 'na.rtp']
    residues = {}
    
    for rtp_file in rtp_files:
        rtp_path = Path(charmm_dir) / rtp_file
        if not rtp_path.exists():
            continue
        
        print(f"  Parsing {rtp_file}...")
        current_residue = None
        in_atoms_section = False
        
        with open(rtp_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith(';') or line.startswith('#'):
                    continue
                
                # Check for residue header
                if line.startswith('[') and not line.startswith('[ atoms ]') and not line.startswith('[ bonds ]'):
                    current_residue = line[1:-1].strip()
                    residues[current_residue] = {'atoms': {}, 'bonds': []}
                    in_atoms_section = False
                
                # Start of atoms section
                elif line == '[ atoms ]':
                    in_atoms_section = True
                
                # Start of bonds section
                elif line == '[ bonds ]':
                    in_atoms_section = False
                
                # Parse atom lines in atoms section
                elif in_atoms_section and current_residue:
                    parts = line.split()
                    if len(parts) >= 4:
                        atom_name = parts[0]
                        atom_type = parts[1]
                        charge = float(parts[2])
                        
                        # Store atom information
                        residues[current_residue]['atoms'][atom_name] = {
                            'type': atom_type,
                            'charge': charge
                        }
                
                # Parse bond lines (optional)
                elif line.startswith('[') and 'bonds' in line.lower():
                    # Skip bond parsing for now
                    pass
    
    return residues

def get_atom_mapping():
    """
    Create a mapping from element names to default CHARMM36 atom types.
    This is used as fallback when exact atom types are not found.
    
    Returns:
    --------
    dict : Mapping from element to default atom type
    """
    # Default mappings for common elements
    return {
        'C': 'CTL2',    # Tetrahedral carbon (sp3)
        'CA': 'CTL2',   # Aromatic carbon
        'CB': 'CTL2',   # Aliphatic carbon
        'N': 'NH1',     # Amide nitrogen
        'O': 'OC',      # Carbonyl oxygen
        'OH': 'OH1',    # Hydroxyl oxygen
        'OG': 'OH1',    # Serine/Threonine oxygen
        'OD1': 'OC',    # Aspartic acid oxygen
        'OD2': 'OC',    # Aspartic acid oxygen
        'OE1': 'OC',    # Glutamic acid oxygen
        'OE2': 'OC',    # Glutamic acid oxygen
        'NZ': 'NH3',    # Lysine nitrogen
        'NE': 'NH2',    # Arginine nitrogen
        'NH1': 'NH2',   # Arginine nitrogen
        'NH2': 'NH2',   # Arginine nitrogen
        'ND2': 'NH2',   # Asparagine nitrogen
        'NE2': 'NH2',   # Glutamine nitrogen
        'SD': 'S',      # Sulfur
        'SG': 'S',      # Cysteine sulfur
        'H': 'HAL2',    # Hydrogen
        'HA': 'HAL2',   # Aliphatic hydrogen
        'HB': 'HAL2',   # Aliphatic hydrogen
        'HD': 'HAL2',   # Aliphatic hydrogen
        'HE': 'HAL2',   # Aliphatic hydrogen
        'HG': 'HAL2',   # Aliphatic hydrogen
        'HH': 'HAL2',   # Aliphatic hydrogen
        'HZ': 'HAL2',   # Aliphatic hydrogen
    }

def get_default_params(element, atom_name=""):
    """
    Get default CHARMM36 parameters based on element and atom name.
    
    Parameters:
    -----------
    element : str
        Element symbol
    atom_name : str
        Atom name
    
    Returns:
    --------
    dict : Default parameters
    """
    # Default parameters for common elements
    defaults = {
        'C': {'mass': 12.011, 'charge': 0.0, 'sigma': 0.355, 'epsilon': 0.2929},
        'N': {'mass': 14.007, 'charge': -0.5, 'sigma': 0.325, 'epsilon': 0.7113},
        'O': {'mass': 15.999, 'charge': -0.5, 'sigma': 0.296, 'epsilon': 0.8786},
        'S': {'mass': 32.06, 'charge': 0.0, 'sigma': 0.356, 'epsilon': 1.046},
        'H': {'mass': 1.008, 'charge': 0.25, 'sigma': 0.242, 'epsilon': 0.126},
        'P': {'mass': 30.974, 'charge': 1.5, 'sigma': 0.374, 'epsilon': 0.8368},
    }
    
    # Special cases based on atom name
    special_cases = {
        # Carbonyl carbons
        'C': {'mass': 12.011, 'charge': 0.5, 'sigma': 0.355, 'epsilon': 0.2929} if any(x in atom_name for x in ['C=', 'CO', 'OC']) 
              else defaults.get('C', {'mass': 12.011, 'charge': 0.0, 'sigma': 0.355, 'epsilon': 0.2929}),
        # Carbonyl oxygens
        'O': {'mass': 15.999, 'charge': -0.5, 'sigma': 0.296, 'epsilon': 0.8786} if any(x in atom_name for x in ['O=', 'OC', 'CO']) 
              else defaults.get('O', {'mass': 15.999, 'charge': -0.5, 'sigma': 0.296, 'epsilon': 0.8786}),
        # Amide nitrogens
        'N': {'mass': 14.007, 'charge': -0.5, 'sigma': 0.325, 'epsilon': 0.7113} if any(x in atom_name for x in ['N=', 'NH', 'ND', 'NE']) 
              else defaults.get('N', {'mass': 14.007, 'charge': -0.5, 'sigma': 0.325, 'epsilon': 0.7113}),
    }
    
    if element in special_cases:
        return special_cases[element]
    elif element in defaults:
        return defaults[element]
    else:
        # Generic defaults for unknown elements
        return {'mass': 1.0, 'charge': 0.0, 'sigma': 0.1, 'epsilon': 0.0}

def find_atom_type_for_atom(atom_info, residues_data, atom_mapping):
    """
    Find the CHARMM36 atom type for a given atom.
    
    Parameters:
    -----------
    atom_info : dict
        Atom information from JSON
    residues_data : dict
        Residue topology data from .rtp files
    atom_mapping : dict
        Element to atom type mapping
    
    Returns:
    --------
    str : Atom type
    """
    residue_name = atom_info['residue_name']
    atom_name = atom_info['atom_name'].strip()
    element = atom_info['element']
    
    # Try to find in residues data
    if residue_name in residues_data:
        residue_atoms = residues_data[residue_name]['atoms']
        
        # Try exact match
        if atom_name in residue_atoms:
            return residue_atoms[atom_name]['type']
        
        # Try removing leading/trailing spaces and numbers
        atom_name_clean = atom_name
        while atom_name_clean and atom_name_clean[-1].isdigit():
            atom_name_clean = atom_name_clean[:-1]
        
        if atom_name_clean in residue_atoms:
            return residue_atoms[atom_name_clean]['type']
        
        # Try common variations
        common_names = {
            'C1': ['C1', 'C1*', 'C1\''],
            'C2': ['C2', 'C2*', 'C2\''],
            'C3': ['C3', 'C3*', 'C3\''],
            'C4': ['C4', 'C4*', 'C4\''],
            'C5': ['C5', 'C5*', 'C5\''],
            'C6': ['C6', 'C6*', 'C6\''],
            'C7': ['C7', 'C7*', 'C7\''],
            'C8': ['C8', 'C8*', 'C8\''],
            'C9': ['C9', 'C9*', 'C9\''],
            'N': ['N', 'N*', 'N\'', 'N1', 'N2'],
            'O1': ['O1', 'O1*', 'O1\'', 'O'],
            'O2': ['O2', 'O2*', 'O2\''],
            'O3': ['O3', 'O3*', 'O3\''],
            'O4': ['O4', 'O4*', 'O4\''],
            'O5': ['O5', 'O5*', 'O5\''],
            'O6': ['O6', 'O6*', 'O6\''],
        }
        
        if atom_name in common_names:
            for variant in common_names[atom_name]:
                if variant in residue_atoms:
                    return residue_atoms[variant]['type']
    
    # Fallback to element-based mapping
    if atom_name in atom_mapping:
        return atom_mapping[atom_name]
    
    # Use element symbol
    if element in atom_mapping:
        return atom_mapping[element]
    
    # Default based on element
    if element == 'C':
        return 'CTL2'
    elif element == 'N':
        return 'NH1'
    elif element == 'O':
        return 'OC'
    elif element == 'S':
        return 'S'
    elif element == 'H':
        return 'HAL2'
    else:
        return 'CTL2'  # Default carbon type

def get_atom_charge(atom_info, residues_data):
    """
    Get the CHARMM36 charge for a given atom.
    
    Parameters:
    -----------
    atom_info : dict
        Atom information from JSON
    residues_data : dict
        Residue topology data from .rtp files
    
    Returns:
    --------
    float : Atom charge
    """
    residue_name = atom_info['residue_name']
    atom_name = atom_info['atom_name'].strip()
    
    # Try to find in residues data
    if residue_name in residues_data:
        residue_atoms = residues_data[residue_name]['atoms']
        
        # Try exact match
        if atom_name in residue_atoms:
            return residue_atoms[atom_name]['charge']
        
        # Try removing numbers
        atom_name_clean = atom_name
        while atom_name_clean and atom_name_clean[-1].isdigit():
            atom_name_clean = atom_name_clean[:-1]
        
        if atom_name_clean in residue_atoms:
            return residue_atoms[atom_name_clean]['charge']
    
    # Default charge based on element
    element = atom_info['element']
    defaults = {
        'C': 0.0,
        'N': -0.5,
        'O': -0.5,
        'S': 0.0,
        'H': 0.25,
        'P': 1.5,
    }
    
    return defaults.get(element, 0.0)

def add_charmm_parameters(data, charmm_dir):
    """
    Add CHARMM36 force field parameters to all atoms in the data.
    
    Parameters:
    -----------
    data : dict
        JSON data structure
    charmm_dir : str
        Path to CHARMM36 force field directory
    
    Returns:
    --------
    dict : Updated data with CHARMM36 parameters
    """
    print("Loading CHARMM36 force field parameters...")
    
    # Parse force field files
    atp_file = Path(charmm_dir) / 'atomtypes.atp'
    itp_file = Path(charmm_dir) / 'ffnonbonded.itp'
    
    # Get atom type parameters
    if atp_file.exists():
        atomtypes = parse_atomtypes(atp_file)
        print(f"  Loaded {len(atomtypes)} atom types from atomtypes.atp")
    else:
        print(f"  Warning: atomtypes.atp not found at {atp_file}")
        atomtypes = {}
    
    # Get additional parameters from ffnonbonded.itp
    if itp_file.exists():
        ff_params = parse_ffnonbonded(itp_file)
        print(f"  Loaded {len(ff_params)} atom types from ffnonbonded.itp")
        
        # Merge with atomtypes
        for atom_type, params in ff_params.items():
            if atom_type not in atomtypes:
                atomtypes[atom_type] = params
    else:
        print(f"  Warning: ffnonbonded.itp not found at {itp_file}")
    
    # Parse residue topologies
    residues_data = parse_rtp_files(charmm_dir)
    print(f"  Loaded {len(residues_data)} residue topologies")
    
    # Create atom mapping for fallback
    atom_mapping = get_atom_mapping()
    
    # Add parameters to protein atoms
    print("Adding parameters to protein atoms...")
    protein_atoms_with_params = []
    
    for atom in data['protein']:
        atom_copy = atom.copy()
        
        # Find atom type
        atom_type = find_atom_type_for_atom(atom, residues_data, atom_mapping)
        atom_copy['charmm_type'] = atom_type
        
        # Get charge
        charge = get_atom_charge(atom, residues_data)
        atom_copy['charmm_charge'] = charge
        
        # Get other parameters
        if atom_type in atomtypes:
            params = atomtypes[atom_type]
            atom_copy['charmm_mass'] = params.get('mass', 1.0)
            atom_copy['charmm_sigma'] = params.get('sigma', 0.1)  # nm
            atom_copy['charmm_epsilon'] = params.get('epsilon', 0.0)  # kJ/mol
            
            # Override charge if different
            if 'charge' in params and params['charge'] != 0.0:
                atom_copy['charmm_charge'] = params['charge']
        else:
            # Use defaults based on element
            defaults = get_default_params(atom['element'], atom['atom_name'])
            atom_copy['charmm_mass'] = defaults['mass']
            atom_copy['charmm_sigma'] = defaults['sigma']
            atom_copy['charmm_epsilon'] = defaults['epsilon']
        
        protein_atoms_with_params.append(atom_copy)
    
    print(f"  Processed {len(protein_atoms_with_params)} protein atoms")
    
    # Add parameters to glycan atoms
    print("Adding parameters to glycan atoms...")
    glycans_with_params = {}
    
    for glycan_id, glycan_data in data['glycans'].items():
        glycan_copy = glycan_data.copy()
        atoms_with_params = []
        
        for atom in glycan_data['atoms']:
            atom_copy = atom.copy()
            
            # Find atom type
            atom_type = find_atom_type_for_atom(atom, residues_data, atom_mapping)
            atom_copy['charmm_type'] = atom_type
            
            # Get charge
            charge = get_atom_charge(atom, residues_data)
            atom_copy['charmm_charge'] = charge
            
            # Get other parameters
            if atom_type in atomtypes:
                params = atomtypes[atom_type]
                atom_copy['charmm_mass'] = params.get('mass', 1.0)
                atom_copy['charmm_sigma'] = params.get('sigma', 0.1)
                atom_copy['charmm_epsilon'] = params.get('epsilon', 0.0)
                
                # Override charge if different
                if 'charge' in params and params['charge'] != 0.0:
                    atom_copy['charmm_charge'] = params['charge']
            else:
                # Use defaults based on element
                defaults = get_default_params(atom['element'], atom['atom_name'])
                atom_copy['charmm_mass'] = defaults['mass']
                atom_copy['charmm_sigma'] = defaults['sigma']
                atom_copy['charmm_epsilon'] = defaults['epsilon']
            
            atoms_with_params.append(atom_copy)
        
        glycan_copy['atoms'] = atoms_with_params
        glycans_with_params[glycan_id] = glycan_copy
    
    print(f"  Processed {len(glycans_with_params)} glycans")
    
    # Add parameters to linkage atoms
    print("Adding parameters to linkage atoms...")
    linkages_with_params = []
    
    for linkage in data['linkages']:
        linkage_copy = linkage.copy()
        
        # Update protein atom in linkage
        if 'protein_atom_complete' in linkage:
            protein_atom = linkage['protein_atom_complete'].copy()
            atom_type = find_atom_type_for_atom(protein_atom, residues_data, atom_mapping)
            protein_atom['charmm_type'] = atom_type
            protein_atom['charmm_charge'] = get_atom_charge(protein_atom, residues_data)
            
            if atom_type in atomtypes:
                params = atomtypes[atom_type]
                protein_atom['charmm_mass'] = params.get('mass', 1.0)
                protein_atom['charmm_sigma'] = params.get('sigma', 0.1)
                protein_atom['charmm_epsilon'] = params.get('epsilon', 0.0)
            else:
                defaults = get_default_params(protein_atom['element'], protein_atom['atom_name'])
                protein_atom['charmm_mass'] = defaults['mass']
                protein_atom['charmm_sigma'] = defaults['sigma']
                protein_atom['charmm_epsilon'] = defaults['epsilon']
            
            linkage_copy['protein_atom_complete'] = protein_atom
        
        # Update glycan atom in linkage
        if 'glycan_atom_complete' in linkage:
            glycan_atom = linkage['glycan_atom_complete'].copy()
            atom_type = find_atom_type_for_atom(glycan_atom, residues_data, atom_mapping)
            glycan_atom['charmm_type'] = atom_type
            glycan_atom['charmm_charge'] = get_atom_charge(glycan_atom, residues_data)
            
            if atom_type in atomtypes:
                params = atomtypes[atom_type]
                glycan_atom['charmm_mass'] = params.get('mass', 1.0)
                glycan_atom['charmm_sigma'] = params.get('sigma', 0.1)
                glycan_atom['charmm_epsilon'] = params.get('epsilon', 0.0)
            else:
                defaults = get_default_params(glycan_atom['element'], glycan_atom['atom_name'])
                glycan_atom['charmm_mass'] = defaults['mass']
                glycan_atom['charmm_sigma'] = defaults['sigma']
                glycan_atom['charmm_epsilon'] = defaults['epsilon']
            
            linkage_copy['glycan_atom_complete'] = glycan_atom
        
        linkages_with_params.append(linkage_copy)
    
    print(f"  Processed {len(linkages_with_params)} linkages")
    
    # Create updated data structure
    updated_data = {
        'metadata': {
            **data['metadata'],
            'force_field': 'CHARMM36',
            'charmm_dir': str(charmm_dir),
            'total_atom_types_loaded': len(atomtypes),
            'total_residues_loaded': len(residues_data)
        },
        'protein': protein_atoms_with_params,
        'glycans': glycans_with_params,
        'linkages': linkages_with_params,
        'other_lines': data['other_lines'],
        'force_field_info': {
            'atomtypes_loaded': len(atomtypes),
            'residues_loaded': len(residues_data),
            'default_parameters_used': len(atomtypes) == 0
        }
    }
    
    return updated_data

def main():
    """Main function to add CHARMM36 parameters to JSON data."""
    parser = argparse.ArgumentParser(description="Add CHARMM36 force field parameters to glycan JSON data")
    parser.add_argument("--input_json", required=True, help="Path to input JSON file (glycan_data.json)")
    parser.add_argument("--charmm_dir", required=True, help="Path to CHARMM36 force field directory")
    parser.add_argument("--output_json", required=True, help="Path to output JSON file with CHARMM36 parameters")
    
    args = parser.parse_args()
    
    input_json = args.input_json
    charmm_dir = args.charmm_dir
    output_json = args.output_json
    
    print("=" * 70)
    print("ADDING CHARMM36 PARAMETERS TO GLYCAN DATA")
    print("=" * 70)
    
    # Check if CHARMM directory exists
    charmm_path = Path(charmm_dir)
    if not charmm_path.exists():
        print(f"Error: CHARMM36 directory not found at {charmm_dir}")
        return
    
    print(f"CHARMM36 directory: {charmm_dir}")
    
    # Load input JSON
    print(f"\n1. LOADING INPUT JSON")
    print("-" * 40)
    with open(input_json, 'r') as f:
        data = json.load(f)
    
    print(f"  Loaded data from {input_json}")
    print(f"  Protein atoms: {len(data['protein'])}")
    print(f"  Glycans: {len(data['glycans'])}")
    print(f"  Linkages: {len(data['linkages'])}")
    
    # Add CHARMM36 parameters
    print(f"\n2. ADDING CHARMM36 PARAMETERS")
    print("-" * 40)
    updated_data = add_charmm_parameters(data, charmm_dir)
    
    # Save updated JSON
    print(f"\n3. SAVING UPDATED JSON")
    print("-" * 40)
    with open(output_json, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    print(f"  Saved to {output_json}")
    
    # Print summary
    print(f"\n" + "=" * 70)
    print("PROCESSING COMPLETE!")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  - Input JSON: {input_json}")
    print(f"  - Output JSON: {output_json}")
    print(f"  - CHARMM36 directory: {charmm_dir}")
    print(f"  - Total protein atoms with parameters: {len(updated_data['protein'])}")
    print(f"  - Total glycans with parameters: {len(updated_data['glycans'])}")
    print(f"  - Total linkages with parameters: {len(updated_data['linkages'])}")
    print(f"  - Atom types loaded: {updated_data['metadata']['total_atom_types_loaded']}")
    print(f"  - Residue topologies loaded: {updated_data['metadata']['total_residues_loaded']}")
    
    # Show example of atom with CHARMM parameters
    if updated_data['protein']:
        sample_atom = updated_data['protein'][0]
        print(f"\nExample atom with CHARMM36 parameters:")
        print(f"  Residue: {sample_atom['residue_name']}{sample_atom['residue_number']}")
        print(f"  Atom: {sample_atom['atom_name']} ({sample_atom['element']})")
        print(f"  CHARMM type: {sample_atom.get('charmm_type', 'N/A')}")
        print(f"  Charge: {sample_atom.get('charmm_charge', 'N/A'):.4f} e")
        print(f"  Mass: {sample_atom.get('charmm_mass', 'N/A'):.4f} Da")
        print(f"  Sigma: {sample_atom.get('charmm_sigma', 'N/A'):.4f} nm")
        print(f"  Epsilon: {sample_atom.get('charmm_epsilon', 'N/A'):.4f} kJ/mol")

if __name__ == "__main__":
    main()
