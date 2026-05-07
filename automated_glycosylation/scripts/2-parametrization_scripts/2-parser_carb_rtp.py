import re

def parse_rtp_file(file_path):
    """
    Parses the .rtp file and extracts residue information into a structured dictionary.
    
    Args:
        file_path: Path to the .rtp file
    
    Returns:
        Dictionary with residue names as keys and structured data as values
    """
    residues = {}
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Split content by residue sections
    lines = content.split('\n')
    i = 0
    total_lines = len(lines)
    
    while i < total_lines:
        line = lines[i].strip()
        
        # Look for residue definitions (e.g., [ A2UDM ])
        if line.startswith('[') and line.endswith(']'):
            # Remove brackets and spaces
            res_name = line[1:-1].strip()
            
            # Skip bondedtypes section
            if res_name.upper() == 'BONDEDTYPES':
                i += 1
                continue
            
            # Initialize residue data structure
            residue_data = {
                'name_charmm': res_name,
                'name_iupac': '',
                'atoms': [],
                'bonds': [],
                'impropers': [],
                'comments': {}
            }
            
            # Read the residue description (next line starting with ';')
            if i + 1 < total_lines:
                desc_line = lines[i+1].strip()
                if desc_line.startswith(';'):
                    residue_data['name_iupac'] = desc_line[1:].strip()
                    i += 1
            
            # Look for various sections
            i += 1
            while i < total_lines:
                current_line = lines[i].strip()
                
                # End of this residue (start of next residue or end of file)
                if current_line.startswith('[') and current_line.endswith(']'):
                    # Check if it's a new residue or just a section header
                    section_match = re.match(r'\[\s*(atoms|bonds|impropers)\s*\]', current_line, re.IGNORECASE)
                    if not section_match:
                        break
                
                # [ atoms ] section
                if current_line.lower() == '[ atoms ]':
                    i += 1
                    while i < total_lines and lines[i].strip() and not lines[i].strip().startswith('['):
                        atom_line = lines[i].strip()
                        if atom_line and not atom_line.startswith(';'):
                            parts = atom_line.split()
                            if len(parts) >= 4:
                                atom_info = {
                                    'name': parts[0],
                                    'type': parts[1],
                                    'charge': float(parts[2]),
                                    'charge_group': int(parts[3])
                                }
                                residue_data['atoms'].append(atom_info)
                        i += 1
                    continue
                
                # [ bonds ] section
                elif current_line.lower() == '[ bonds ]':
                    i += 1
                    while i < total_lines and lines[i].strip() and not lines[i].strip().startswith('['):
                        bond_line = lines[i].strip()
                        if bond_line and not bond_line.startswith(';'):
                            parts = bond_line.split()
                            # Process pairs of atoms
                            for j in range(0, len(parts)-1, 2):
                                if j+1 < len(parts):
                                    residue_data['bonds'].append((parts[j], parts[j+1]))
                        i += 1
                    continue
                
                # [ impropers ] section
                elif current_line.lower() == '[ impropers ]':
                    i += 1
                    while i < total_lines and lines[i].strip() and not lines[i].strip().startswith('['):
                        improper_line = lines[i].strip()
                        if improper_line and not improper_line.startswith(';'):
                            parts = improper_line.split()
                            if len(parts) >= 4:
                                residue_data['impropers'].append(parts[:4])
                        i += 1
                    continue
                
                i += 1
            
            # Only add residue if it has atoms (to avoid empty entries)
            if residue_data['atoms']:
                residues[res_name] = residue_data
        
        else:
            i += 1
    
    return residues


def display_residue_info(residues_dict):
    """
    Display information about parsed residues.
    
    Args:
        residues_dict: Dictionary of parsed residues
    """
    print(f"Total residues parsed: {len(residues_dict)}")
    print("\nResidue list:")
    print("-" * 80)
    
    # Sort residues by name for consistent display
    sorted_residues = sorted(residues_dict.items(), key=lambda x: x[0])
    
    for res_name, res_data in sorted_residues[:10]:  # Show first 10
        print(f"{res_name:10} : {res_data['name_iupac'][:60]}...")
        print(f"           Atoms: {len(res_data['atoms']):3d} | "
              f"Bonds: {len(res_data['bonds']):3d} | "
              f"Impropers: {len(res_data['impropers']):2d}")
    
    if len(sorted_residues) > 10:
        print(f"... and {len(sorted_residues) - 10} more residues")
    
    print("\n" + "="*80)
    

def get_residue_details(residues_dict, residue_name):
    """
    Get detailed information about a specific residue.
    
    Args:
        residues_dict: Dictionary of parsed residues
        residue_name: Name of the residue to query
    
    Returns:
        Formatted string with residue details
    """
    if residue_name not in residues_dict:
        return f"Residue '{residue_name}' not found."
    
    res_data = residues_dict[residue_name]
    
    output = []
    output.append(f"\n{'='*80}")
    output.append(f"RESIDUE: {residue_name}")
    output.append(f"IUPAC Name: {res_data['name_iupac']}")
    output.append(f"{'='*80}")
    
    # Atoms section
    output.append(f"\nATOMS ({len(res_data['atoms'])}):")
    output.append("-" * 80)
    output.append(f"{'Name':<8} {'Type':<10} {'Charge':>8} {'Charge Group':>12}")
    output.append("-" * 80)
    
    for atom in res_data['atoms'][:20]:  # Show first 20 atoms
        output.append(f"{atom['name']:<8} {atom['type']:<10} {atom['charge']:>8.3f} {atom['charge_group']:>12}")
    
    if len(res_data['atoms']) > 20:
        output.append(f"... and {len(res_data['atoms']) - 20} more atoms")
    
    # Bonds section
    output.append(f"\nBONDS ({len(res_data['bonds'])}):")
    output.append("-" * 80)
    bonds_per_line = 5
    bonds_displayed = 0
    
    for i in range(0, min(25, len(res_data['bonds']))):  # Show up to 25 bonds
        if i % bonds_per_line == 0:
            if i > 0:
                output.append(current_line)
            current_line = "  "
        bond = res_data['bonds'][i]
        current_line += f"{bond[0]}-{bond[1]:<12}"
        bonds_displayed += 1
    
    if current_line.strip():
        output.append(current_line)
    
    if len(res_data['bonds']) > 25:
        output.append(f"... and {len(res_data['bonds']) - 25} more bonds")
    
    # Impropers section
    if res_data['impropers']:
        output.append(f"\nIMPROPERS ({len(res_data['impropers'])}):")
        output.append("-" * 80)
        for improper in res_data['impropers'][:10]:  # Show first 10 impropers
            output.append(f"  {' - '.join(improper)}")
        
        if len(res_data['impropers']) > 10:
            output.append(f"... and {len(res_data['impropers']) - 10} more impropers")
    
    output.append(f"{'='*80}")
    
    return "\n".join(output)


def export_to_pkl(residues_dict, output_file):
    """
    Export residues dictionary to PKL file.
    
    Args:
        residues_dict: Dictionary of parsed residues
        output_file: Output PKL file path
    """
    import pickle

    with open(output_file, 'wb') as f:
        pickle.dump(residues_dict, f)

    print(f"\nData exported to PKL file: {output_file}")


import sys

def main():
    # Use user-provided file path if given
    if len(sys.argv) > 1:
        carb_rtp_path = sys.argv[1]
    else:
        # default path if no argument is provided
        carb_rtp_path = "/home/anacleto/Desktop/testing_glycosilation/CHARMM_files/Delta/charmm36/carb.rtp"
    
    try:
        print(f"Parsing carb.rtp file at: {carb_rtp_path}")
        residues_dict = parse_rtp_file(carb_rtp_path)
        
        if not residues_dict:
            print("No residues found in the file.")
            return
        
        # Display summary
        display_residue_info(residues_dict)
        
        # Example: Get details for specific residues
        example_residues = ['A2UDM', 'ABAC', 'BOG', 'BOM']
        for res_name in example_residues:
            if res_name in residues_dict:
                print(get_residue_details(residues_dict, res_name))
                break  # Just show first found example
        
        export_to_pkl(residues_dict, "carb_residues.pkl")
        
        # Demonstrate dictionary access
        print("\n" + "="*80)
        print("DICTIONARY ACCESS EXAMPLES:")
        print("="*80)
        
        if 'BOG' in residues_dict:
            bog_data = residues_dict['BOG']
            print(f"\nAccessing BOG data directly:")
            print(f"  Name: {bog_data['name_charmm']}")
            print(f"  IUPAC: {bog_data['name_iupac']}")
            print(f"  First atom: {bog_data['atoms'][0]['name']} ({bog_data['atoms'][0]['type']})")
            print(f"  Charge of first atom: {bog_data['atoms'][0]['charge']}")
            print(f"  First bond: {bog_data['bonds'][0][0]} - {bog_data['bonds'][0][1]}")
        
        # Count statistics
        total_atoms = sum(len(res['atoms']) for res in residues_dict.values())
        total_bonds = sum(len(res['bonds']) for res in residues_dict.values())
        total_impropers = sum(len(res['impropers']) for res in residues_dict.values())
        
        print(f"\nOverall statistics:")
        print(f"  Total residues: {len(residues_dict)}")
        print(f"  Total atoms: {total_atoms}")
        print(f"  Total bonds: {total_bonds}")
        print(f"  Total impropers: {total_impropers}")
        
    except FileNotFoundError:
        print(f"Error: File not found at {carb_rtp_path}")
    except Exception as e:
        print(f"Error parsing file: {e}")

if __name__ == "__main__":
    main()
