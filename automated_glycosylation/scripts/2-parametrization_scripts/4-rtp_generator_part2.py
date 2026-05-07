import pickle
import math
from collections import defaultdict, Counter
import re
import json
import argparse

class RTPMatcherGenerator:
    def __init__(self, pdb_analysis_path, rtp_dict_path, json_mapping=None):
        """
        Initialize with PDB analysis, RTP dictionary, and optional JSON mapping
        """
        with open(pdb_analysis_path, 'rb') as f:
            self.pdb_data = pickle.load(f)
        
        with open(rtp_dict_path, 'rb') as f:
            self.rtp_data = pickle.load(f)
        
        # Reference sequence - can be overridden by JSON
        self.sequence = "GlcNAc(b1-2)Man(a1-6)[GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
        
        # Process JSON mapping if provided
        self.json_mapping = json_mapping
        self.pdb_to_charmm = {}
        
        if json_mapping:
            print(f"Loading JSON mapping: {json_mapping}")
            with open(json_mapping, 'r') as f:
                mapping_data = json.load(f)
            
            # Override sequence if provided in JSON
            if 'sequence' in mapping_data:
                self.sequence = mapping_data['sequence']
                print(f"Sequence overridden from JSON: {self.sequence}")
            
            # Create PDB to CHARMM mapping
            if 'residues' in mapping_data:
                for residue_map in mapping_data['residues']:
                    pdb_name = residue_map.get('pdb_name')
                    charmm_name = residue_map.get('charmm_name')
                    if pdb_name and charmm_name:
                        self.pdb_to_charmm[pdb_name] = charmm_name
                
                print(f"JSON mapping loaded: {self.pdb_to_charmm}")
        
        # IUPAC element mapping
        self.iupac_elements = {
            'C', 'N', 'O', 'H', 'S', 'P', 'F', 'Cl', 'Br', 'I',
        }
        
        # Map common atom prefixes to IUPAC elements
        self.element_prefix_map = {
            'C': 'C', 'CA': 'C', 'CB': 'C', 'CG': 'C', 'CD': 'C', 'CE': 'C', 'CZ': 'C',
            'N': 'N', 'NA': 'N', 'NB': 'N', 'NT': 'N', 'NE': 'N', 'N2': 'N', 
            'O': 'O', 'OA': 'O', 'OB': 'O', 'OG': 'O', 'OD': 'O', 'OE': 'O', 'OH': 'O',
            'H': 'H', 'HA': 'H', 'HB': 'H', 'HG': 'H', 'HD': 'H', 'HE': 'H', 'HZ': 'H',
            'S': 'S', 'SD': 'S', 'SG': 'S',
            'P': 'P',
        }
        
        # Default mapping (used if no JSON or for unmapped residues)
        self.default_pdb_to_expected = {
            'NDG': 'GlcNAc',
            'NAG': 'GlcNAc',
            'MAN': 'Man',
            'BMA': 'Man',  
            'FUC': 'Fuc',
            'FCA': 'Fuc',  # Added for JSON compatibility
            'SIA': 'Neu5Ac',  
            'GLC': 'Glc',
            'GAL': 'Gal',
            'XYP': 'Xyl',
        }
    
    def get_expected_type(self, pdb_code):
        """
        Get expected type for PDB code, using JSON mapping first, then default
        """
        # First try JSON mapping
        if self.pdb_to_charmm and pdb_code in self.pdb_to_charmm:
            charmm_name = self.pdb_to_charmm[pdb_code]
            # Try to find expected type based on CHARMM name
            for default_pdb, expected in self.default_pdb_to_expected.items():
                if expected.lower() in charmm_name.lower():
                    return expected
            # If not found, return a simplified version of CHARMM name
            return charmm_name.replace('B', '').replace('A', '').replace('G', '')
        
        # If no JSON or not mapped, use default
        return self.default_pdb_to_expected.get(pdb_code, 'Unknown')
    
    def extract_element_from_atom_name(self, atom_name):
        """
        Extract IUPAC element from atom name
        """
        base_name = ''.join([c for c in atom_name if not c.isdigit()])
        
        for prefix, element in self.element_prefix_map.items():
            if base_name.startswith(prefix):
                return element
        
        if base_name and base_name[0] in self.iupac_elements:
            return base_name[0]
        
        if base_name and base_name[0].isupper():
            return base_name[0]
        
        return 'Unknown'
    
    def get_iupac_composition(self, atoms):
        """
        Get IUPAC element composition from atoms
        """
        composition = Counter()
        
        for atom in atoms:
            if isinstance(atom, dict):
                element = atom.get('element', '')
                if element and element in self.iupac_elements:
                    composition[element] += 1
                else:
                    element = self.extract_element_from_atom_name(atom.get('name', ''))
                    if element in self.iupac_elements:
                        composition[element] += 1
        
        return dict(composition)
    
    def analyze_pdb_residues(self):
        """
        Analyze all PDB residues
        """
        pdb_residues_info = []
        
        for residue in self.pdb_data['residues']:
            composition = self.get_iupac_composition(residue['atoms'])
            total_atoms = sum(composition.values())
            
            pdb_code = residue['name_pdb']
            expected_type = self.get_expected_type(pdb_code)
            
            residue_info = {
                'pdb_code': pdb_code,
                'chain': residue['chain_carb'],
                'sequence_number': residue['res_seq'],
                'iupac_composition': composition,
                'total_atoms': total_atoms,
                'atom_names': [atom['name'] for atom in residue['atoms']],
                'expected_type': expected_type,
                'original_atoms': residue['atoms']  
            }
            
            # Add specific counts
            for element in ['C', 'N', 'O', 'H', 'S']:
                residue_info[f'{element}_count'] = composition.get(element, 0)
            
            pdb_residues_info.append(residue_info)
        
        return pdb_residues_info
    
    def analyze_rtp_residues(self):
        """
        Analyze all RTP residues
        """
        rtp_residues_info = []
        
        for res_name, res_data in self.rtp_data.items():
            atoms = []
            for atom_info in res_data.get('atoms', []):
                atoms.append({
                    'name': atom_info['name'],
                    'type': atom_info['type'],
                    'charge': atom_info['charge'],
                    'charge_group': atom_info['charge_group'],
                    'element': self.extract_element_from_atom_name(atom_info['type'])
                })
            
            composition = self.get_iupac_composition(atoms)
            
            residue_info = {
                'charmm_name': res_name,
                'iupac_name': res_data.get('name_iupac', ''),
                'iupac_composition': composition,
                'total_atoms': sum(composition.values()),
                'atom_names': [atom['name'] for atom in atoms],
                'bonds': res_data.get('bonds', []),
                'impropers': res_data.get('impropers', []),
                'atoms': atoms,  
                'original_data': res_data  
            }
            
            # Add specific counts
            for element in ['C', 'N', 'O', 'H', 'S']:
                residue_info[f'{element}_count'] = composition.get(element, 0)
            
            rtp_residues_info.append(residue_info)
        
        return rtp_residues_info
    
    def calculate_composition_match_percentage(self, pdb_comp, rtp_comp):
        """
        Calculate match percentage based on IUPAC element composition
        """
        all_elements = set(list(pdb_comp.keys()) + list(rtp_comp.keys()))
        
        if not all_elements:
            return 0.0
        
        total_matching_atoms = 0
        total_atoms_in_pdb = sum(pdb_comp.values())
        
        for element in all_elements:
            pdb_count = pdb_comp.get(element, 0)
            rtp_count = rtp_comp.get(element, 0)
            matching = min(pdb_count, rtp_count)
            total_matching_atoms += matching
        
        if total_atoms_in_pdb == 0:
            return 0.0
        
        match_percentage = (total_matching_atoms / total_atoms_in_pdb) * 100
        return match_percentage
    
    def find_all_matches(self):
        """
        Find all possible matches for each PDB residue
        """
        pdb_residues_info = self.analyze_pdb_residues()
        rtp_residues_info = self.analyze_rtp_residues()
        
        all_matches = []
        
        for pdb_res in pdb_residues_info:
            pdb_comp = pdb_res['iupac_composition']
            matches_for_residue = []
            
            for rtp_res in rtp_residues_info:
                rtp_comp = rtp_res['iupac_composition']
                match_percentage = self.calculate_composition_match_percentage(pdb_comp, rtp_comp)
                
                # Check if it matches expected type
                expected_type = pdb_res['expected_type'].lower()
                rtp_iupac_name = rtp_res['iupac_name'].lower()
                matches_expected = expected_type in rtp_iupac_name if expected_type != 'unknown' else False
                
                # Check direct JSON mapping
                json_mapped = False
                if self.pdb_to_charmm and pdb_res['pdb_code'] in self.pdb_to_charmm:
                    json_mapped = (rtp_res['charmm_name'] == self.pdb_to_charmm[pdb_res['pdb_code']])
                
                match_info = {
                    'pdb_residue': pdb_res,
                    'rtp_residue': rtp_res,
                    'match_percentage': match_percentage,
                    'matches_expected': matches_expected,
                    'json_mapped': json_mapped,
                    'pdb_composition': pdb_comp,
                    'rtp_composition': rtp_comp,
                }
                
                matches_for_residue.append(match_info)
            
            # Sort matches: first JSON, then expected, then percentage
            matches_for_residue.sort(key=lambda x: (
                -x['json_mapped'],        # JSON mapping first
                -x['matches_expected'],   # Expected then
                -x['match_percentage']    # Finally by percentage
            ))
            
            all_matches.append({
                'pdb_residue': pdb_res,
                'possible_matches': matches_for_residue
            })
        
        return all_matches
    
    def auto_select_matches(self, all_matches):
        """
        Automatically select matches based on JSON mapping or best matches
        """
        selected_matches = []
        
        for residue_matches in all_matches:
            pdb_res = residue_matches['pdb_residue']
            matches = residue_matches['possible_matches']
            
            # Check if there's JSON mapping for this residue
            pdb_code = pdb_res['pdb_code']
            if self.pdb_to_charmm and pdb_code in self.pdb_to_charmm:
                target_charmm = self.pdb_to_charmm[pdb_code]
                # Look for exact match by CHARMM name
                for match in matches:
                    if match['rtp_residue']['charmm_name'] == target_charmm:
                        selected_matches.append(match)
                        print(f"✓ Auto-selection via JSON: {pdb_code} → {target_charmm}")
                        break
                else:
                    # If not found, use first match
                    if matches:
                        selected_matches.append(matches[0])
                        print(f"⚠️ JSON mapping not found, using best match: {pdb_code} → {matches[0]['rtp_residue']['charmm_name']}")
            else:
                # Without JSON, use best match
                if matches:
                    selected_matches.append(matches[0])
                    print(f"✓ Auto-selection: {pdb_code} → {matches[0]['rtp_residue']['charmm_name']} ({matches[0]['match_percentage']:.1f}%)")
        
        return selected_matches
    
    def display_match_selection_menu(self, all_matches):
        """
        Display interactive menu for user to select matches
        """
        print("=" * 100)
        print("CARBOHYDRATE RTP MATCH SELECTION")
        print("=" * 100)
        
        selected_matches = []
        
        for i, residue_matches in enumerate(all_matches):
            pdb_res = residue_matches['pdb_residue']
            matches = residue_matches['possible_matches']
            
            print(f"\n{'='*80}")
            print(f"RESIDUE {i+1}: {pdb_res['pdb_code']} (Chain {pdb_res['chain']}, Position {pdb_res['sequence_number']})")
            print(f"Expected type: {pdb_res['expected_type']}")
            
            # Show JSON mapping if exists
            if self.pdb_to_charmm and pdb_res['pdb_code'] in self.pdb_to_charmm:
                print(f"JSON mapping: → {self.pdb_to_charmm[pdb_res['pdb_code']]}")
            
            # IUPAC composition
            pdb_comp = pdb_res['iupac_composition']
            comp_str = ", ".join([f"{k}{v}" for k, v in sorted(pdb_comp.items())])
            print(f"Composition: {pdb_res['total_atoms']} atoms [{comp_str}]")
            
            # Available options
            print(f"\nAvailable options:")
            print("-" * 80)
            print(f"{'#':<3} {'RTP Name':<12} {'IUPAC Name':<40} {'Match %':<10} {'Atoms':<8} {'JSON':<6} {'Exp':<4}")
            print("-" * 80)
            
            for j, match in enumerate(matches[:10]):  
                rtp_res = match['rtp_residue']
                json_mark = "✓" if match['json_mapped'] else ""
                expected_mark = "✓" if match['matches_expected'] else ""
                
                # Truncate long IUPAC name
                iupac_name = rtp_res['iupac_name']
                if len(iupac_name) > 38:
                    iupac_name = iupac_name[:35] + "..."
                
                print(f"{j+1:<3} {rtp_res['charmm_name']:<12} {iupac_name:<40} "
                      f"{match['match_percentage']:>8.1f}%  "
                      f"{rtp_res['total_atoms']:>7}  "
                      f"{json_mark:>5}   {expected_mark:>3}")
            
            # Manual option
            print(f"{'0':<3} {'MANUAL':<12} {'Enter CHARMM name manually':<40} {'':<10} {'':<8} {'':<6} {'':<4}")
            
            # Selection
            while True:
                try:
                    choice = input(f"\nSelect option for {pdb_res['pdb_code']} (1-{min(10, len(matches))}, 0=manual, Enter=auto/JSON): ").strip()
                    
                    if choice == "":
                        # Auto-selection based on JSON or best match
                        if match.get('json_mapped'):
                            selected = match
                        else:
                            selected = matches[0] if matches else None
                        break
                    elif choice == "0":
                        # Manual input
                        manual_name = input("Enter CHARMM residue name: ").strip().upper()
                        
                        # Check if name exists in RTP data
                        if manual_name in self.rtp_data:
                            # Look for existing match
                            for match in matches:
                                if match['rtp_residue']['charmm_name'] == manual_name:
                                    selected = match
                                    break
                            else:
                                # Create new manual match
                                rtp_res_data = self.rtp_data[manual_name]
                                selected = {
                                    'pdb_residue': pdb_res,
                                    'rtp_residue': {
                                        'charmm_name': manual_name,
                                        'iupac_name': rtp_res_data.get('name_iupac', ''),
                                        'original_data': rtp_res_data
                                    },
                                    'match_percentage': 0.0,
                                    'matches_expected': False,
                                    'json_mapped': False,
                                    'manual_selection': True
                                }
                        else:
                            print(f"Residue {manual_name} was not found in RTP.")
                            continue
                        
                        break
                    else:
                        choice_idx = int(choice) - 1
                        if 0 <= choice_idx < len(matches[:10]):
                            selected = matches[choice_idx]
                            break
                        else:
                            print("Invalid option. Try again.")
                except ValueError:
                    print("Invalid input. Enter a number.")
            
            if selected:
                selected_matches.append(selected)
                print(f"✓ Selected: {selected['rtp_residue']['charmm_name']}")
            else:
                print(f"✗ {pdb_res['pdb_code']} skipped")
        
        return selected_matches
    
    def generate_modified_rtp(self, selected_matches, output_filename="carb_modified.rtp"):
        """
        Generate modified RTP file with selected matches
        """
        print(f"\n{'='*80}")
        print(f"GENERATING MODIFIED RTP FILE: {output_filename}")
        print(f"{'='*80}")
        
        # Atom name mapping
        def map_atom_names(pdb_atom_name, rtp_atom_names, residue_type):
            """
            Mapping atom names from PDB to CHARMM
            """
            # Common mappings for carbohydrates
            common_mappings = {
                'C1': 'C1', 'C2': 'C2', 'C3': 'C3', 'C4': 'C4', 'C5': 'C5', 'C6': 'C6',
                'O1': 'O1', 'O2': 'O2', 'O3': 'O3', 'O4': 'O4', 'O5': 'O5', 'O6': 'O6',
                'N2': 'N', 'N': 'N',
                'H1': 'H1', 'H2': 'H2', 'H3': 'H3', 'H4': 'H4', 'H5': 'H5',
                'H61': 'H61', 'H62': 'H62', 'H11': 'H11', 'H12': 'H12',
            }
            
            # Try common mapping first
            if pdb_atom_name in common_mappings:
                mapped_name = common_mappings[pdb_atom_name]
                if mapped_name in rtp_atom_names:
                    return mapped_name
            
            # For hydrogens
            if pdb_atom_name.startswith('H'):
                base_name = pdb_atom_name[1:]  
                for rtp_name in rtp_atom_names:
                    if rtp_name.startswith('H') and base_name in rtp_name:
                        return rtp_name
                    elif 'H' + base_name == rtp_name:
                        return rtp_name
            
            # For carbons
            if pdb_atom_name.startswith('C'):
                for rtp_name in rtp_atom_names:
                    if rtp_name.startswith('C') and pdb_atom_name[1:] in rtp_name:
                        return rtp_name
            
            # For oxygens
            if pdb_atom_name.startswith('O'):
                for rtp_name in rtp_atom_names:
                    if rtp_name.startswith('O') and pdb_atom_name[1:] in rtp_name:
                        return rtp_name
            
            return pdb_atom_name
        
        # Generate RTP content
        rtp_content = []
        rtp_content.append("; =============================================================================")
        rtp_content.append("; MODIFIED CARBOHYDRATE RTP FILE")
        rtp_content.append("; Generated by RTPMatcherGenerator")
        rtp_content.append("; Based on PDB analysis and CHARMM carbohydrate parameters")
        if self.json_mapping:
            rtp_content.append(f"; JSON mapping file: {self.json_mapping}")
        rtp_content.append("; =============================================================================")
        rtp_content.append("")
        
        # Add selected residues
        for i, match in enumerate(selected_matches):
            pdb_res = match['pdb_residue']
            rtp_res = match['rtp_residue']
            rtp_data = rtp_res.get('original_data', {})
            
            if not rtp_data:
                print(f"Warning: RTP data not found for {rtp_res['charmm_name']}")
                continue
            
            residue_name = f"{pdb_res['pdb_code']}_{pdb_res['chain']}{pdb_res['sequence_number']}"
            
            rtp_content.append(f"[ {residue_name} ]")
            rtp_content.append(f"; Original CHARMM name: {rtp_res['charmm_name']}")
            rtp_content.append(f"; IUPAC name: {rtp_res.get('iupac_name', 'Unknown')}")
            rtp_content.append(f"; PDB code: {pdb_res['pdb_code']}, Chain: {pdb_res['chain']}, Position: {pdb_res['sequence_number']}")
            rtp_content.append(f"; Match percentage: {match.get('match_percentage', 0):.1f}%")
            
            # Atoms
            rtp_content.append(" [ atoms ]")
            
            # List of CHARMM atom names
            charmm_atom_names = [atom['name'] for atom in rtp_data.get('atoms', [])]
            
            # Add atoms
            for atom in rtp_data.get('atoms', []):
                name = atom['name']
                atom_type = atom['type']
                charge = atom['charge']
                charge_group = atom['charge_group']
                
                rtp_content.append(f"  {name:6} {atom_type:8} {charge:8.4f} {charge_group}")
            
            # Bonds
            if rtp_data.get('bonds'):
                rtp_content.append(" [ bonds ]")
                for bond in rtp_data['bonds']:
                    rtp_content.append(f"  {bond[0]:6} {bond[1]:6}")
            
            # Impropers
            if rtp_data.get('impropers'):
                rtp_content.append(" [ impropers ]")
                for improper in rtp_data['impropers']:
                    rtp_content.append(f"  {improper[0]:6} {improper[1]:6} {improper[2]:6} {improper[3]:6}")
            
            # Exceptions
            if rtp_data.get('exceptions'):
                rtp_content.append(" [ exclusions ]")
                for excl in rtp_data['exceptions']:
                    rtp_content.append(f"  {excl[0]:6} {excl[1]:6}")
            
            rtp_content.append("")
        
        # Add original reference residues
        rtp_content.append("; =============================================================================")
        rtp_content.append("; ORIGINAL CHARMM CARBOHYDRATE RESIDUES FOR REFERENCE")
        rtp_content.append("; =============================================================================")
        rtp_content.append("")
        
        # Common residues
        common_residues = ['NAG', 'MAN', 'FUC', 'GLC', 'GAL', 'BMA']
        for res_name in common_residues:
            # Find corresponding entry in RTP data
            for rtp_name, rtp_data in self.rtp_data.items():
                if res_name in rtp_name or res_name in rtp_data.get('name_iupac', ''):
                    rtp_content.append(f"[ {rtp_name} ]")
                    rtp_content.append(f"; IUPAC name: {rtp_data.get('name_iupac', 'Unknown')}")
                    
                    # Atoms
                    rtp_content.append(" [ atoms ]")
                    for atom in rtp_data.get('atoms', []):
                        rtp_content.append(f"  {atom['name']:6} {atom['type']:8} {atom['charge']:8.4f} {atom['charge_group']}")
                    
                    # Bonds
                    if rtp_data.get('bonds'):
                        rtp_content.append(" [ bonds ]")
                        for bond in rtp_data['bonds']:
                            rtp_content.append(f"  {bond[0]:6} {bond[1]:6}")
                    
                    rtp_content.append("")
                    break
        
        # Write file
        with open(output_filename, 'w') as f:
            f.write("\n".join(rtp_content))
        
        print(f"\n✓ File {output_filename} generated successfully!")
        print(f"  Total residues included: {len(selected_matches)}")
        print(f"  Original reference residues: {len(common_residues)}")
        
        # Also generate a mapping file
        self.generate_mapping_file(selected_matches, "carb_mapping.txt")
        
        return output_filename
    
    def generate_mapping_file(self, selected_matches, filename="carb_mapping.txt"):
        """
        Generate a mapping file showing PDB to CHARMM residue mapping
        """
        with open(filename, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("PDB TO CHARMM CARBOHYDRATE RESIDUE MAPPING\n")
            f.write("=" * 80 + "\n\n")
            
            # Information about data source
            if self.json_mapping:
                f.write(f"JSON mapping source: {self.json_mapping}\n")
                f.write(f"Sequence: {self.sequence}\n\n")
            
            f.write(f"{'PDB Residue':<15} {'Chain':<8} {'Position':<10} {'CHARMM Name':<15} {'IUPAC Name':<40} {'Match %':<10}\n")
            f.write("-" * 100 + "\n")
            
            for match in selected_matches:
                pdb_res = match['pdb_residue']
                rtp_res = match['rtp_residue']
                
                # Truncate long IUPAC name
                iupac_name = rtp_res.get('iupac_name', '')
                if len(iupac_name) > 38:
                    iupac_name = iupac_name[:35] + "..."
                
                f.write(f"{pdb_res['pdb_code']:<15} {pdb_res['chain']:<8} {pdb_res['sequence_number']:<10} "
                       f"{rtp_res['charmm_name']:<15} {iupac_name:<40} {match.get('match_percentage', 0):>8.1f}%\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("NOTES:\n")
            f.write("-" * 80 + "\n")
            f.write("1. Residue names in the RTP file follow the pattern: PDBcode_ChainPosition\n")
            f.write("2. Original CHARMM names are preserved in comments\n")
            f.write("3. Atom names and parameters are from CHARMM carbohydrate force field\n")
            if self.json_mapping:
                f.write("4. JSON mapping was used for residue selection\n")
            f.write("5. Use this mapping when setting up your simulation system\n")
        
        print(f"✓ Mapping file generated: {filename}")
    
    def run_interactive_selection(self):
        """
        Run complete interactive selection and RTP generation
        """
        print("=" * 100)
        print("INTERACTIVE CARBOHYDRATE RTP GENERATOR")
        print("=" * 100)
        
        # Inform about operation mode
        if self.json_mapping:
            print(f"Mode: Using JSON mapping ({self.json_mapping})")
            print(f"Sequence: {self.sequence}")
        else:
            print("Mode: Interactive (no JSON)")
        
        # Find all matches
        print("\nAnalyzing compositions and finding matches...")
        all_matches = self.find_all_matches()
        
        print(f"✓ Analysis complete!")
        print(f"  PDB residues found: {len(all_matches)}")
        
        # Show quick statistics
        total_matches = sum(len(res['possible_matches']) for res in all_matches)
        print(f"  Total possible matches: {total_matches}")
        
        # Select matches
        if self.json_mapping:
            print("\nUsing JSON mapping for automatic selection...")
            selected_matches = self.auto_select_matches(all_matches)
        else:
            # Interactive menu
            selected_matches = self.display_match_selection_menu(all_matches)
        
        if not selected_matches:
            print("\nNo residues selected. Operation cancelled.")
            return
        
        # Generate RTP file
        output_file = "carb_modified.rtp"
        
        self.generate_modified_rtp(selected_matches, output_file)
        
        # Show final summary
        print(f"\n{'='*80}")
        print("SELECTION SUMMARY")
        print(f"{'='*80}")
        
        for i, match in enumerate(selected_matches):
            pdb_res = match['pdb_residue']
            rtp_res = match['rtp_residue']
            
            json_info = " (JSON)" if match.get('json_mapped') else ""
            print(f"{i+1:2d}. {pdb_res['pdb_code']}({pdb_res['chain']}{pdb_res['sequence_number']}) → "
                  f"{rtp_res['charmm_name']} ({match.get('match_percentage', 0):.1f}%){json_info}")
        
        print(f"\nProcess completed successfully!")


# Main execution
def main():
    parser = argparse.ArgumentParser(description='Carbohydrate RTP Matcher Generator')
    parser.add_argument('--pdb', default='D_1_parser.pkl', help='Path to PDB analysis pickle file')
    parser.add_argument('--rtp', default='carb_residues.pkl', help='Path to RTP dictionary pickle file')
    parser.add_argument('--json', help='Path to JSON mapping file (optional)')
    parser.add_argument('--output', default='carb_modified.rtp', help='Output RTP file name')
    
    args = parser.parse_args()
    
    # Initialize generator
    generator = RTPMatcherGenerator(
        pdb_analysis_path=args.pdb,
        rtp_dict_path=args.rtp,
        json_mapping=args.json
    )
    
    # Run interactive selection
    generator.run_interactive_selection()

if __name__ == "__main__":
    main()
