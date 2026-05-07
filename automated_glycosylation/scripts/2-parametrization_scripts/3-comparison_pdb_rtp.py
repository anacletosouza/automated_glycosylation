import pickle
import math
from collections import defaultdict, Counter
import re

class CarbohydrateCompositionMatcher:
    def __init__(self, pdb_analysis_path, rtp_dict_path):
        """
        Initialize with PDB analysis and RTP dictionary
        """
        with open(pdb_analysis_path, 'rb') as f:
            self.pdb_data = pickle.load(f)
        
        with open(rtp_dict_path, 'rb') as f:
            self.rtp_data = pickle.load(f)
        
        # Reference sequence
        self.sequence = "GlcNAc(b1-2)Man(a1-6)[GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
        
        # IUPAC element mapping - standard elements only
        self.iupac_elements = {
            'C', 'N', 'O', 'H', 'S', 'P', 'F', 'Cl', 'Br', 'I',
        }
        
        # Map common atom prefixes to IUPAC elements
        self.element_prefix_map = {
            'C': 'C', 'CA': 'C', 'CB': 'C', 'CG': 'C', 'CD': 'C', 'CE': 'C', 'CZ': 'C',
            'N': 'N', 'NA': 'N', 'NB': 'N', 'NT': 'N', 'NE': 'N', #'N2' : 'N',
            'O': 'O', 'OA': 'O', 'OB': 'O', 'OG': 'O', 'OD': 'O', 'OE': 'O', 'OH': 'O',
            'H': 'H', 'HA': 'H', 'HB': 'H', 'HG': 'H', 'HD': 'H', 'HE': 'H', 'HZ': 'H',
            'S': 'S', 'SD': 'S', 'SG': 'S',
            'P': 'P',
        }
    
    def extract_element_from_atom_name(self, atom_name):
        """
        Extract IUPAC element from atom name
        Follows IUPAC naming conventions
        """
        # Remove numbers and special characters
        base_name = ''.join([c for c in atom_name if not c.isdigit()])
        
        # Check common prefixes
        for prefix, element in self.element_prefix_map.items():
            if base_name.startswith(prefix):
                return element
        
        # If starts with known element letter
        if base_name and base_name[0] in self.iupac_elements:
            return base_name[0]
        
        # Default to first character if it's a capital letter
        if base_name and base_name[0].isupper():
            return base_name[0]
        
        return 'Unknown'
    
    def get_iupac_composition(self, atoms):
        """
        Get IUPAC element composition from atoms
        atoms: list of atom dictionaries with 'name' and optionally 'element' fields
        """
        composition = Counter()
        
        for atom in atoms:
            if isinstance(atom, dict):
                # Try to get element from 'element' field first
                element = atom.get('element', '')
                if element and element in self.iupac_elements:
                    composition[element] += 1
                else:
                    # Infer from atom name
                    element = self.extract_element_from_atom_name(atom.get('name', ''))
                    if element in self.iupac_elements:
                        composition[element] += 1
        
        return dict(composition)
    
    def analyze_pdb_residues(self):
        """
        Analyze all PDB residues for IUPAC element composition
        """
        pdb_residues_info = []
        
        for residue in self.pdb_data['residues']:
            # Get IUPAC composition
            composition = self.get_iupac_composition(residue['atoms'])
            
            # Calculate total atoms for percentage calculation
            total_atoms = sum(composition.values())
            
            residue_info = {
                'pdb_code': residue['name_pdb'],
                'chain': residue['chain_carb'],
                'sequence_number': residue['res_seq'],
                'iupac_composition': composition,
                'total_atoms': total_atoms,
                'atom_names': [atom['name'] for atom in residue['atoms']]
            }
            
            # Add specific counts for common elements
            for element in ['C', 'N', 'O', 'H', 'S']:
                residue_info[f'{element}_count'] = composition.get(element, 0)
            
            pdb_residues_info.append(residue_info)
        
        return pdb_residues_info
    
    def analyze_rtp_residues(self):
        """
        Analyze all RTP residues for IUPAC element composition
        """
        rtp_residues_info = []
        
        for res_name, res_data in self.rtp_data.items():
            # Get atoms from RTP residue
            atoms = []
            for atom_info in res_data.get('atoms', []):
                atoms.append({
                    'name': atom_info['name'],
                    'element': self.extract_element_from_atom_name(atom_info['type'])
                })
            
            # Get IUPAC composition
            composition = self.get_iupac_composition(atoms)
            
            residue_info = {
                'charmm_name': res_name,
                'iupac_name': res_data.get('name_iupac', ''),
                'iupac_composition': composition,
                'total_atoms': sum(composition.values()),
                'atom_names': [atom['name'] for atom in atoms]
            }
            
            # Add specific counts
            for element in ['C', 'N', 'O', 'H', 'S']:
                residue_info[f'{element}_count'] = composition.get(element, 0)
            
            rtp_residues_info.append(residue_info)
        
        return rtp_residues_info
    
    def calculate_composition_match_percentage(self, pdb_comp, rtp_comp):
        """
        Calculate match percentage based on IUPAC element composition
        Returns: percentage (0-100) of matching element counts
        """
        # Get all elements present in either composition
        all_elements = set(list(pdb_comp.keys()) + list(rtp_comp.keys()))
        
        if not all_elements:
            return 0.0
        
        total_matching_atoms = 0
        total_atoms_in_pdb = sum(pdb_comp.values())
        
        # For each element, calculate match
        for element in all_elements:
            pdb_count = pdb_comp.get(element, 0)
            rtp_count = rtp_comp.get(element, 0)
            
            # Atoms that match (min of the two counts)
            matching = min(pdb_count, rtp_count)
            total_matching_atoms += matching
        
        # Avoid division by zero
        if total_atoms_in_pdb == 0:
            return 0.0
        
        # Percentage of PDB atoms that have a match in RTP
        match_percentage = (total_matching_atoms / total_atoms_in_pdb) * 100
        
        return match_percentage
    
    def find_all_possible_matches(self, pdb_residues_info, rtp_residues_info):
        """
        Find all possible matches for each PDB residue with match percentages
        """
        all_matches = []
        
        for pdb_res in pdb_residues_info:
            pdb_comp = pdb_res['iupac_composition']
            matches_for_residue = []
            
            for rtp_res in rtp_residues_info:
                rtp_comp = rtp_res['iupac_composition']
                
                # Calculate match percentage
                match_percentage = self.calculate_composition_match_percentage(pdb_comp, rtp_comp)
                
                # Calculate complementary metrics
                exact_element_match = all(pdb_comp.get(e, 0) == rtp_comp.get(e, 0) 
                                        for e in set(list(pdb_comp.keys()) + list(rtp_comp.keys())))
                
                # Check for common carbohydrate patterns
                has_nitrogen = (pdb_comp.get('N', 0) > 0) == (rtp_comp.get('N', 0) > 0)
                carbon_oxygen_ratio_similar = abs(pdb_comp.get('C', 0) - rtp_comp.get('C', 0)) <= 2 and \
                                            abs(pdb_comp.get('O', 0) - rtp_comp.get('O', 0)) <= 2
                
                match_info = {
                    'pdb_residue': pdb_res['pdb_code'],
                    'pdb_chain': pdb_res['chain'],
                    'pdb_seq': pdb_res['sequence_number'],
                    'rtp_charmm_name': rtp_res['charmm_name'],
                    'rtp_iupac_name': rtp_res['iupac_name'],
                    'match_percentage': match_percentage,
                    'exact_element_match': exact_element_match,
                    'has_nitrogen_match': has_nitrogen,
                    'carbon_oxygen_similar': carbon_oxygen_ratio_similar,
                    'pdb_composition': pdb_comp,
                    'rtp_composition': rtp_comp,
                    'pdb_total_atoms': pdb_res['total_atoms'],
                    'rtp_total_atoms': rtp_res['total_atoms'],
                }
                
                matches_for_residue.append(match_info)
            
            # Sort by match percentage
            matches_for_residue.sort(key=lambda x: x['match_percentage'], reverse=True)
            all_matches.append({
                'pdb_residue': pdb_res,
                'possible_matches': matches_for_residue
            })
        
        return all_matches
    
    def filter_by_expected_residues(self, all_matches):
        """
        Filter matches based on expected residues from sequence
        """
        # Expected residues from the sequence
        expected_residues = ['GlcNAc', 'Man', 'Fuc']
        
        filtered_matches = []
        
        for residue_matches in all_matches:
            pdb_res = residue_matches['pdb_residue']
            filtered = []
            
            for match in residue_matches['possible_matches']:
                rtp_iupac_name = match['rtp_iupac_name'].lower()
                
                # Check if RTP name contains any expected residue
                matches_expected = False
                matched_residue = None
                
                for expected in expected_residues:
                    if expected.lower() in rtp_iupac_name:
                        matches_expected = True
                        matched_residue = expected
                        break
                
                match['matches_expected_residue'] = matches_expected
                match['expected_residue_type'] = matched_residue
                
                # Keep all matches, but mark expected ones
                filtered.append(match)
            
            # Resort with expected matches first
            filtered.sort(key=lambda x: (
                -x['matches_expected_residue'],  # Expected first
                -x['match_percentage']           # Then by percentage
            ))
            
            filtered_matches.append({
                'pdb_residue': pdb_res,
                'possible_matches': filtered
            })
        
        return filtered_matches
    
    def display_all_matches_for_residue(self, residue_data):
        """
        Display ALL possible matches for a single PDB residue
        """
        pdb_res = residue_data['pdb_residue']
        matches = residue_data['possible_matches']
        
        print(f"\n{'='*120}")
        print(f"PDB RESIDUE: {pdb_res['pdb_code']} (Chain {pdb_res['chain']}, "
              f"Position {pdb_res['sequence_number']})")
        print(f"{'='*120}")
        
        # Display PDB residue composition
        pdb_comp = pdb_res['iupac_composition']
        comp_str = ", ".join([f"{k}{v}" for k, v in sorted(pdb_comp.items())])
        print(f"Composition: {pdb_res['total_atoms']} atoms [{comp_str}]")
        
        # Display ALL matches
        print(f"\nALL POSSIBLE RTP MATCHES ({len(matches)} total):")
        print("-" * 120)
        print(f"{'RTP Name':<12} {'IUPAC Name':<55} {'Match %':<10} {'Atoms':<10} {'Expected':<10}")
        print("-" * 120)
        
        match_count = 0
        for match in matches:
            match_count += 1
            expected_mark = "✓" if match['matches_expected_residue'] else ""
            expected_type = match['expected_residue_type'] or ""
            
            # Truncate long IUPAC names if necessary
            iupac_name = match['rtp_iupac_name']
            if len(iupac_name) > 52:
                iupac_name = iupac_name[:49] + "..."
            
            print(f"{match['rtp_charmm_name']:<12} {iupac_name:<55} "
                  f"{match['match_percentage']:>8.1f}%  "
                  f"{match['rtp_total_atoms']:>8}  "
                  f"{expected_mark:>3} {expected_type:<7}")
        
        # Show statistics for this residue
        exact_matches = [m for m in matches if m['exact_element_match']]
        high_matches = [m for m in matches if m['match_percentage'] >= 90]
        expected_matches = [m for m in matches if m['matches_expected_residue']]
        
        print(f"\nMatch Statistics for {pdb_res['pdb_code']}:")
        print(f"  • Total possible matches: {len(matches)}")
        print(f"  • Exact element matches: {len(exact_matches)}")
        print(f"  • Matches with ≥90% similarity: {len(high_matches)}")
        print(f"  • Matches with expected residue type: {len(expected_matches)}")
        
        if exact_matches:
            print(f"\nExact element matches:")
            for match in exact_matches:
                print(f"  {match['rtp_charmm_name']}: {match['rtp_iupac_name'][:60]}...")
        
        if expected_matches:
            print(f"\nMatches with expected residue types:")
            for match in expected_matches[:10]:  # Show top 10 expected
                print(f"  {match['rtp_charmm_name']} ({match['match_percentage']:.1f}%): "
                      f"{match['rtp_iupac_name'][:50]}...")
        
        return match_count
    
    def generate_comprehensive_report(self):
        """
        Generate comprehensive report with ALL possibilities
        """
        print("=" * 120)
        print("COMPREHENSIVE CARBOHYDRATE RESIDUE COMPOSITION ANALYSIS")
        print("=" * 120)
        
        print(f"\nREFERENCE SEQUENCE:")
        print(f"{self.sequence}")
        
        # Analyze compositions
        pdb_residues_info = self.analyze_pdb_residues()
        rtp_residues_info = self.analyze_rtp_residues()
        
        print(f"\nPDB RESIDUES FOUND ({len(pdb_residues_info)}):")
        print("-" * 120)
        for res in pdb_residues_info:
            comp_str = ", ".join([f"{k}{v}" for k, v in sorted(res['iupac_composition'].items())])
            print(f"  {res['pdb_code']} (Chain {res['chain']}, Pos {res['sequence_number']}): "
                  f"{res['total_atoms']} atoms [{comp_str}]")
        
        print(f"\nRTP RESIDUES AVAILABLE ({len(rtp_residues_info)} total)")
        print("-" * 120)
        
        # Find all possible matches
        all_matches = self.find_all_possible_matches(pdb_residues_info, rtp_residues_info)
        filtered_matches = self.filter_by_expected_residues(all_matches)
        
        # Generate detailed report for EACH PDB residue
        total_possible_matches = 0
        
        for residue_data in filtered_matches:
            total_possible_matches += self.display_all_matches_for_residue(residue_data)
        
        # Summary statistics
        print(f"\n{'='*120}")
        print("GLOBAL SUMMARY STATISTICS")
        print(f"{'='*120}")
        
        total_pdb = len(pdb_residues_info)
        total_rtp = len(rtp_residues_info)
        
        # Calculate statistics
        exact_match_counts = []
        high_match_counts = []
        expected_match_counts = []
        avg_best_match = 0
        
        for residue_data in filtered_matches:
            matches = residue_data['possible_matches']
            
            if matches:
                # Best match percentage
                avg_best_match += matches[0]['match_percentage']
                
                # Counts for this residue
                exact_matches = len([m for m in matches if m['exact_element_match']])
                high_matches = len([m for m in matches if m['match_percentage'] >= 90])
                expected_matches = len([m for m in matches if m['matches_expected_residue']])
                
                exact_match_counts.append(exact_matches)
                high_match_counts.append(high_matches)
                expected_match_counts.append(expected_matches)
        
        avg_best_match = avg_best_match / total_pdb if total_pdb > 0 else 0
        avg_exact_matches = sum(exact_match_counts) / total_pdb if total_pdb > 0 else 0
        avg_high_matches = sum(high_match_counts) / total_pdb if total_pdb > 0 else 0
        
        print(f"\nAnalysis Results:")
        print(f"  Total PDB residues analyzed: {total_pdb}")
        print(f"  Total RTP residues in database: {total_rtp}")
        print(f"  Total possible matches considered: {total_possible_matches}")
        print(f"  Average best match percentage: {avg_best_match:.1f}%")
        print(f"  Average exact element matches per residue: {avg_exact_matches:.1f}")
        print(f"  Average high-similarity matches (≥90%) per residue: {avg_high_matches:.1f}")
        
        # Composition analysis
        print(f"\nComposition Analysis:")
        unique_compositions = set()
        for res in pdb_residues_info:
            comp_tuple = tuple(sorted(res['iupac_composition'].items()))
            unique_compositions.add(comp_tuple)
        
        print(f"  Unique element compositions in PDB: {len(unique_compositions)}")
        
        # Show each unique composition
        print(f"\nUnique PDB Compositions Found:")
        for i, comp_tuple in enumerate(sorted(unique_compositions), 1):
            comp_str = ", ".join([f"{k}{v}" for k, v in comp_tuple])
            # Find residues with this composition
            residues_with_comp = []
            for res in pdb_residues_info:
                if tuple(sorted(res['iupac_composition'].items())) == comp_tuple:
                    residues_with_comp.append(f"{res['pdb_code']}({res['chain']}{res['sequence_number']})")
            
            print(f"  {i:2d}. [{comp_str:30}] → {', '.join(residues_with_comp)}")
        
        # Expected residue analysis
        print(f"\nExpected Residue Analysis (from sequence):")
        expected_residues = ['GlcNAc', 'Man', 'Fuc']
        
        for expected in expected_residues:
            total_for_type = 0
            best_matches_for_type = []
            
            for residue_data in filtered_matches:
                for match in residue_data['possible_matches']:
                    if match['expected_residue_type'] == expected:
                        total_for_type += 1
                        if match['match_percentage'] >= 80:  # Good matches
                            best_matches_for_type.append(
                                f"{match['pdb_residue']}({match['pdb_chain']}{match['pdb_seq']})"
                                f"→{match['rtp_charmm_name']}({match['match_percentage']:.1f}%)"
                            )
            
            print(f"  {expected}: {total_for_type} possible matches in RTP")
            if best_matches_for_type:
                print(f"    Best matches: {', '.join(best_matches_for_type[:3])}")
                if len(best_matches_for_type) > 3:
                    print(f"    ... and {len(best_matches_for_type) - 3} more")
        
        return filtered_matches
    
    def export_all_matches_to_file(self, filtered_matches, filename="all_possible_matches.txt"):
        """
        Export ALL matches to a text file for reference
        """
        with open(filename, 'w') as f:
            f.write("=" * 120 + "\n")
            f.write("ALL POSSIBLE CARBOHYDRATE MATCHES - COMPREHENSIVE LIST\n")
            f.write("=" * 120 + "\n\n")
            
            f.write(f"Reference sequence: {self.sequence}\n\n")
            
            for residue_data in filtered_matches:
                pdb_res = residue_data['pdb_residue']
                matches = residue_data['possible_matches']
                
                f.write(f"\n{'='*120}\n")
                f.write(f"PDB RESIDUE: {pdb_res['pdb_code']} (Chain {pdb_res['chain']}, "
                       f"Position {pdb_res['sequence_number']})\n")
                f.write(f"{'='*120}\n")
                
                # Composition
                pdb_comp = pdb_res['iupac_composition']
                comp_str = ", ".join([f"{k}{v}" for k, v in sorted(pdb_comp.items())])
                f.write(f"Composition: {pdb_res['total_atoms']} atoms [{comp_str}]\n\n")
                
                # Header
                f.write(f"{'RTP Name':<12} {'IUPAC Name':<55} {'Match %':<10} {'Atoms':<10} {'Expected':<10}\n")
                f.write("-" * 120 + "\n")
                
                # All matches
                for match in matches:
                    expected_mark = "✓" if match['matches_expected_residue'] else ""
                    expected_type = match['expected_residue_type'] or ""
                    
                    iupac_name = match['rtp_iupac_name']
                    if len(iupac_name) > 52:
                        iupac_name = iupac_name[:49] + "..."
                    
                    f.write(f"{match['rtp_charmm_name']:<12} {iupac_name:<55} "
                           f"{match['match_percentage']:>8.1f}%  "
                           f"{match['rtp_total_atoms']:>8}  "
                           f"{expected_mark:>3} {expected_type:<7}\n")
                
                # Statistics for this residue
                exact_matches = [m for m in matches if m['exact_element_match']]
                high_matches = [m for m in matches if m['match_percentage'] >= 90]
                expected_matches = [m for m in matches if m['matches_expected_residue']]
                
                f.write(f"\nStatistics:\n")
                f.write(f"  • Total possible matches: {len(matches)}\n")
                f.write(f"  • Exact element matches: {len(exact_matches)}\n")
                f.write(f"  • Matches with ≥90% similarity: {len(high_matches)}\n")
                f.write(f"  • Matches with expected residue type: {len(expected_matches)}\n")
                
                if exact_matches:
                    f.write(f"\nExact element matches:\n")
                    for em in exact_matches[:5]:  # Show top 5 exact matches
                        f.write(f"  {em['rtp_charmm_name']}: {em['rtp_iupac_name']}\n")
                
                f.write("\n")
        
        print(f"\nAll possible matches exported to: {filename}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run CarbohydrateCompositionMatcher with custom pickle files.")
    parser.add_argument("--pdb", default="D_1_parser.pkl", help="Path to the PDB analysis pickle file")
    parser.add_argument("--rtp", default="carb_residues.pkl", help="Path to the RTP residues pickle file")
    args = parser.parse_args()

    # Initialize matcher with user-provided paths
    matcher = CarbohydrateCompositionMatcher(
        pdb_analysis_path=args.pdb,
        rtp_dict_path=args.rtp
    )
    
    # Generate comprehensive report showing ALL matches
    filtered_matches = matcher.generate_comprehensive_report()
    
    # Export ALL matches to text file
    matcher.export_all_matches_to_file(filtered_matches)
    
    # Condensed summary by residue type
    print(f"\n{'='*120}")
    print("CONDENSED SUMMARY BY RESIDUE TYPE")
    print(f"{'='*120}")
    
    # Group by PDB residue code
    pdb_code_groups = defaultdict(list)
    for residue_data in filtered_matches:
        pdb_res = residue_data['pdb_residue']
        pdb_code = pdb_res['pdb_code']
        pdb_code_groups[pdb_code].append(residue_data)
    
    for pdb_code, residues in sorted(pdb_code_groups.items()):
        print(f"\n{pdb_code} residues ({len(residues)} total):")
        for residue_data in residues:
            pdb_res = residue_data['pdb_residue']
            matches = residue_data['possible_matches']
            expected_matches = [m for m in matches if m['matches_expected_residue']]
            exact_matches = [m for m in matches if m['exact_element_match']]
            print(f"  Position {pdb_res['chain']}{pdb_res['sequence_number']}: "
                  f"{len(matches)} total matches, "
                  f"{len(expected_matches)} expected, "
                  f"{len(exact_matches)} exact")
            if expected_matches:
                print(f"    Top expected: ", end="")
                top_expected = expected_matches[:3]
                for match in top_expected:
                    print(f"{match['rtp_charmm_name']}({match['match_percentage']:.1f}%) ", end="")
                print()

