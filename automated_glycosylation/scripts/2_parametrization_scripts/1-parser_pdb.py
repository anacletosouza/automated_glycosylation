import numpy as np
import itertools
from collections import defaultdict
import math
import pickle

class CarbohydratePDBParser:
    def __init__(self, pdb_content):
        """
        Initializes the parser with PDB content
        """
        self.pdb_content = pdb_content
        self.residues = []
        self.atoms = []
        self.parse_pdb()
    
    def parse_pdb(self):
        """
        Parses the PDB file
        """
        lines = self.pdb_content.strip().split('\n')
        
        for line in lines:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                atom_data = self.parse_atom_line(line)
                if atom_data:
                    self.atoms.append(atom_data)
        
        # Groups atoms by residue
        self.group_atoms_by_residue()
    
    def parse_atom_line(self, line):
        """
        Parses a PDB line in HETATM/ATOM format
        """
        try:
            # Parse according to standard PDB format
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
                'x': x,
                'y': y,
                'z': z,
                'element': line[76:78].strip() if len(line) > 76 else ''
            }
        except (ValueError, IndexError) as e:
            print(f"Error parsing line: {line}")
            return None
    
    def group_atoms_by_residue(self):
        """
        Groups atoms by residue
        """
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
        
        # Converts to list ordered by res_seq
        self.residues = sorted(residues_dict.values(), key=lambda x: x['res_seq'])
    
    def calculate_distance(self, atom1, atom2):
        """
        Calculates Euclidean distance between two atoms
        """
        dx = atom1['x'] - atom2['x']
        dy = atom1['y'] - atom2['y']
        dz = atom1['z'] - atom2['z']
        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    def generate_distance_matrix(self):
        """
        Generates distance matrix between all atoms
        """
        n_atoms = len(self.atoms)
        distance_matrix = np.zeros((n_atoms, n_atoms))
        
        for i in range(n_atoms):
            for j in range(n_atoms):
                if i != j:
                    distance_matrix[i][j] = self.calculate_distance(self.atoms[i], self.atoms[j])
        
        return distance_matrix
    
    def identify_bonds(self, max_distance=1.8):
        """
        Identifies bonds between atoms based on distance
        """
        bonds = []
        
        for i, atom1 in enumerate(self.atoms):
            for j, atom2 in enumerate(self.atoms):
                if i < j:  # Avoids duplicates
                    distance = self.calculate_distance(atom1, atom2)
                    
                    # Criteria for covalent bonds (adjustable)
                    if distance < max_distance:
                        bonds.append({
                            'atom1': atom1['name'],
                            'atom1_serial': atom1['serial'],
                            'atom1_res': atom1['res_name'],
                            'atom2': atom2['name'],
                            'atom2_serial': atom2['serial'],
                            'atom2_res': atom2['res_name'],
                            'distance': distance
                        })
        
        return bonds
    
    def get_intra_residue_bonds(self, max_distance=1.8):
        """
        Returns bonds within the same residue
        """
        bonds = self.identify_bonds(max_distance)
        intra_bonds = []
        
        for bond in bonds:
            if (bond['atom1_res'] == bond['atom2_res'] and
                self.get_residue_for_atom(bond['atom1_serial'])['res_seq'] == 
                self.get_residue_for_atom(bond['atom2_serial'])['res_seq']):
                intra_bonds.append(bond)
        
        return intra_bonds
    
    def get_inter_residue_bonds(self, max_distance=1.8):
        """
        Returns bonds between different residues
        """
        bonds = self.identify_bonds(max_distance)
        inter_bonds = []
        
        for bond in bonds:
            atom1_res = self.get_residue_for_atom(bond['atom1_serial'])
            atom2_res = self.get_residue_for_atom(bond['atom2_serial'])
            
            if (atom1_res['name_pdb'] != atom2_res['name_pdb'] or
                atom1_res['res_seq'] != atom2_res['res_seq']):
                inter_bonds.append(bond)
        
        return inter_bonds
    
    def get_residue_for_atom(self, atom_serial):
        """
        Finds the residue for a specific atom
        """
        for residue in self.residues:
            for atom in residue['atoms']:
                if atom['serial'] == atom_serial:
                    return residue
        return None
    
    def get_atom_coordinates(self):
        """
        Returns coordinates of all atoms
        """
        coordinates = []
        for atom in self.atoms:
            coordinates.append({
                'serial': atom['serial'],
                'name': atom['name'],
                'residue': atom['res_name'],
                'x': atom['x'],
                'y': atom['y'],
                'z': atom['z']
            })
        return coordinates
    
    def generate_analysis(self, max_bond_distance=1.8):
        """
        Generates complete carbohydrate analysis
        """
        analysis = {
            'residues': [],
            'atoms_pdb': self.get_atom_coordinates(),
            'coordinates': [(atom['x'], atom['y'], atom['z']) for atom in self.atoms],
            'distance_matrix_bounds': self.generate_distance_matrix().tolist(),
            'bonds_pdb': self.identify_bonds(max_bond_distance),
            'intra_residue_bounds': self.get_intra_residue_bonds(max_bond_distance),
            'inter_residue_bounds': self.get_inter_residue_bonds(max_bond_distance)
        }
        
        # Adds residue information
        for residue in self.residues:
            residue_info = {
                'name_pdb': residue['name_pdb'],
                'chain_carb': residue['chain_carb'],
                'res_seq': residue['res_seq'],
                'atoms': [{
                    'name': atom['name'],
                    'serial': atom['serial'],
                    'x': atom['x'],
                    'y': atom['y'],
                    'z': atom['z'],
                    'element': atom['element']
                } for atom in residue['atoms']]
            }
            analysis['residues'].append(residue_info)
        
        return analysis
    
    def print_summary(self):
        """
        Prints an analysis summary
        """
        print("=" * 60)
        print("PDB CARBOHYDRATE ANALYSIS")
        print("=" * 60)
        
        print(f"\nTotal residues: {len(self.residues)}")
        print(f"Total atoms: {len(self.atoms)}")
        
        print("\nResidues found:")
        for residue in self.residues:
            print(f"  {residue['name_pdb']} ({residue['chain_carb']}{residue['res_seq']}): "
                  f"{len(residue['atoms'])} atoms")
        
        analysis = self.generate_analysis()
        
        print(f"\nTotal bonds identified: {len(analysis['bonds_pdb'])}")
        print(f"Intra-residue bonds: {len(analysis['intra_residue_bounds'])}")
        print(f"Inter-residue bonds: {len(analysis['inter_residue_bounds'])}")
        
        print("\nDistance matrix:")
        print(f"  Dimensions: {len(analysis['distance_matrix_bounds'])} x "
              f"{len(analysis['distance_matrix_bounds'][0])}")
        
        # Example of some inter-residue bonds
        if analysis['inter_residue_bounds']:
            print("\nExample of inter-residue bonds (first 5):")
            for bond in analysis['inter_residue_bounds'][:5]:
                print(f"  {bond['atom1']} ({bond['atom1_res']}) -- "
                      f"{bond['atom2']} ({bond['atom2_res']}): "
                      f"{bond['distance']:.3f} Å")



if __name__ == "__main__":
    import argparse

    # Create argument parser
    parser_arg = argparse.ArgumentParser(description="Parse a PDB file and analyze carbohydrates.")
    parser_arg.add_argument(
        "pdb_file", 
        help="Path to the input PDB file"
    )
    parser_arg.add_argument(
        "--output", "-o", 
        default=None, 
        help="Optional output pickle file name"
    )
    
    args = parser_arg.parse_args()

    # Read PDB content from user-specified file
    with open(args.pdb_file, 'r') as f:
        pdb_content = f.read()
        
    # Create parser and execute analysis
    parser = CarbohydratePDBParser(pdb_content)
    
    # Print summary
    parser.print_summary()
    
    # Get complete analysis
    analysis = parser.generate_analysis()
    
    # Save pickle file
    output_file = args.output if args.output else args.pdb_file.replace(".pdb", "_parser.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(analysis, f)
    
    # Example of data access (unchanged)
    print("\n" + "=" * 60)
    print("EXAMPLE OF AVAILABLE DATA:")
    print("=" * 60)
    
    print(f"\n1. First residue:")
    first_res = analysis['residues'][0]
    print(f"   Name: {first_res['name_pdb']}")
    print(f"   Chain: {first_res['chain_carb']}")
    print(f"   Number: {first_res['res_seq']}")
    print(f"   Atoms: {len(first_res['atoms'])}")
    
    print(f"\n2. Coordinates of first atom:")
    first_atom = analysis['atoms_pdb'][0]
    print(f"   {first_atom['name']}: ({first_atom['x']:.3f}, "
          f"{first_atom['y']:.3f}, {first_atom['z']:.3f})")
    
    print(f"\n3. First intra-residue bond:")
    if analysis['intra_residue_bounds']:
        first_bond = analysis['intra_residue_bounds'][0]
        print(f"   {first_bond['atom1']} -- {first_bond['atom2']}: "
              f"{first_bond['distance']:.3f} Å")
    
    print(f"\n4. First inter-residue bond:")
    if analysis['inter_residue_bounds']:
        first_inter = analysis['inter_residue_bounds'][0]
        print(f"   {first_inter['atom1']} ({first_inter['atom1_res']}) -- "
              f"{first_inter['atom2']} ({first_inter['atom2_res']}): "
              f"{first_inter['distance']:.3f} Å")

