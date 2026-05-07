#!/usr/bin/env python3
"""
Optimize HD22 position in ASN residues by rotating around CA-CB axis.
Goal: Maximize distance between HD22 and neighboring atoms within sphere centered at CA.
Rotation considers the plane defined by [CB, CG, ND2] atoms.
"""

import sys
import math
import numpy as np
from collections import defaultdict
import multiprocessing as mp
from typing import List, Tuple, Dict, Any
import argparse
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Optimize HD22 position in ASN residues')
    parser.add_argument('input_pdb', type=str, 
                       help='Input PDB file path')
    parser.add_argument('-o', '--output', type=str, default='optimized.pdb',
                       help='Output PDB file name (default: optimized.pdb)')
    parser.add_argument('--rotate-atoms', type=str, 
                       default='OD1,CG,ND2,HD22,HD21,HB2,HB3',
                       help='Comma-separated list of atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)')
    parser.add_argument('--fixed-atom', type=str, default='CB',
                       help='Fixed atom for rotation (default: CB)')
    parser.add_argument('--center-atom', type=str, default='CA',
                       help='Center atom for neighbor sphere (default: CA)')
    parser.add_argument('--radius', type=float, default=30.0,
                       help='Radius of neighbor sphere in Å (default: 30.0)')
    parser.add_argument('--rotation-step', type=int, default=1,
                       help='Rotation step in degrees (default: 1)')
    
    return parser.parse_args()

def read_pdb(filename: str) -> Tuple[List[Dict], Dict]:
    """
    Read PDB file and return list of atoms and residue information.
    
    Returns:
        atoms: List of atom dictionaries with keys: line, serial, name, resName, chainID, resSeq, x, y, z
        residues: Dictionary mapping residue numbers to list of atom indices
    """
    atoms = []
    residues = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_info = {
                    'line': line,
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'resName': line[17:20].strip(),
                    'chainID': line[21],
                    'resSeq': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54])
                }
                atoms.append(atom_info)
                res_id = (atom_info['chainID'], atom_info['resSeq'], atom_info['resName'])
                residues[res_id].append(len(atoms) - 1)
    
    return atoms, residues

def find_asn_residues(atoms: List[Dict], residues: Dict, rotate_atoms: List[str], 
                      fixed_atom: str, center_atom: str) -> List[Tuple]:
    """
    Find all ASN residues that have the required atoms for rotation.
    
    Returns:
        List of residue identifiers (chainID, resSeq, resName) for ASN residues
    """
    asn_residues = []
    
    for res_id, atom_indices in residues.items():
        chainID, resSeq, resName = res_id
        
        if resName != "ASN":
            continue
        
        # Check if this ASN has all required atoms
        atom_names = [atoms[idx]['name'] for idx in atom_indices]
        required_atoms = rotate_atoms + [fixed_atom, center_atom]
        
        if all(atom in atom_names for atom in required_atoms):
            asn_residues.append(res_id)
    
    return asn_residues

def extract_coordinates(atoms: List[Dict], atom_indices: List[int], atom_names: List[str]) -> Dict[str, np.ndarray]:
    """
    Extract coordinates for specific atoms in a residue.
    
    Returns:
        Dictionary mapping atom names to coordinates (numpy arrays)
    """
    coords = {}
    for idx in atom_indices:
        atom = atoms[idx]
        if atom['name'] in atom_names:
            coords[atom['name']] = np.array([atom['x'], atom['y'], atom['z']])
    return coords

def get_neighbors(atoms: List[Dict], center: np.ndarray, radius: float, 
                  exclude_indices: List[int]) -> List[np.ndarray]:
    """
    Get coordinates of all atoms within radius of center, excluding specified indices.
    
    Returns:
        List of neighbor atom coordinates
    """
    neighbors = []
    for i, atom in enumerate(atoms):
        if i in exclude_indices:
            continue
        coord = np.array([atom['x'], atom['y'], atom['z']])
        if np.linalg.norm(coord - center) <= radius:
            neighbors.append(coord)
    return neighbors

def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """
    Create a rotation matrix for rotation around an axis by given angle (in radians).
    
    Using Rodrigues' rotation formula.
    """
    axis = axis / np.linalg.norm(axis)
    a = math.cos(angle / 2.0)
    b, c, d = -axis * math.sin(angle / 2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    
    return np.array([
        [aa + bb - cc - dd, 2*(bc + ad), 2*(bd - ac)],
        [2*(bc - ad), aa + cc - bb - dd, 2*(cd + ab)],
        [2*(bd + ac), 2*(cd - ab), aa + dd - bb - cc]
    ])

def apply_rotation(coords: Dict[str, np.ndarray], fixed_point: np.ndarray, 
                   axis: np.ndarray, angle: float) -> Dict[str, np.ndarray]:
    """
    Apply rotation around axis through fixed_point to given coordinates.
    
    Returns:
        Dictionary of rotated coordinates
    """
    rot_mat = rotation_matrix(axis, angle)
    rotated_coords = {}
    
    for atom_name, coord in coords.items():
        # Translate so fixed_point is at origin
        translated = coord - fixed_point
        # Rotate
        rotated = np.dot(rot_mat, translated)
        # Translate back
        rotated_coords[atom_name] = rotated + fixed_point
    
    return rotated_coords

def calculate_average_distance(hd22_coord: np.ndarray, neighbor_coords: List[np.ndarray]) -> float:
    """
    Calculate the average distance from HD22 to all neighbor atoms.
    
    Args:
        hd22_coord: Coordinates of HD22 atom
        neighbor_coords: List of neighbor atom coordinates
        
    Returns:
        Average distance to neighbors
    """
    if not neighbor_coords:
        return 0.0
    
    total_distance = 0.0
    for neighbor in neighbor_coords:
        total_distance += np.linalg.norm(hd22_coord - neighbor)
    
    return total_distance / len(neighbor_coords)

def optimize_asn_residue(args: Tuple) -> Dict:
    """
    Optimize HD22 position for a single ASN residue.
    
    Args:
        args: Tuple containing (res_id, atoms, all_atom_indices, rotate_atoms, 
              fixed_atom, center_atom, radius, rotation_step)
    
    Returns:
        Dictionary with optimization results
    """
    res_id, all_atoms, atom_indices, rotate_atoms, fixed_atom, center_atom, radius, rotation_step = args
    chainID, resSeq, resName = res_id
    
    # Get atom coordinates for this residue
    atom_coords = extract_coordinates(all_atoms, atom_indices, rotate_atoms + [fixed_atom, center_atom])
    
    # Extract key points
    ca_coord = atom_coords[center_atom]
    cb_coord = atom_coords[fixed_atom]
    hd22_original = atom_coords["HD22"]
    
    # Get neighbor atoms (within radius of CA)
    neighbor_coords = get_neighbors(all_atoms, ca_coord, radius, atom_indices)
    
    # Rotation axis: CA -> CB
    axis = cb_coord - ca_coord
    axis = axis / np.linalg.norm(axis)
    
    best_distance = -float('inf')
    best_angle = 0
    best_rotated_coords = None
    
    # Try rotations from 0 to 360 degrees
    for angle_deg in range(0, 360, rotation_step):
        angle_rad = math.radians(angle_deg)
        
        # Rotate the specified atoms
        rotated = apply_rotation({k: v for k, v in atom_coords.items() if k in rotate_atoms}, 
                                 cb_coord, axis, angle_rad)
        
        # Calculate average distance from HD22 to neighbors
        hd22_new = rotated["HD22"]
        
        # Calculate average distance to all neighbors
        avg_dist = calculate_average_distance(hd22_new, neighbor_coords)
        
        # Maximize the average distance
        if avg_dist > best_distance:
            best_distance = avg_dist
            best_angle = angle_deg
            best_rotated_coords = rotated
    
    # Calculate original average distance for comparison
    original_avg_dist = calculate_average_distance(hd22_original, neighbor_coords)
    
    return {
        'residue': res_id,
        'original_coord': hd22_original,
        'new_coord': best_rotated_coords["HD22"],
        'original_avg_dist': original_avg_dist,
        'new_avg_dist': best_distance,
        'best_angle': best_angle,
        'atom_indices': atom_indices,
        'rotated_coords': best_rotated_coords
    }

def update_atom_coordinates(atoms: List[Dict], atom_indices: List[int], 
                           rotated_coords: Dict[str, np.ndarray]) -> None:
    """
    Update atom coordinates in the atoms list with rotated coordinates.
    """
    for idx in atom_indices:
        atom = atoms[idx]
        if atom['name'] in rotated_coords:
            coord = rotated_coords[atom['name']]
            atom['x'], atom['y'], atom['z'] = coord[0], coord[1], coord[2]
            
            # Update the PDB line
            line = atom['line']
            new_line = (line[:30] + 
                       f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}" +
                       line[54:])
            atom['line'] = new_line

def write_pdb(filename: str, atoms: List[Dict]) -> None:
    """
    Write atoms to a PDB file.
    """
    with open(filename, 'w') as f:
        for atom in atoms:
            f.write(atom['line'])

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Parse rotate atoms from comma-separated string
    rotate_atoms = [atom.strip() for atom in args.rotate_atoms.split(',')]
    
    print("="*80)
    print("PARAMETERS:")
    print(f"Input PDB: {args.input_pdb}")
    print(f"Output PDB: {args.output}")
    print(f"Atoms to rotate: {', '.join(rotate_atoms)}")
    print(f"Fixed atom: {args.fixed_atom}")
    print(f"Center atom: {args.center_atom}")
    print(f"Sphere radius: {args.radius} Å")
    print(f"Rotation step: {args.rotation_step}°")
    print("="*80)
    
    # Read PDB file
    print(f"\nReading PDB file: {args.input_pdb}")
    atoms, residues = read_pdb(args.input_pdb)
    
    # Find all ASN residues
    asn_residues = find_asn_residues(atoms, residues, rotate_atoms, args.fixed_atom, args.center_atom)
    print(f"Found {len(asn_residues)} ASN residues to optimize")
    
    if not asn_residues:
        print("No ASN residues found with all required atoms.")
        sys.exit(1)
    
    # Prepare arguments for parallel processing
    args_list = []
    for res_id in asn_residues:
        args_list.append((res_id, atoms, residues[res_id], rotate_atoms, 
                         args.fixed_atom, args.center_atom, args.radius, args.rotation_step))
    
    # Use all available CPUs
    num_cpus = mp.cpu_count()
    print(f"Using {num_cpus} CPUs for parallel processing")
    
    # Process residues in parallel
    with mp.Pool(processes=num_cpus) as pool:
        results = pool.map(optimize_asn_residue, args_list)
    
    # Apply optimizations and print results
    print("\nOptimization results:")
    print("="*120)
    print(f"{'Residue':<10} {'Original HD22':<30} {'New HD22':<30} {'Orig Avg Dist':<12} {'New Avg Dist':<12} {'Angle':<6}")
    print("-"*120)
    
    for result in results:
        res_id = result['residue']
        chainID, resSeq, resName = res_id
        res_str = f"{chainID}:{resName}{resSeq}"
        
        # Print old and new coordinates
        old_coord = result['original_coord']
        new_coord = result['new_coord']
        
        print(f"{res_str:<10} "
              f"({old_coord[0]:7.3f}, {old_coord[1]:7.3f}, {old_coord[2]:7.3f})  "
              f"({new_coord[0]:7.3f}, {new_coord[1]:7.3f}, {new_coord[2]:7.3f})  "
              f"{result['original_avg_dist']:12.3f}  "
              f"{result['new_avg_dist']:12.3f}  "
              f"{result['best_angle']:6}°")
        
        # Update atom coordinates
        update_atom_coordinates(atoms, result['atom_indices'], result['rotated_coords'])
    
    # Write optimized structure
    write_pdb(args.output, atoms)
    print(f"\nOptimized structure written to: {args.output}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY:")
    print(f"Total ASN residues optimized: {len(asn_residues)}")
    print(f"Rotation step used: {args.rotation_step}°")
    print(f"Neighbor sphere radius: {args.radius} Å (centered at {args.center_atom})")
    print(f"Atoms rotated around {args.center_atom}-{args.fixed_atom} axis: {', '.join(rotate_atoms)}")
    print(f"Fixed point for rotation: {args.fixed_atom}")
    print("Objective: Maximize AVERAGE distance from HD22 to atoms within sphere")

if __name__ == "__main__":
    main()
