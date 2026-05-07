"""
Script to optimize glycan orientations using single-axis rotation followed by MCMC refinement.
For each glycan, rotate around the axis from protein atom (ND2 for N-linked, OG/OG1 for O-linked)
to C1 atom of first glycan residue to find best orientation.

Usage:
------

python3 "$SCRIPTS/4-optimize_glycans_mcmc.py" \
    --input_json "$STEP3_RESULTS/JSON_FILES/glycan_data_charmm36.json" \
    --output_json "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/glycan_optimized.json" \
    --output_pdb "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/spike_glycosylated_final_optimized.pdb" \
    --glycans_output_dir "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PBD_CARB_ONLY" \
    --theta_step 10 \
    --n_steps 2 \
    --max_cycles 3 \
    --radius 300 \
    --use_coulomb no \
    --n_workers 12 \
    --report_file "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/report.txt" \
    --save_individual_glycans \
    --save_before_after \
    --verbose

"""

import numpy as np
import json
import argparse
import random
import math
from multiprocessing import Pool, cpu_count, Manager
from tqdm import tqdm
import copy
import sys
import os
from datetime import datetime
from itertools import product
import warnings
warnings.filterwarnings('ignore')

# Global variable for report file
REPORT_FILE = None

def print_and_save(message):
    """Print to console and save to report file."""
    print(message)
    if REPORT_FILE:
        with open(REPORT_FILE, 'a') as f:
            f.write(message + '\n')

# Constants
KB = 0.008314462618  # Boltzmann constant in kJ/mol·K
EPSILON_0 = 8.854187817e-12 * 1e9  # Vacuum permittivity in C²/(N·m²) * 1e9 for nm
E_CHARGE = 1.60217662e-19  # Elementary charge in C
AVOGADRO = 6.02214076e23  # Avogadro's number
ANGSTROM_TO_NM = 0.1  # Conversion factor

def save_pdb_file(data, filename):
    """
    Save all atoms (protein + glycans) to a PDB file with updated coordinates.
    
    Parameters:
    -----------
    data : dict
        Complete system data with updated coordinates
    filename : str
        Path to output PDB file
    """
    print_and_save(f"Saving PDB file to {filename}...")
    
    with open(filename, 'w') as f:
        atom_id = 1
        
        # Write protein atoms
        for atom in data['protein']:
            line = f"ATOM  {atom_id:5d} {atom['atom_name']:4s} {atom['residue_name']:3s} {atom['chain_id']:1s}{atom['residue_number']:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00  0.00           {atom['element']:2s}"
            f.write(line + '\n')
            atom_id += 1
        
        # Write glycan atoms
        for glycan_id, glycan in data['glycans'].items():
            for atom in glycan['atoms']:
                line = f"HETATM{atom_id:5d} {atom['atom_name']:4s} {atom['residue_name']:3s} {atom['chain_id']:1s}{atom['residue_number']:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00  0.00           {atom['element']:2s}"
                f.write(line + '\n')
                atom_id += 1
        
        f.write("END\n")
    
    print_and_save(f"PDB file saved with {atom_id-1} atoms.")

def save_individual_glycans_pdb(data, output_dir):
    """
    Save each glycan separately as individual PDB files with ALL carbohydrate residues.
    
    Parameters:
    -----------
    data : dict
        Complete system data with updated coordinates
    output_dir : str
        Directory to save individual glycan PDB files
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print_and_save(f"Created directory: {output_dir}")
    
    print_and_save(f"\nSaving individual glycans to {output_dir}/")
    
    for glycan_id, glycan in data['glycans'].items():
        filename = os.path.join(output_dir, f"{glycan_id}.pdb")
        
        with open(filename, 'w') as f:
            atom_id = 1
            
            # Write header
            f.write(f"REMARK Individual glycan file: {glycan_id}\n")
            f.write(f"REMARK Generated from optimized coordinates\n")
            f.write(f"REMARK Date: {datetime.now().isoformat()}\n")
            
            # Group atoms by residue to preserve residue numbering
            residue_atoms = {}
            for atom in glycan['atoms']:
                res_key = f"{atom['residue_name']}_{atom['residue_number']}"
                if res_key not in residue_atoms:
                    residue_atoms[res_key] = []
                residue_atoms[res_key].append(atom)
            
            # Write all residues for this glycan
            for res_key in sorted(residue_atoms.keys()):
                residue_name = res_key.split('_')[0]
                residue_number = int(res_key.split('_')[1])
                
                for atom in residue_atoms[res_key]:
                    line = f"HETATM{atom_id:5d} {atom['atom_name']:4s} {residue_name:3s} {atom['chain_id']:1s}{residue_number:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00  0.00           {atom['element']:2s}"
                    f.write(line + '\n')
                    atom_id += 1
            
            f.write("TER\n")
            f.write("END\n")
        
        # Count residues
        residue_count = len(residue_atoms)
        print_and_save(f"  Saved: {glycan_id}.pdb ({atom_id-1} atoms, {residue_count} residues)")

def save_glycan_pdb(glycan_atoms, filename, glycan_id, cycle, stage="before"):
    """
    Save glycan atoms to a PDB file.
    
    Parameters:
    -----------
    glycan_atoms : list
        Atoms of the glycan
    filename : str
        Path to output PDB file
    glycan_id : str
        ID of the glycan
    cycle : int
        Current optimization cycle
    stage : str
        "before" or "after" optimization
    """
    # Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    with open(filename, 'w') as f:
        atom_id = 1
        f.write(f"REMARK Glycan: {glycan_id}\n")
        f.write(f"REMARK Optimization cycle: {cycle}\n")
        f.write(f"REMARK Stage: {stage}\n")
        f.write(f"REMARK Date: {datetime.now().isoformat()}\n")
        
        # Group atoms by residue
        residue_atoms = {}
        for atom in glycan_atoms:
            res_key = f"{atom['residue_name']}_{atom['residue_number']}_{atom['chain_id']}"
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append(atom)
        
        # Write atoms grouped by residue
        for res_key in sorted(residue_atoms.keys()):
            residue_name, residue_number, chain_id = res_key.split('_')
            residue_number = int(residue_number)
            
            for atom in residue_atoms[res_key]:
                line = f"HETATM{atom_id:5d} {atom['atom_name']:4s} {residue_name:3s} {chain_id:1s}{residue_number:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00  0.00           {atom['element']:2s}"
                f.write(line + '\n')
                atom_id += 1
        
        f.write("TER\n")
        f.write("END\n")
    
    print_and_save(f"    Saved {stage} coordinates to: {filename}")

class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for numpy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return super(NumpyEncoder, self).default(obj)

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two points."""
    return np.sqrt(np.sum((np.array(coord1) - np.array(coord2))**2))

def check_collisions(glycan_atoms, all_fixed_atoms, threshold=0.1):
    """
    Check for collisions between glycan atoms and other atoms.
    
    Parameters:
    -----------
    glycan_atoms : list
        Atoms of the glycan
    all_fixed_atoms : list
        All other atoms (protein + other glycans)
    threshold : float
        Distance threshold for collision in nm
    
    Returns:
    --------
    tuple : (has_collisions, collision_count, min_distance)
    """
    collision_count = 0
    min_distance = float('inf')
    
    for atom1 in glycan_atoms:
        coord1 = np.array([atom1['x'], atom1['y'], atom1['z']])
        
        for atom2 in all_fixed_atoms:
            # Skip if same atom or same glycan (if comparing within glycan)
            if atom1 is atom2:
                continue
                
            coord2 = np.array([atom2['x'], atom2['y'], atom2['z']])
            distance = calculate_distance(coord1, coord2)
            
            if distance < min_distance:
                min_distance = distance
            
            if distance < threshold:
                collision_count += 1
    
    has_collisions = collision_count > 0
    return has_collisions, collision_count, min_distance

def rotation_matrix_from_axis_angle(axis, angle):
    """
    Create rotation matrix from axis-angle representation.
    
    Parameters:
    -----------
    axis : np.array
        Rotation axis (unit vector)
    angle : float
        Rotation angle in radians
    
    Returns:
    --------
    np.array : 3x3 rotation matrix
    """
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)
    
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    one_minus_cos = 1 - cos_a
    
    x, y, z = axis
    
    return np.array([
        [cos_a + x*x*one_minus_cos, x*y*one_minus_cos - z*sin_a, x*z*one_minus_cos + y*sin_a],
        [y*x*one_minus_cos + z*sin_a, cos_a + y*y*one_minus_cos, y*z*one_minus_cos - x*sin_a],
        [z*x*one_minus_cos - y*sin_a, z*y*one_minus_cos + x*sin_a, cos_a + z*z*one_minus_cos]
    ])

def rotate_around_axis(points, pivot, axis, angle):
    """
    Rotate points around a single axis through a pivot point.
    
    Parameters:
    -----------
    points : list of np.array
        Points to rotate
    pivot : np.array
        Pivot point (fixed point)
    axis : np.array
        Rotation axis (unit vector)
    angle : float
        Rotation angle in radians
    
    Returns:
    --------
    list : Rotated points
    """
    # Ensure axis is normalized
    axis = axis / np.linalg.norm(axis)
    
    # Create rotation matrix
    R = rotation_matrix_from_axis_angle(axis, angle)
    
    # Apply rotation
    points_array = np.array(points) - pivot
    rotated = np.dot(points_array, R.T)
    
    return (rotated + pivot).tolist()

def get_fixed_point_and_axis(linkage, protein_atoms, glycan_atoms):
    """
    Get the fixed point and rotation axis for glycan rotation.
    For N-linked: axis from ND2 of ASN to C1 of first glycan residue (ND1, ND2, etc.)
    For O-linked (SER): axis from OG of SER to C1 of first glycan residue (A21, A22, etc.)
    For O-linked (THR): axis from OG1 of THR to C1 of first glycan residue (A21, A22, etc.)
    
    Parameters:
    -----------
    linkage : dict
        Linkage information
    protein_atoms : list
        All protein atoms
    glycan_atoms : list
        Glycan atoms
    
    Returns:
    --------
    tuple : (fixed_point, rotation_axis, c1_atom, c1_point, fixed_atom)
    """
    protein_res_num = linkage['protein_residue_number']
    protein_chain = linkage['protein_chain']
    protein_atom_name = linkage['protein_atom']
    linking_type = linkage['linking_type']
    
    # Find the protein atom (fixed point)
    fixed_atom = None
    for atom in protein_atoms:
        if (atom['residue_number'] == protein_res_num and
            atom['chain_id'] == protein_chain and
            atom['atom_name'].strip() == protein_atom_name):
            
            fixed_atom = atom
            fixed_point = np.array([atom['x'], atom['y'], atom['z']])
            break
    
    if fixed_atom is None:
        # If not found, use the coordinates from linkage
        if 'protein_atom_complete' in linkage:
            atom = linkage['protein_atom_complete']
            fixed_point = np.array([atom['x'], atom['y'], atom['z']])
            fixed_atom = atom
        else:
            raise ValueError(f"Fixed point atom not found for linkage: {linkage['glycan_binding']}")
    
    # Find C1 atom of first glycan residue
    c1_atom = None
    
    # First, identify the first residue of the glycan
    # Sort atoms by residue number to find the first residue
    sorted_atoms = sorted(glycan_atoms, key=lambda x: x['residue_number'])
    
    if sorted_atoms:
        first_residue_number = sorted_atoms[0]['residue_number']
        # Find C1 atom in the first residue
        for atom in sorted_atoms:
            if atom['residue_number'] == first_residue_number and atom['atom_name'].strip() == 'C1':
                c1_atom = atom
                c1_point = np.array([atom['x'], atom['y'], atom['z']])
                break
    
    if c1_atom is None:
        # If C1 not found, use first atom as proxy
        if sorted_atoms:
            c1_atom = sorted_atoms[0]
            c1_point = np.array([c1_atom['x'], c1_atom['y'], c1_atom['z']])
        else:
            raise ValueError(f"C1 atom not found for glycan")
    
    # Rotation axis: from fixed point to C1 atom
    rotation_axis = c1_point - fixed_point
    if np.linalg.norm(rotation_axis) < 0.1:
        # Default axis if C1 too close to fixed point
        rotation_axis = np.array([1, 0, 0])
    else:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    
    return fixed_point, rotation_axis, c1_atom, c1_point, fixed_atom

def calculate_vdw_energy(atom1, atom2, distance):
    """
    Calculate Van der Waals energy (Lennard-Jones) between two atoms.
    
    Parameters:
    -----------
    atom1, atom2 : dict
        Atom information with CHARMM sigma and epsilon
    distance : float
        Distance between atoms in nm
    
    Returns:
    --------
    float : Van der Waals energy in kJ/mol
    """
    # Penalidade maior para colisões fortes
    if distance < 0.05:  # Colisão muito forte (0.5 Å)
        return 10000.0
    elif distance < 0.1:  # Colisão moderada (1.0 Å)
        return 5000.0
    elif distance < 0.15:  # Colisão leve (1.5 Å)
        return 1000.0
    
    sigma1 = atom1.get('charmm_sigma', 0.1)
    epsilon1 = atom1.get('charmm_epsilon', 0.0)
    sigma2 = atom2.get('charmm_sigma', 0.1)
    epsilon2 = atom2.get('charmm_epsilon', 0.0)
    
    # Lorentz-Berthelot mixing rules
    sigma_ij = (sigma1 + sigma2) / 2
    epsilon_ij = np.sqrt(epsilon1 * epsilon2)
    
    if epsilon_ij == 0 or sigma_ij == 0:
        return 0.0
    
    # Repulsão mais forte em distâncias curtas
    if distance < sigma_ij * 0.8:  # Dentro de 80% do sigma combinado
        repulsion_factor = (sigma_ij * 0.8 / distance) ** 12
        return 10.0 * epsilon_ij * repulsion_factor
    
    r6 = (sigma_ij / distance) ** 6
    r12 = r6 * r6
    
    # Lennard-Jones potential: 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
    energy = 4 * epsilon_ij * (r12 - r6)
    
    return energy

def calculate_coulomb_energy(atom1, atom2, distance):
    """
    Calculate Coulomb energy between two atoms.
    
    Parameters:
    -----------
    atom1, atom2 : dict
        Atom information with CHARMM charge
    distance : float
        Distance between atoms in nm
    
    Returns:
    --------
    float : Coulomb energy in kJ/mol
    """
    if distance < 0.1:  # Avoid singularity
        distance = 0.1
    
    q1 = atom1.get('charmm_charge', 0.0)
    q2 = atom2.get('charmm_charge', 0.0)
    
    # Convert charges to Coulomb
    q1_c = q1 * E_CHARGE
    q2_c = q2 * E_CHARGE
    
    # Calculate energy in J, then convert to kJ/mol
    # Formula: E = (q1 * q2) / (4 * π * ε0 * r)
    energy_j = (q1_c * q2_c) / (4 * math.pi * EPSILON_0 * distance * 1e-9)
    energy_kjmol = energy_j * AVOGADRO / 1000
    
    return energy_kjmol

def get_atoms_within_sphere(center, radius_nm, all_atoms):
    """
    Get all atoms within a sphere of given radius from center.
    
    Parameters:
    -----------
    center : np.array
        Center of sphere in nm
    radius_nm : float
        Radius of sphere in nm
    all_atoms : list
        List of all atoms to check
    
    Returns:
    --------
    list : Atoms within the sphere
    """
    atoms_within = []
    for atom in all_atoms:
        atom_coord = np.array([atom['x'], atom['y'], atom['z']])
        if calculate_distance(center, atom_coord) <= radius_nm:
            atoms_within.append(atom)
    return atoms_within

def calculate_local_energy(glycan_atoms, all_fixed_atoms, center, radius_nm, use_coulomb=False):
    """
    Calculate energy considering ALL atoms within the spherical region.
    
    Parameters:
    -----------
    glycan_atoms : list
        Atoms of the glycan being optimized
    all_fixed_atoms : list
        All atoms in the system (protein + other glycans)
    center : np.array
        Center of the sphere (usually glycan COM)
    radius_nm : float
        Radius of the spherical region in nm
    use_coulomb : bool
        Whether to include Coulomb energy
    
    Returns:
    --------
    tuple : (total_energy, vdw_energy, coulomb_energy)
    """
    if not glycan_atoms or not all_fixed_atoms:
        return 0.0, 0.0, 0.0
    
    total_vdw = 0.0
    total_coulomb = 0.0
    
    # Get all atoms within sphere (including current glycan atoms)
    all_atoms_within = get_atoms_within_sphere(center, radius_nm, all_fixed_atoms)
    
    if not all_atoms_within:
        return 0.0, 0.0, 0.0
    
    # Counters for debugging
    coulomb_interactions = 0
    nonzero_coulomb = 0
    
    # Calculate energies between glycan atoms and all other atoms within sphere
    for atom1 in glycan_atoms:
        coord1 = np.array([atom1['x'], atom1['y'], atom1['z']])
        
        for atom2 in all_atoms_within:
            coord2 = np.array([atom2['x'], atom2['y'], atom2['z']])
            distance = calculate_distance(coord1, coord2)
            
            # Skip very small distances (same atom or covalent bonds)
            if distance < 0.01:
                continue
            
            # Van der Waals energy (always calculated)
            vdw_energy = calculate_vdw_energy(atom1, atom2, distance)
            total_vdw += vdw_energy
            
            # Coulomb energy (optional)
            if use_coulomb:
                coulomb_energy = calculate_coulomb_energy(atom1, atom2, distance)
                total_coulomb += coulomb_energy
                
                # Debug counters
                coulomb_interactions += 1
                if abs(coulomb_energy) > 1e-10:  # Small threshold for "nonzero"
                    nonzero_coulomb += 1
    
    # Debug information
    if use_coulomb and total_coulomb != 0:
        print_and_save(f"    Coulomb energy calculation:")
        print_and_save(f"      Total Coulomb interactions: {coulomb_interactions}")
        print_and_save(f"      Non-zero Coulomb interactions: {nonzero_coulomb}")
        print_and_save(f"      Total Coulomb energy: {total_coulomb:.6f} kJ/mol")
    
    total_energy = total_vdw + total_coulomb
    
    return total_energy, total_vdw, total_coulomb

def calculate_center_of_mass(glycan_atoms):
    """Calculate center of mass of ALL glycan atoms."""
    total_mass = 0.0
    com = np.zeros(3)
    
    for atom in glycan_atoms:
        mass = atom.get('charmm_mass', 1.0)
        coord = np.array([atom['x'], atom['y'], atom['z']])
        com += mass * coord
        total_mass += mass
    
    if total_mass > 0:
        com /= total_mass
    
    return com

def get_glycan_atom_coordinates(glycan_atoms):
    """Extract coordinates from glycan atoms."""
    return [np.array([atom['x'], atom['y'], atom['z']]) for atom in glycan_atoms]

def update_glycan_atom_coordinates(glycan_atoms, new_coords):
    """Update coordinates in glycan atoms."""
    for i, atom in enumerate(glycan_atoms):
        atom['x'] = new_coords[i][0]
        atom['y'] = new_coords[i][1]
        atom['z'] = new_coords[i][2]
    return glycan_atoms

def single_axis_grid_search_worker(args):
    """
    Worker function for parallel single-axis grid search.
    
    Parameters:
    -----------
    args : tuple
        (glycan_atoms, all_fixed_atoms, fixed_point, rotation_axis,
         c1_point, theta_angles, center, radius_nm, use_coulomb)
    
    Returns:
    --------
    tuple : (theta, energy, vdw_energy, coulomb_energy)
    """
    (glycan_atoms, all_fixed_atoms, fixed_point, rotation_axis,
     c1_point, theta_angles, center, radius_nm, use_coulomb) = args
    
    best_energy = float('inf')
    best_theta = 0.0
    best_vdw = 0.0
    best_coulomb = 0.0
    
    original_coords = get_glycan_atom_coordinates(glycan_atoms)
    
    for theta in theta_angles:
        # Convert to radians
        theta_rad = np.deg2rad(theta)
        
        # Apply rotation around single axis
        new_coords = rotate_around_axis(
            original_coords, c1_point, rotation_axis, theta_rad
        )
        
        # Create temporary atoms with new coordinates
        temp_atoms = copy.deepcopy(glycan_atoms)
        temp_atoms = update_glycan_atom_coordinates(temp_atoms, new_coords)
        
        # Calculate energy with ALL atoms within sphere
        energy, vdw_energy, coulomb_energy = calculate_local_energy(
            temp_atoms, all_fixed_atoms, center, radius_nm, use_coulomb
        )
        
        if energy < best_energy:
            best_energy = energy
            best_theta = theta
            best_vdw = vdw_energy
            best_coulomb = coulomb_energy
    
    return best_theta, best_energy, best_vdw, best_coulomb

def parallel_single_axis_grid_search(glycan_atoms, all_fixed_atoms, fixed_point, 
                                     rotation_axis, c1_point, theta_range,
                                     center, radius_nm, use_coulomb, n_workers):
    """
    Perform parallel grid search for best theta angle (single axis rotation).
    
    Parameters:
    -----------
    glycan_atoms : list
        Glycan atoms
    all_fixed_atoms : list
        All fixed atoms in the system (protein + other glycans)
    fixed_point : np.array
        Pivot point for rotation
    rotation_axis : np.array
        Rotation axis (from protein atom to C1 of first glycan residue)
    c1_point : np.array
        C1 atom coordinates (rotation pivot)
    theta_range : list
        List of theta angles to test (0-360 degrees)
    center : np.array
        Center of spherical region
    radius_nm : float
        Radius of spherical region in nm
    use_coulomb : bool
        Whether to include Coulomb energy
    n_workers : int
        Number of parallel workers
    
    Returns:
    --------
    tuple : (best_theta, best_energy, best_vdw, best_coulomb)
    """
    # Split theta_range into chunks for parallel processing
    chunk_size = max(1, len(theta_range) // n_workers)
    theta_chunks = [theta_range[i:i + chunk_size] for i in range(0, len(theta_range), chunk_size)]
    
    # Prepare arguments for workers
    args_list = []
    for theta_chunk in theta_chunks:
        args_list.append((
            glycan_atoms, all_fixed_atoms, fixed_point, rotation_axis,
            c1_point, theta_chunk, center, radius_nm, use_coulomb
        ))
    
    # Run parallel grid search
    with Pool(min(n_workers, len(args_list))) as pool:
        results = pool.map(single_axis_grid_search_worker, args_list)
    
    # Find overall best result
    best_energy = float('inf')
    best_theta = 0.0
    best_vdw = 0.0
    best_coulomb = 0.0
    
    for theta, energy, vdw, coulomb in results:
        if energy < best_energy:
            best_energy = energy
            best_theta = theta
            best_vdw = vdw
            best_coulomb = coulomb
    
    return best_theta, best_energy, best_vdw, best_coulomb

def refine_with_mcmc_single_axis(glycan_atoms, all_fixed_atoms, fixed_point, 
                                 rotation_axis, c1_point, initial_theta,
                                 center, radius_nm, use_coulomb, n_steps=10000, 
                                 temperature=300, step_size=10, verbose=True):
    """
    Refine orientation using MCMC around initial theta angle (single axis).
    
    Parameters:
    -----------
    glycan_atoms : list
        Glycan atoms
    all_fixed_atoms : list
        All fixed atoms in the system (protein + other glycans)
    fixed_point : np.array
        Pivot point for rotation
    rotation_axis : np.array
        Rotation axis
    c1_point : np.array
        C1 atom coordinates (rotation pivot)
    initial_theta : float
        Initial theta angle in degrees
    center : np.array
        Center of spherical region
    radius_nm : float
        Radius of spherical region in nm
    use_coulomb : bool
        Whether to include Coulomb energy
    n_steps : int
        Number of MCMC steps
    temperature : float
        Temperature for Boltzmann acceptance
    step_size : float
        Maximum step size in degrees
    verbose : bool
        Print progress information
    
    Returns:
    --------
    tuple : (final_theta, final_energy, energy_history, acceptance_rate)
    """
    current_theta = initial_theta
    
    # Get original coordinates
    original_coords = get_glycan_atom_coordinates(glycan_atoms)
    
    # Calculate initial energy
    theta_rad = np.deg2rad(current_theta)
    current_coords = rotate_around_axis(
        original_coords, c1_point, rotation_axis, theta_rad
    )
    temp_atoms = copy.deepcopy(glycan_atoms)
    temp_atoms = update_glycan_atom_coordinates(temp_atoms, current_coords)
    
    current_energy, current_vdw, current_coulomb = calculate_local_energy(
        temp_atoms, all_fixed_atoms, center, radius_nm, use_coulomb
    )
    
    best_theta = current_theta
    best_energy = current_energy
    
    energy_history = [current_energy]
    acceptance_history = []
    
    step_size_rad = np.deg2rad(step_size)
    
    if verbose:
        print_and_save(f"    MCMC refinement: {n_steps} steps, step_size={step_size}°")
        pbar = tqdm(total=n_steps, desc="MCMC refinement", leave=False)
    
    for step in range(n_steps):
        # Propose new angle
        delta_theta = (np.random.random() - 0.5) * 2 * step_size_rad
        
        new_theta_rad = np.deg2rad(current_theta) + delta_theta
        
        # Apply rotation
        new_coords = rotate_around_axis(
            original_coords, c1_point, rotation_axis, new_theta_rad
        )
        
        # Create temporary atoms
        temp_atoms = copy.deepcopy(glycan_atoms)
        temp_atoms = update_glycan_atom_coordinates(temp_atoms, new_coords)
        
        # Calculate new energy
        new_energy, new_vdw, new_coulomb = calculate_local_energy(
            temp_atoms, all_fixed_atoms, center, radius_nm, use_coulomb
        )
        
        # Metropolis acceptance criterion
        delta_energy = new_energy - current_energy
        
        if delta_energy < 0:
            accept = True
        else:
            probability = np.exp(-delta_energy / (KB * temperature))
            accept = np.random.random() < probability
        
        if accept:
            current_theta = np.rad2deg(new_theta_rad)
            current_energy = new_energy
            current_coords = new_coords
            acceptance_history.append(1)
            
            # Update best if improved
            if new_energy < best_energy:
                best_theta = current_theta
                best_energy = new_energy
        else:
            acceptance_history.append(0)
        
        energy_history.append(current_energy)
        
        if verbose:
            pbar.update(1)
            pbar.set_postfix({
                'Energy': f'{current_energy:.2f}',
                'Best': f'{best_energy:.2f}',
                'Accept': f'{np.mean(acceptance_history[-100:] if len(acceptance_history) >= 100 else acceptance_history):.2f}'
            })
    
    if verbose:
        pbar.close()
    
    acceptance_rate = np.mean(acceptance_history) if acceptance_history else 0.0
    
    return best_theta, best_energy, energy_history, acceptance_rate

def optimize_glycan_single_axis(data, glycan_id, linkage, theta_step=30,
                                n_steps=10000, radius_angstrom=80, use_coulomb=False,
                                n_workers=4, verbose=True, cycle=1, 
                                save_before_after=False, output_dir=None):
    """
    Optimize a single glycan using single-axis grid search + MCMC.
    
    Parameters:
    -----------
    data : dict
        Complete system data
    glycan_id : str
        ID of glycan to optimize
    linkage : dict
        Linkage information
    theta_step : int
        Step size for theta grid search (single axis)
    n_steps : int
        Number of MCMC refinement steps
    radius_angstrom : float
        Radius of spherical region in Angstroms
    use_coulomb : bool
        Whether to include Coulomb energy
    n_workers : int
        Number of parallel workers
    verbose : bool
        Print progress information
    cycle : int
        Current optimization cycle number
    save_before_after : bool
        Whether to save PDB files before and after optimization
    output_dir : str
        Directory to save PDB files
    
    Returns:
    --------
    tuple : (optimized_data, optimization_info)
    """
    if verbose:
        print_and_save(f"\n[{cycle}] Optimizing glycan: {glycan_id}")
        print_and_save(f"  Linkage: {linkage['site_protein_residue']} -> {glycan_id}")
        print_and_save(f"  Type: {linkage['linking_type']}")
    
    # Get glycan data
    glycan_data = data['glycans'][glycan_id]
    glycan_atoms = glycan_data['atoms']
    
    # Count residues
    residues = set()
    for atom in glycan_atoms:
        residues.add((atom['residue_name'], atom['residue_number']))
    
    if verbose:
        print_and_save(f"  Glycan has {len(residues)} residues: {sorted(residues)[:5]}{'...' if len(residues) > 5 else ''}")
    
    # Check if atoms have CHARMM parameters
    missing_params = []
    for atom in glycan_atoms:
        if 'charmm_charge' not in atom:
            missing_params.append(f"{atom['atom_name']}")
    
    if missing_params and verbose:
        print_and_save(f"  WARNING: {len(missing_params)} atoms missing CHARMM parameters: {missing_params[:5]}{'...' if len(missing_params) > 5 else ''}")
        print_and_save(f"  Using default values (charge=0.0, epsilon=0.0, sigma=0.1)")
    
    # Get ALL fixed atoms (protein + other glycans) - IMPORTANT: all atoms
    protein_atoms = data['protein']
    all_fixed_atoms = protein_atoms.copy()
    for other_id, other_glycan in data['glycans'].items():
        if other_id != glycan_id:
            all_fixed_atoms.extend(other_glycan['atoms'])
    
    # Check collisions before optimization
    has_collisions_before, collision_count_before, min_dist_before = check_collisions(
        glycan_atoms, all_fixed_atoms, threshold=0.1  # 1.0 Å threshold
    )
    
    if verbose:
        if has_collisions_before:
            print_and_save(f"  WARNING: {collision_count_before} collisions detected before optimization!")
            print_and_save(f"  Minimum distance: {min_dist_before*10:.2f} Å")
        else:
            print_and_save(f"  No collisions detected before optimization")
            print_and_save(f"  Minimum distance: {min_dist_before*10:.2f} Å")
    
    # Save before PDB if requested
    if save_before_after and output_dir:
        before_filename = os.path.join(output_dir, f"{glycan_id}_before_{cycle}.pdb")
        save_glycan_pdb(glycan_atoms, before_filename, glycan_id, cycle, "before")
    
    # Get fixed point and rotation axis
    fixed_point, rotation_axis, c1_atom, c1_point, fixed_atom = get_fixed_point_and_axis(
        linkage, protein_atoms, glycan_atoms
    )
    
    if verbose:
        print_and_save(f"  Fixed point: {fixed_atom['residue_name']}{fixed_atom['residue_number']}"
              f".{fixed_atom['atom_name']} at ({fixed_point[0]:.2f}, {fixed_point[1]:.2f}, {fixed_point[2]:.2f})")
        print_and_save(f"  First glycan residue C1: {c1_atom['residue_name']}{c1_atom['residue_number']}"
              f".{c1_atom['atom_name']} at ({c1_point[0]:.2f}, {c1_point[1]:.2f}, {c1_point[2]:.2f})")
    
    # Calculate center of mass of ALL glycan atoms
    com = calculate_center_of_mass(glycan_atoms)
    radius_nm = radius_angstrom * ANGSTROM_TO_NM
    
    if verbose:
        print_and_save(f"  Center of mass of glycan: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})")
        print_and_save(f"  Search radius: {radius_angstrom} Å ({radius_nm:.2f} nm)")
        print_and_save(f"  Coulomb energy: {'INCLUDED' if use_coulomb else 'EXCLUDED'}")
    
    # Calculate and print distance and vectors before optimization
    distance_fixed_to_com = calculate_distance(fixed_point, com)
    distance_fixed_to_c1 = calculate_distance(fixed_point, c1_point)
    vector_c1_to_com = com - c1_point
    
    print_and_save(f"\n  VECTORS AND DISTANCES BEFORE OPTIMIZATION:")
    print_and_save(f"    Distance fixed point to C1: {distance_fixed_to_c1:.2f} Å")
    print_and_save(f"    Distance fixed point to COM: {distance_fixed_to_com:.2f} Å")
    print_and_save(f"    Vector C1 -> COM: [{vector_c1_to_com[0]:.2f}, {vector_c1_to_com[1]:.2f}, {vector_c1_to_com[2]:.2f}]")
    print_and_save(f"    Vector magnitude (C1->COM): {np.linalg.norm(vector_c1_to_com):.2f} Å")
    print_and_save(f"    Rotation axis (fixed->C1): [{rotation_axis[0]:.2f}, {rotation_axis[1]:.2f}, {rotation_axis[2]:.2f}]")
    
    # Define angle range for single-axis grid search
    theta_range = list(range(0, 360, theta_step))
    
    if verbose:
        print_and_save(f"\n  Single-axis grid search: {len(theta_range)} angles")
        print_and_save(f"  Theta step: {theta_step}°")
        print_and_save(f"  Parallel workers: {n_workers}")
    
    # Initial coordinates
    initial_coords = get_glycan_atom_coordinates(glycan_atoms)
    initial_energy, initial_vdw, initial_coulomb = calculate_local_energy(
        glycan_atoms, all_fixed_atoms, com, radius_nm, use_coulomb
    )
    
    if verbose:
        print_and_save(f"\n  Initial energy: {initial_energy:.2f} kJ/mol")
        print_and_save(f"    VdW: {initial_vdw:.2f} kJ/mol")
        if use_coulomb:
            print_and_save(f"    Coulomb: {initial_coulomb:.2f} kJ/mol")
    
    # Perform parallel single-axis grid search
    best_theta, grid_energy, grid_vdw, grid_coulomb = parallel_single_axis_grid_search(
        glycan_atoms, all_fixed_atoms, fixed_point, rotation_axis, c1_point,
        theta_range, com, radius_nm, use_coulomb, n_workers
    )
    
    if verbose:
        print_and_save(f"\n  Grid search results:")
        print_and_save(f"    Best theta angle: {best_theta}°")
        print_and_save(f"    Best energy: {grid_energy:.2f} kJ/mol")
        print_and_save(f"      VdW: {grid_vdw:.2f} kJ/mol")
        if use_coulomb:
            print_and_save(f"      Coulomb: {grid_coulomb:.2f} kJ/mol")
        print_and_save(f"    Energy improvement: {initial_energy - grid_energy:.2f} kJ/mol")
    
    # MCMC refinement
    final_theta, final_energy, energy_history, acceptance_rate = refine_with_mcmc_single_axis(
        glycan_atoms, all_fixed_atoms, fixed_point, rotation_axis, c1_point,
        best_theta, com, radius_nm, use_coulomb, n_steps=n_steps,
        temperature=300, step_size=10, verbose=verbose
    )
    
    if verbose:
        print_and_save(f"\n  MCMC refinement results:")
        print_and_save(f"    Final theta angle: {final_theta:.1f}°")
        print_and_save(f"    Final energy: {final_energy:.2f} kJ/mol")
        print_and_save(f"    Acceptance rate: {acceptance_rate:.3f}")
        print_and_save(f"    Total improvement: {initial_energy - final_energy:.2f} kJ/mol")
    
    # Apply final rotation ONLY if final energy is lower than initial energy
    if final_energy < initial_energy:
        final_theta_rad = np.deg2rad(final_theta)
        final_coords = rotate_around_axis(
            initial_coords, c1_point, rotation_axis, final_theta_rad
        )
        
        # Update glycan coordinates
        glycan_data['atoms'] = update_glycan_atom_coordinates(glycan_atoms, final_coords)
        data['glycans'][glycan_id] = glycan_data
        
        # Check collisions after optimization
        has_collisions_after, collision_count_after, min_dist_after = check_collisions(
            glycan_data['atoms'], all_fixed_atoms, threshold=0.1
        )
        
        if verbose:
            if has_collisions_after:
                print_and_save(f"  WARNING: {collision_count_after} collisions detected after optimization!")
                print_and_save(f"  Minimum distance: {min_dist_after*10:.2f} Å")
            else:
                print_and_save(f"  No collisions detected after optimization")
                print_and_save(f"  Minimum distance: {min_dist_after*10:.2f} Å")
            
            # Report improvement in collisions
            if has_collisions_before or has_collisions_after:
                collision_diff = collision_count_before - collision_count_after
                if collision_diff > 0:
                    print_and_save(f"  ✓ Reduced collisions by {collision_diff}")
                elif collision_diff < 0:
                    print_and_save(f"  ✗ Increased collisions by {-collision_diff}")
                else:
                    print_and_save(f"  → Collision count unchanged")
        
        # Calculate new COM and C1 position after coordinate update
        new_com = calculate_center_of_mass(glycan_data['atoms'])
        
        # Find new C1 position
        new_c1_point = None
        for atom in glycan_data['atoms']:
            if atom['residue_number'] == c1_atom['residue_number'] and atom['atom_name'].strip() == 'C1':
                new_c1_point = np.array([atom['x'], atom['y'], atom['z']])
                break
        
        if new_c1_point is None:
            # Fallback: use first atom
            new_c1_point = np.array([glycan_data['atoms'][0]['x'], 
                                     glycan_data['atoms'][0]['y'], 
                                     glycan_data['atoms'][0]['z']])
        
        # Calculate distances and vectors after rotation
        distance_fixed_to_com_after = calculate_distance(fixed_point, new_com)
        distance_fixed_to_c1_after = calculate_distance(fixed_point, new_c1_point)
        new_vector_c1_to_com = new_com - new_c1_point
        
        print_and_save(f"\n  COORDINATE CHANGE APPLIED (final energy < initial energy)")
        print_and_save(f"  VECTORS AND DISTANCES AFTER COORDINATE CHANGE:")
        print_and_save(f"    Distance fixed point to C1: {distance_fixed_to_c1_after:.2f} Å")
        print_and_save(f"    Distance fixed point to COM: {distance_fixed_to_com_after:.2f} Å")
        print_and_save(f"    Vector C1 -> COM: [{new_vector_c1_to_com[0]:.2f}, {new_vector_c1_to_com[1]:.2f}, {new_vector_c1_to_com[2]:.2f}]")
        print_and_save(f"    Vector magnitude (C1->COM): {np.linalg.norm(new_vector_c1_to_com):.2f} Å")
        print_and_save(f"    Distance change (fixed->C1): {distance_fixed_to_c1_after - distance_fixed_to_c1:.2f} Å")
        print_and_save(f"    Vector change (C1->COM): [{new_vector_c1_to_com[0]-vector_c1_to_com[0]:.2f}, {new_vector_c1_to_com[1]-vector_c1_to_com[1]:.2f}, {new_vector_c1_to_com[2]-vector_c1_to_com[2]:.2f}]")
        
        # Verify that distance from fixed point to C1 doesn't change
        distance_change = abs(distance_fixed_to_c1_after - distance_fixed_to_c1)
        if distance_change > 0.01:
            print_and_save(f"  WARNING: Distance from fixed point to C1 changed by {distance_change:.4f} Å")
        else:
            print_and_save(f"  OK: Distance from fixed point to C1 maintained (change: {distance_change:.4f} Å)")
    else:
        print_and_save(f"\n  COORDINATE CHANGE NOT APPLIED (final energy >= initial energy)")
        print_and_save(f"  Keeping original coordinates.")
        new_com = com  # No change in COM
        new_vector_c1_to_com = vector_c1_to_com
        has_collisions_after = has_collisions_before
        collision_count_after = collision_count_before
        min_dist_after = min_dist_before
    
    # Save after PDB if requested
    if save_before_after and output_dir:
        after_filename = os.path.join(output_dir, f"{glycan_id}_after_{cycle}.pdb")
        save_glycan_pdb(glycan_data['atoms'], after_filename, glycan_id, cycle, "after")
    
    # Store optimization info with proper type conversions
    optimization_info = {
        'glycan_id': glycan_id,
        'cycle': int(cycle),
        'initial_energy': float(initial_energy),
        'initial_vdw': float(initial_vdw),
        'initial_coulomb': float(initial_coulomb),
        'grid_best_theta': float(best_theta),
        'grid_energy': float(grid_energy),
        'grid_vdw': float(grid_vdw),
        'grid_coulomb': float(grid_coulomb),
        'final_theta': float(final_theta),
        'final_energy': float(final_energy),
        'acceptance_rate': float(acceptance_rate),
        'improvement': float(initial_energy - final_energy),
        'collisions_before': int(collision_count_before),
        'collisions_after': int(collision_count_after),
        'min_distance_before': float(min_dist_before * 10),  # Convert to Å
        'min_distance_after': float(min_dist_after * 10),    # Convert to Å
        'collision_reduction': int(collision_count_before - collision_count_after),
        'fixed_point': [float(x) for x in fixed_point.tolist()],
        'rotation_axis': [float(x) for x in rotation_axis.tolist()],
        'c1_point_initial': [float(x) for x in c1_point.tolist()],
        'initial_center_of_mass': [float(x) for x in com.tolist()],
        'final_center_of_mass': [float(x) for x in new_com.tolist()],
        'initial_vector_c1_to_com': [float(x) for x in vector_c1_to_com.tolist()],
        'initial_vector_magnitude': float(np.linalg.norm(vector_c1_to_com)),
        'final_vector_c1_to_com': [float(x) for x in new_vector_c1_to_com.tolist()],
        'final_vector_magnitude': float(np.linalg.norm(new_vector_c1_to_com)),
        'distance_fixed_to_c1_initial': float(distance_fixed_to_c1),
        'distance_fixed_to_c1_final': float(distance_fixed_to_c1_after if final_energy < initial_energy else distance_fixed_to_c1),
        'coordinate_change_applied': bool(final_energy < initial_energy),
        'linking_type': linkage['linking_type'],
        'residue_count': len(residues)
    }
    
    # Print coordinate comparison
    if verbose:
        print_and_save(f"\n  Coordinate comparison:")
        print_and_save(f"    Initial COM: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})")
        print_and_save(f"    Final COM:   ({new_com[0]:.2f}, {new_com[1]:.2f}, {new_com[2]:.2f})")
        
        # Print first atom coordinates
        if glycan_data['atoms']:
            first_atom = glycan_data['atoms'][0]
            print_and_save(f"    First atom {first_atom['atom_name']}:")
            print_and_save(f"      Initial: ({initial_coords[0][0]:.2f}, {initial_coords[0][1]:.2f}, {initial_coords[0][2]:.2f})")
            print_and_save(f"      Final:   ({first_atom['x']:.2f}, {first_atom['y']:.2f}, {first_atom['z']:.2f})")
    
    return data, optimization_info

def optimize_all_glycans_single_axis(data, theta_step=30, n_steps=10000,
                                     max_cycles=5, radius_angstrom=80, use_coulomb=False,
                                     n_workers=4, verbose=True, report_file=None,
                                     save_before_after=True, output_pdb=None,
                                     glycans_output_dir=None):
    """
    Optimize all glycans in the system using single-axis rotation.
    
    Parameters:
    -----------
    data : dict
        Complete system data
    theta_step : int
        Step size for theta grid search (single axis)
    n_steps : int
        Number of MCMC refinement steps
    max_cycles : int
        Maximum optimization cycles
    radius_angstrom : float
        Radius of spherical region in Angstroms
    use_coulomb : bool
        Whether to include Coulomb energy
    n_workers : int
        Number of parallel workers
    verbose : bool
        Print progress information
    report_file : str
        Path to report file
    save_before_after : bool
        Whether to save PDB files before and after optimization
    output_pdb : str
        Path to output PDB file for the complex
    glycans_output_dir : str
        Directory to save individual glycan PDB files
    
    Returns:
    --------
    dict : Optimized system data
    """
    global REPORT_FILE
    REPORT_FILE = report_file
    
    # Create necessary directories
    if save_before_after and glycans_output_dir:
        os.makedirs(glycans_output_dir, exist_ok=True)
        print_and_save(f"Created directory for individual glycans: {glycans_output_dir}")
    
    # Save initial complex PDB
    if output_pdb:
        initial_pdb = output_pdb.replace('.pdb', '_initial.pdb')
        save_pdb_file(data, initial_pdb)
        print_and_save(f"Saved initial complex to: {initial_pdb}")
    
    # Initialize report file
    if REPORT_FILE:
        with open(REPORT_FILE, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("GLYCAN ORIENTATION OPTIMIZATION - SINGLE AXIS ROTATION\n")
            f.write("=" * 80 + "\n")
            f.write(f"Report file created: {datetime.now().isoformat()}\n")
            f.write("=" * 80 + "\n\n")
    
    print_and_save("=" * 80)
    print_and_save("GLYCAN ORIENTATION OPTIMIZATION - SINGLE AXIS ROTATION")
    print_and_save("=" * 80)
    print_and_save(f"Number of glycans: {len(data['glycans'])}")
    print_and_save(f"Single-axis rotation: theta_step={theta_step}°")
    print_and_save(f"MCMC refinement steps: {n_steps}")
    print_and_save(f"Maximum cycles: {max_cycles}")
    print_and_save(f"Search radius: {radius_angstrom} Å")
    print_and_save(f"Include Coulomb: {use_coulomb}")
    print_and_save(f"Parallel workers: {n_workers}")
    print_and_save(f"Save before/after PDBs: {save_before_after}")
    if save_before_after:
        print_and_save(f"Glycans output directory: {glycans_output_dir}")
    print_and_save(f"Report file: {report_file if report_file else 'Not specified'}")
    print_and_save("=" * 80)
    
    # Create a copy of data
    optimized_data = copy.deepcopy(data)
    
    # Get linkage information
    glycan_linkages = {}
    for linkage in optimized_data['linkages']:
        glycan_id = linkage['glycan_binding']
        if glycan_id in optimized_data['glycans']:
            glycan_linkages[glycan_id] = linkage
    
    # Sort glycans by size (largest first)
    glycan_ids = list(optimized_data['glycans'].keys())
    glycan_sizes = []
    for glycan_id in glycan_ids:
        if glycan_id in glycan_linkages:
            size = len(optimized_data['glycans'][glycan_id]['atoms'])
            glycan_sizes.append((glycan_id, size))
    
    glycan_sizes.sort(key=lambda x: x[1], reverse=True)
    sorted_glycan_ids = [x[0] for x in glycan_sizes]
    
    # Store all optimization results
    all_results = {}
    unconverged_glycans = []
    
    # Multiple optimization cycles
    for cycle in range(max_cycles):
        print_and_save(f"\n{'='*40}")
        print_and_save(f"CYCLE {cycle + 1}/{max_cycles}")
        print_and_save(f"{'='*40}")
        
        cycle_results = {}
        cycle_unconverged = []
        
        for i, glycan_id in enumerate(sorted_glycan_ids):
            if glycan_id not in glycan_linkages:
                if verbose:
                    print_and_save(f"Skipping {glycan_id}: No linkage information")
                continue
            
            linkage = glycan_linkages[glycan_id]
            
            # Optimize this glycan with current cycle number using single-axis rotation
            optimized_data, result = optimize_glycan_single_axis(
                optimized_data, glycan_id, linkage, theta_step,
                n_steps, radius_angstrom, use_coulomb, n_workers, 
                verbose, cycle=cycle+1,
                save_before_after=save_before_after,
                output_dir=glycans_output_dir
            )
            
            cycle_results[glycan_id] = result
            
            # Check convergence (energy improvement > 1 kJ/mol)
            if result['improvement'] < 1.0:
                cycle_unconverged.append(glycan_id)
                if verbose:
                    print_and_save(f"  Glycan {glycan_id}: Did not converge (improvement: {result['improvement']:.2f} kJ/mol)")
            else:
                if verbose:
                    print_and_save(f"  Glycan {glycan_id}: Converged (improvement: {result['improvement']:.2f} kJ/mol)")
        
        # Store results for this cycle
        all_results[f'cycle_{cycle+1}'] = cycle_results
        
        # Update unconverged list
        unconverged_glycans = cycle_unconverged
        
        # Check if all converged
        if not unconverged_glycans:
            print_and_save(f"\n✓ All glycans converged in cycle {cycle + 1}")
            break
        
        # Update sorted list for next cycle
        sorted_glycan_ids = unconverged_glycans + [
            g for g in sorted_glycan_ids if g not in unconverged_glycans
        ]
        
        # For next cycles, use updated coordinates (already stored in optimized_data)
        # No need to explicitly update as we're using the same data object
        
        if cycle < max_cycles - 1:
            if verbose:
                print_and_save(f"\nNext cycle parameters:")
                print_and_save(f"  theta_step: {theta_step}°")
                print_and_save(f"  MCMC steps: {n_steps}")
    
    # Save final complex PDB
    if output_pdb:
        save_pdb_file(optimized_data, output_pdb)
        print_and_save(f"Saved optimized complex to: {output_pdb}")
    
    # Add optimization metadata
    optimized_data['optimization'] = {
        'parameters': {
            'theta_step': theta_step,
            'mcmc_steps': n_steps,
            'max_cycles': max_cycles,
            'radius_angstrom': radius_angstrom,
            'use_coulomb': use_coulomb,
            'n_workers': n_workers,
            'save_before_after': save_before_after,
            'glycans_output_dir': glycans_output_dir
        },
        'results': all_results,
        'summary': {
            'total_glycans': len(glycan_linkages),
            'cycles_completed': min(cycle + 1, max_cycles),
            'unconverged_glycans': unconverged_glycans,
            'unconverged_count': len(unconverged_glycans)
        },
        'timestamp': datetime.now().isoformat()
    }
    
    # Print final summary
    print_and_save("\n" + "=" * 80)
    print_and_save("OPTIMIZATION SUMMARY")
    print_and_save("=" * 80)
    
    total_glycans = len(glycan_linkages)
    converged_count = total_glycans - len(unconverged_glycans)
    
    print_and_save(f"Total glycans processed: {total_glycans}")
    print_and_save(f"Glycans converged: {converged_count} ({converged_count/total_glycans*100 if total_glycans > 0 else 0:.1f}%)")
    
    # Calculate collision statistics
    total_collisions_before = 0
    total_collisions_after = 0
    collision_reductions = []
    
    for cycle_key, cycle_results in all_results.items():
        for glycan_id, result in cycle_results.items():
            total_collisions_before += result.get('collisions_before', 0)
            total_collisions_after += result.get('collisions_after', 0)
            collision_reductions.append(result.get('collision_reduction', 0))
    
    print_and_save(f"\nCOLLISION STATISTICS:")
    print_and_save(f"  Total collisions before optimization: {total_collisions_before}")
    print_and_save(f"  Total collisions after optimization: {total_collisions_after}")
    
    if total_collisions_before > 0:
        collision_reduction_pct = (total_collisions_before - total_collisions_after) / total_collisions_before * 100
        print_and_save(f"  Collision reduction: {total_collisions_before - total_collisions_after} ({collision_reduction_pct:.1f}%)")
    
    if collision_reductions:
        print_and_save(f"  Average collision reduction per glycan: {np.mean(collision_reductions):.1f}")
    
    if unconverged_glycans:
        print_and_save(f"\nUnconverged glycans ({len(unconverged_glycans)}):")
        for glycan_id in unconverged_glycans:
            # Get best improvement from all cycles
            best_improvement = 0
            for cycle_key, cycle_results in all_results.items():
                if glycan_id in cycle_results:
                    improvement = cycle_results[glycan_id].get('improvement', 0)
                    if improvement > best_improvement:
                        best_improvement = improvement
            
            print_and_save(f"  {glycan_id}: Best improvement = {best_improvement:.2f} kJ/mol")
    
    # Print average improvements
    all_improvements = []
    for cycle_key, cycle_results in all_results.items():
        for glycan_id, result in cycle_results.items():
            all_improvements.append(result.get('improvement', 0))
    
    if all_improvements:
        print_and_save(f"\nENERGY IMPROVEMENT STATISTICS:")
        print_and_save(f"  Average energy improvement: {np.mean(all_improvements):.2f} ± {np.std(all_improvements):.2f} kJ/mol")
        print_and_save(f"  Maximum improvement: {np.max(all_improvements):.2f} kJ/mol")
        print_and_save(f"  Minimum improvement: {np.min(all_improvements):.2f} kJ/mol")
    
    # Print summary of coordinate changes
    print_and_save("\nCOORDINATE CHANGE SUMMARY:")
    coordinate_changes_applied = 0
    for cycle_key, cycle_results in all_results.items():
        for glycan_id, result in cycle_results.items():
            if result.get('coordinate_change_applied', False):
                coordinate_changes_applied += 1
    
    print_and_save(f"Coordinate changes applied: {coordinate_changes_applied}/{total_glycans} glycans")
    
    # Print file output summary
    print_and_save("\nOUTPUT FILES:")
    if output_pdb:
        print_and_save(f"  Initial complex PDB: {output_pdb.replace('.pdb', '_initial.pdb')}")
        print_and_save(f"  Optimized complex PDB: {output_pdb}")
    
    if save_before_after and glycans_output_dir:
        print_and_save(f"  Individual glycan PDBs: {glycans_output_dir}/")
        print_and_save(f"    Format: {{glycan_id}}_{{before/after}}_{{cycle}}.pdb")
    
    print_and_save("=" * 80)
    
    return optimized_data

def main():
    """Main function for glycan optimization."""
    parser = argparse.ArgumentParser(
        description="Optimize glycan orientations using single-axis rotation followed by MCMC refinement"
    )
    
    parser.add_argument("--input_json", required=True,
                       help="Path to input JSON file with CHARMM36 parameters")
    parser.add_argument("--output_json", required=True,
                       help="Path to output JSON file with optimized coordinates")
    parser.add_argument("--output_pdb", required=False, default=None,
                       help="Path to output PDB file with updated coordinates")
    parser.add_argument("--save_individual_glycans", action="store_true", default=False,
                       help="Save individual glycans as separate PDB files")
    parser.add_argument("--glycans_output_dir", required=False, default="PDB_CARB_ONLY",
                       help="Directory to save individual glycan PDB files (default: PDB_CARB_ONLY)")
    parser.add_argument("--theta_step", type=int, default=30,
                       help="Step size for theta grid search in degrees (single axis, default: 30)")
    parser.add_argument("--n_steps", type=int, default=10000,
                       help="Number of MCMC refinement steps (default: 10000)")
    parser.add_argument("--max_cycles", type=int, default=5,
                       help="Maximum optimization cycles (default: 5)")
    parser.add_argument("--radius", type=float, default=80.0,
                       help="Radius for local energy calculation in Angstroms (default: 80)")
    parser.add_argument("--use_coulomb", type=str, default="false",
                       choices=["true", "false", "yes", "no"],
                       help="Include Coulomb energy in calculation (default: false)")
    parser.add_argument("--n_workers", type=int, default=4,
                       help="Number of parallel workers (default: 4)")
    parser.add_argument("--verbose", action="store_true", default=True,
                       help="Print progress information (default: True)")
    parser.add_argument("--report_file", type=str, default="report.txt",
                       help="Path to report file where all output will be saved (default: report.txt)")
    parser.add_argument("--save_before_after", action="store_true", default=True,
                       help="Save PDB files before and after optimization for each glycan")
    
    args = parser.parse_args()
    
    # Convert use_coulomb to boolean
    use_coulomb = args.use_coulomb.lower() in ["true", "yes"]
    
    # Set default output_pdb if not provided
    if args.output_pdb is None:
        args.output_pdb = args.output_json.replace('.json', '.pdb')
    
    # Validate arguments
    if args.n_workers > cpu_count():
        print_and_save(f"Warning: Requested {args.n_workers} workers but only {cpu_count()} CPUs available")
        args.n_workers = max(1, cpu_count() - 1)
    
    print_and_save("Loading input data...")
    with open(args.input_json, 'r') as f:
        data = json.load(f)
    
    print_and_save(f"Loaded data from {args.input_json}")
    print_and_save(f"Protein atoms: {len(data['protein'])}")
    print_and_save(f"Glycans: {len(data['glycans'])}")
    print_and_save(f"Linkages: {len(data['linkages'])}")
    
    # Print glycan IDs for reference
    if args.save_individual_glycans or args.verbose:
        print_and_save("\nGlycan IDs found in input:")
        for glycan_id in sorted(data['glycans'].keys()):
            num_atoms = len(data['glycans'][glycan_id]['atoms'])
            # Count residues
            residues = set()
            for atom in data['glycans'][glycan_id]['atoms']:
                residues.add((atom['residue_name'], atom['residue_number']))
            print_and_save(f"  {glycan_id}: {num_atoms} atoms, {len(residues)} residues")
    
    # Optimize glycan orientations using single-axis rotation
    optimized_data = optimize_all_glycans_single_axis(
        data,
        theta_step=args.theta_step,
        n_steps=args.n_steps,
        max_cycles=args.max_cycles,
        radius_angstrom=args.radius,
        use_coulomb=use_coulomb,
        n_workers=args.n_workers,
        verbose=args.verbose,
        report_file=args.report_file,
        save_before_after=args.save_before_after,
        output_pdb=args.output_pdb,
        glycans_output_dir=args.glycans_output_dir
    )
    
    # Save optimized data
    print_and_save(f"\nSaving optimized data to {args.output_json}...")
    with open(args.output_json, 'w') as f:
        json.dump(optimized_data, f, indent=2, cls=NumpyEncoder)
    
    # Save PDB file with updated coordinates (already done in optimize_all_glycans_single_axis)
    # save_pdb_file(optimized_data, args.output_pdb)
    
    # Save individual glycans if requested
    if args.save_individual_glycans:
        save_individual_glycans_pdb(optimized_data, args.glycans_output_dir)
    
    print_and_save("Optimization complete!")
    
    # Create summary file
    summary_file = args.output_json.replace('.json', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("GLYCAN OPTIMIZATION SUMMARY - SINGLE AXIS ROTATION\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("PARAMETERS:\n")
        f.write(f"  Input file: {args.input_json}\n")
        f.write(f"  Output JSON: {args.output_json}\n")
        f.write(f"  Output PDB: {args.output_pdb}\n")
        f.write(f"  Save individual glycans: {args.save_individual_glycans}\n")
        if args.save_individual_glycans:
            f.write(f"  Glycans output directory: {args.glycans_output_dir}\n")
        f.write(f"  Save before/after PDBs: {args.save_before_after}\n")
        f.write(f"  Report file: {args.report_file}\n")
        f.write(f"  Theta step: {args.theta_step}°\n")
        f.write(f"  MCMC steps: {args.n_steps}\n")
        f.write(f"  Max cycles: {args.max_cycles}\n")
        f.write(f"  Search radius: {args.radius} Å\n")
        f.write(f"  Include Coulomb: {use_coulomb}\n")
        f.write(f"  Parallel workers: {args.n_workers}\n\n")
        
        f.write("RESULTS:\n")
        summary = optimized_data.get('optimization', {}).get('summary', {})
        
        total = summary.get('total_glycans', 0)
        unconverged = summary.get('unconverged_count', 0)
        converged = total - unconverged
        
        f.write(f"  Total glycans: {total}\n")
        f.write(f"  Converged: {converged} ({converged/total*100 if total > 0 else 0:.1f}%)\n")
        f.write(f"  Cycles completed: {summary.get('cycles_completed', 0)}\n\n")
        
        # Collision statistics
        total_collisions_before = 0
        total_collisions_after = 0
        
        all_results = optimized_data.get('optimization', {}).get('results', {})
        for cycle_key, cycle_results in all_results.items():
            for glycan_id, result in cycle_results.items():
                total_collisions_before += result.get('collisions_before', 0)
                total_collisions_after += result.get('collisions_after', 0)
        
        f.write(f"COLLISION STATISTICS:\n")
        f.write(f"  Total collisions before optimization: {total_collisions_before}\n")
        f.write(f"  Total collisions after optimization: {total_collisions_after}\n")
        
        if total_collisions_before > 0:
            collision_reduction_pct = (total_collisions_before - total_collisions_after) / total_collisions_before * 100
            f.write(f"  Collision reduction: {total_collisions_before - total_collisions_after} ({collision_reduction_pct:.1f}%)\n\n")
        
        if unconverged > 0:
            f.write(f"  Unconverged glycans ({unconverged}):\n")
            for glycan_id in summary.get('unconverged_glycans', []):
                f.write(f"    - {glycan_id}\n")
        
        # Calculate average improvement
        all_improvements = []
        for cycle_key, cycle_results in all_results.items():
            for glycan_id, result in cycle_results.items():
                improvement = result.get('improvement', 0)
                all_improvements.append(improvement)
        
        if all_improvements:
            f.write(f"\n  Energy improvements:\n")
            f.write(f"    Average: {np.mean(all_improvements):.2f} kJ/mol\n")
            f.write(f"    Std Dev: {np.std(all_improvements):.2f} kJ/mol\n")
            f.write(f"    Maximum: {np.max(all_improvements):.2f} kJ/mol\n")
            f.write(f"    Minimum: {np.min(all_improvements):.2f} kJ/mol\n")
        
        # PDB file info
        f.write(f"\nOUTPUT FILES:\n")
        f.write(f"  JSON file: {args.output_json}\n")
        if args.output_pdb:
            f.write(f"  Initial complex PDB: {args.output_pdb.replace('.pdb', '_initial.pdb')}\n")
            f.write(f"  Optimized complex PDB: {args.output_pdb}\n")
        
        if args.save_individual_glycans:
            f.write(f"  Individual glycan PDBs: {args.glycans_output_dir}/ directory\n")
        
        if args.save_before_after:
            f.write(f"  Before/after PDBs for each glycan: {args.glycans_output_dir}/ directory\n")
            f.write(f"    Format: {{glycan_id}}_{{before/after}}_{{cycle}}.pdb\n")
        
        f.write(f"  Summary file: {summary_file}\n")
        f.write(f"  Report file: {args.report_file}\n")
    
    print_and_save(f"Summary saved to {summary_file}")
    print_and_save(f"Report saved to {args.report_file}")
    
    if args.output_pdb:
        print_and_save(f"Initial complex PDB saved to: {args.output_pdb.replace('.pdb', '_initial.pdb')}")
        print_and_save(f"Optimized complex PDB saved to: {args.output_pdb}")
    
    if args.save_individual_glycans:
        print_and_save(f"Individual glycan PDB files saved to {args.glycans_output_dir}/ directory")
    
    if args.save_before_after:
        print_and_save(f"Before/after PDB files for each glycan saved to {args.glycans_output_dir}/ directory")

if __name__ == "__main__":
    main()
