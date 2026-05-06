#!/usr/bin/env python3
"""
Script to process glycosylated PDB files (distance-based, robust)

- DOES NOT change sugar labels
- Detects N-glycosylation via ASN ND2 to NAG/NDG C1 distance
- Extracts ONLY the glycan block connected to each specific glycosylation site
"""

import os
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict
import argparse
import MDAnalysis as mda

# =========================
# CONFIG
# =========================

WATER_RESNAMES = {
    "HOH", "WAT", "TIP3", "TIP3P", "SOL"
}

ION_RESNAMES = {
    "NA", "CL", "K", "CA", "MG", "ZN", "MN",
    "FE", "CU", "CO", "NI", "CD"
}

# N-glycosidic bond cutoff distance (Å)
N_GLY_CUTOFF = 2.0  # Å

# Specific N-glycan residue names (first residue attached to ASN)
N_GLYCAN_FIRST_RESIDUES = {
    "NDG",  # N-acetylglucosamine (GROMACS naming)
    "NAG",  # N-acetylglucosamine (alternative naming)
    "NGN",  # N-glycan
    "GLCNAC",  # N-acetylglucosamine (full name)
}

# =========================
# DETECT GLYCAN RESNAMES AUTOMATICALLY
# =========================

def detect_glycan_resnames(pdb_file):
    """Automatically detect glycan residue names using MDAnalysis"""
    u = mda.Universe(pdb_file)
    
    glycan_resnames = set()
    
    for res in u.residues:
        if res.atoms.select_atoms("protein").n_atoms > 0:
            continue
        
        resname = res.resname.strip()
        
        if resname in WATER_RESNAMES:
            continue
        
        if resname in ION_RESNAMES:
            continue
        
        if len(res.atoms) > 1:
            glycan_resnames.add(resname)
    
    print(f"Detected {len(glycan_resnames)} glycan residue types:")
    for name in sorted(glycan_resnames):
        print(f"  - {name}")
    
    return glycan_resnames


# =========================
# BASIC UTILITIES
# =========================

def extract_coordinates_from_line(line):
    try:
        return (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        )
    except:
        return None


def distance(a, b):
    if a is None or b is None:
        return float('inf')
    return np.linalg.norm(np.array(a) - np.array(b))


# =========================
# REMOVE HYDROGENS
# =========================

def remove_hydrogens(input_pdb, output_pdb):
    """Remove hydrogen atoms from PDB file"""
    with open(input_pdb) as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                atom = line[12:16].strip()
                element = line[76:78].strip()
                if element == "H" or atom.startswith("H"):
                    continue
            fout.write(line)
    print(f"  Removed hydrogens -> {output_pdb}")


# =========================
# PARSE PDB AND BUILD CONNECTIVITY
# =========================

def parse_pdb_and_build_connectivity(pdb_file):
    """
    Parse PDB and build residue connectivity based on sequential numbering
    and chain breaks
    """
    all_lines = []
    with open(pdb_file) as f:
        all_lines = f.readlines()
    
    # Group lines by chain and residue
    residues = {}  # (chain, resnum) -> list of lines
    for line in all_lines:
        if line.startswith(("ATOM", "HETATM")):
            chain = line[21]
            resnum = int(line[22:26])
            key = (chain, resnum)
            if key not in residues:
                residues[key] = []
            residues[key].append(line)
    
    return residues, all_lines


def find_connected_glycan_block(residues, start_chain, start_resnum, glycan_resnames):
    """
    Find all residues connected to a starting glycan residue
    by traversing sequentially until a break is found.
    A break occurs when:
      1. Next residue number is not consecutive (gap > 1)
      2. Next residue is from a different chain
      3. Next residue is not a glycan (protein, water, ion)
    """
    connected = []
    current_resnum = start_resnum
    
    # Sort residues in this chain by residue number
    chain_residues = sorted([r for r in residues.keys() if r[0] == start_chain], key=lambda x: x[1])
    
    # Find the index of the starting residue
    start_idx = None
    for i, (chain, resnum) in enumerate(chain_residues):
        if chain == start_chain and resnum == start_resnum:
            start_idx = i
            break
    
    if start_idx is None:
        return []
    
    # Traverse forward from start residue
    for i in range(start_idx, len(chain_residues)):
        chain, resnum = chain_residues[i]
        
        # Check if this is still the same chain
        if chain != start_chain:
            break
        
        # For residues beyond the first, check if consecutive
        if i > start_idx:
            prev_resnum = chain_residues[i-1][1]
            if resnum - prev_resnum > 1:
                # Break in numbering - new glycan block
                break
        
        # Get residue name
        res_lines = residues[(chain, resnum)]
        if not res_lines:
            continue
        
        resname = res_lines[0][17:20].strip()
        
        # Check if this is a glycan residue
        if resname in glycan_resnames:
            connected.append((chain, resnum, res_lines))
        else:
            # Hit a non-glycan residue (protein, water, ion) - stop
            break
    
    return connected


# =========================
# DETECT GLYCAN ROOTS
# =========================

def detect_glycan_roots(residues, all_lines, glycan_resnames):
    """
    Identify N-glycosylation sites by detecting ASN ND2 to NAG/NDG C1 distance
    """
    # First, collect all ASN ND2 atoms
    asn_nd2_atoms = []
    for (chain, resnum), lines in residues.items():
        for line in lines:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if resname == "ASN" and atomname == "ND2":
                    coords = extract_coordinates_from_line(line)
                    if coords:
                        asn_nd2_atoms.append({
                            "coords": coords,
                            "chain": chain,
                            "resnum": resnum,
                            "line": line
                        })
    
    print(f"  Found {len(asn_nd2_atoms)} ASN ND2 atoms")
    
    # Collect all C1 atoms from glycan residues
    glycan_c1_atoms = []
    for (chain, resnum), lines in residues.items():
        for line in lines:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if resname in glycan_resnames and atomname == "C1":
                    coords = extract_coordinates_from_line(line)
                    if coords:
                        glycan_c1_atoms.append({
                            "coords": coords,
                            "chain": chain,
                            "resnum": resnum,
                            "resname": resname,
                            "line": line
                        })
    
    print(f"  Found {len(glycan_c1_atoms)} glycan C1 atoms")
    
    # Detect N-glycosidic bonds
    roots = []
    for glycan in glycan_c1_atoms:
        for protein in asn_nd2_atoms:
            dist = distance(glycan["coords"], protein["coords"])
            if dist <= N_GLY_CUTOFF:
                roots.append({
                    "glycan_chain": glycan["chain"],
                    "glycan_resnum": glycan["resnum"],
                    "glycan_resname": glycan["resname"],
                    "protein_chain": protein["chain"],
                    "protein_resnum": protein["resnum"],
                    "distance": round(dist, 3)
                })
                print(f"  N-glycosylation: ASN{protein['resnum']}{protein['chain']} -> "
                      f"{glycan['resname']}{glycan['resnum']}{glycan['chain']} (dist: {dist:.3f} Å)")
    
    return roots


# =========================
# EXTRACT GLYCAN BLOCKS
# =========================

def extract_glycan_blocks(residues, roots, glycan_resnames):
    """
    Extract each glycan block separately (not all glycans together)
    """
    glycans = {}
    
    for i, root in enumerate(roots, 1):
        glycan_id = f"GLYCAN_{i}"
        
        # Find the connected block starting from this root
        connected = find_connected_glycan_block(
            residues, 
            root["glycan_chain"], 
            root["glycan_resnum"], 
            glycan_resnames
        )
        
        if not connected:
            continue
        
        # Extract all lines for this glycan block
        lines = []
        for chain, resnum, res_lines in connected:
            lines.extend(res_lines)
        
        # Build residue sequence
        residue_info = {}
        for line in lines:
            resn = line[17:20].strip()
            resi = int(line[22:26])
            if resi not in residue_info:
                residue_info[resi] = resn
        
        seq = [residue_info[r] for r in sorted(residue_info)]
        
        glycans[glycan_id] = {
            "lines": lines,
            "chain": root["glycan_chain"],
            "residue_numbers": [r for _, r, _ in connected],
            "start_residue": min([r for _, r, _ in connected]),
            "end_residue": max([r for _, r, _ in connected]),
            "residues": seq,
            "residue_sequence": "_".join(f"{r}{i+1}" for i, r in enumerate(seq)),
            "simple_residue_sequence": "_".join(seq),
            "unit_number": len(seq),
            "protein_resnum": root["protein_resnum"],
            "protein_chain": root["protein_chain"],
            "distance": root["distance"]
        }
        
        print(f"  Glycan {glycan_id}: {len(seq)} residues -> {'-'.join(seq)}")
    
    return glycans


# =========================
# SAVE GLYCANS
# =========================

def save_glycan_pdbs(glycans, output_dir):
    pdb_dir = Path(output_dir) / "PDB"
    pdb_dir.mkdir(parents=True, exist_ok=True)

    for gid, data in glycans.items():
        out = pdb_dir / f"{gid}.pdb"
        with open(out, "w") as f:
            for l in data["lines"]:
                f.write(l)
            f.write("TER\n")

    print(f"  Saved {len(glycans)} glycan PDBs to {pdb_dir}")
    return pdb_dir


# =========================
# SAVE TSV
# =========================

def save_linkages_tsv(glycans, output_dir):
    out = Path(output_dir) / "topol_carb.tsv"
    with open(out, "w") as f:
        f.write(
            "site_protein_residue\tglycan_binding\tsequence_poly\tunit_number\t"
            "protein_residue_number\tprotein_chain\tprotein_atom\tglycan_atom\t"
            "glycan_resname\tdistance\tlinking_type\n"
        )
        for gid, data in glycans.items():
            f.write(f"ASN{data['protein_resnum']}\t")
            f.write(f"{gid}\t")
            f.write(f"{data['residue_sequence']}\t")
            f.write(f"{data['unit_number']}\t")
            f.write(f"{data['protein_resnum']}\t")
            f.write(f"{data['protein_chain']}\t")
            f.write(f"ND2\t")
            f.write(f"C1\t")
            f.write(f"{data['residues'][0] if data['residues'] else 'UNK'}\t")
            f.write(f"{data['distance']}\t")
            f.write(f"N-linked\n")
    
    print(f"  Saved {len(glycans)} linkages to {out}")


# =========================
# MAIN
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="Extract glycan coordinates from glycosylated PDB structure",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input_pdb", required=True,
                        help="Input PDB file (glycosylated protein)")
    parser.add_argument("--output_noH", required=True,
                        help="Output PDB file without hydrogens")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for glycans and TSV file")
    args = parser.parse_args()
    
    # Convert to absolute paths
    input_pdb = Path(args.input_pdb).resolve()
    output_noH = Path(args.output_noH).resolve()
    output_dir = Path(args.output_dir).resolve()
    
    if not input_pdb.exists():
        print(f"Error: Input PDB not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("EXTRACTING GLYCANS FROM GLYCOSYLATED PROTEIN")
    print("=" * 70)
    print(f"Input:  {input_pdb}")
    print(f"Output: {output_dir}")
    print()
    
    # Detect glycan residue names
    print("-" * 70)
    print("Step 1: Detecting glycan residue names...")
    print("-" * 70)
    glycan_resnames = detect_glycan_resnames(str(input_pdb))
    
    # Remove hydrogens
    print("\n" + "-" * 70)
    print("Step 2: Removing hydrogens...")
    print("-" * 70)
    remove_hydrogens(str(input_pdb), str(output_noH))
    
    # Parse PDB and build connectivity
    print("\n" + "-" * 70)
    print("Step 3: Parsing PDB and building connectivity...")
    print("-" * 70)
    residues, all_lines = parse_pdb_and_build_connectivity(str(input_pdb))
    print(f"  Found {len(residues)} unique residues")
    
    # Detect glycosylation roots
    print("\n" + "-" * 70)
    print("Step 4: Detecting N-glycosylation sites...")
    print("-" * 70)
    roots = detect_glycan_roots(residues, all_lines, glycan_resnames)
    
    if not roots:
        print("  No N-glycosylation sites found!")
        sys.exit(0)
    
    # Extract glycan blocks
    print("\n" + "-" * 70)
    print("Step 5: Extracting individual glycan blocks...")
    print("-" * 70)
    glycans = extract_glycan_blocks(residues, roots, glycan_resnames)
    print(f"  Extracted {len(glycans)} independent glycan block(s)")
    
    # Save glycan PDBs
    print("\n" + "-" * 70)
    print("Step 6: Saving glycan PDB files...")
    print("-" * 70)
    save_glycan_pdbs(glycans, str(output_dir))
    
    # Save TSV
    print("\n" + "-" * 70)
    print("Step 7: Saving linkage TSV...")
    print("-" * 70)
    save_linkages_tsv(glycans, str(output_dir))
    
    print("\n" + "=" * 70)
    print("EXTRACTION COMPLETED SUCCESSFULLY!")
    print("=" * 70)
    print(f"  No-H PDB:      {output_noH}")
    print(f"  Glycans PDBs:  {output_dir}/PDB/")
    print(f"  Linkages TSV:  {output_dir}/topol_carb.tsv")
    print("=" * 70)


if __name__ == "__main__":
    main()
