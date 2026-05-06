#!/usr/bin/env python3
"""
Script to process glycosylated PDB files (distance-based, robust)

- Detects N-glycosylation via ASN ND2/ND1? to NAG/NDG C1 distance
- Detects O-glycosylation via SER OG or THR OG1 to A2G C1 distance
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

# Protein residue names for detection
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "HIE", "HID", "HIP"
}

# N-glycosidic bond cutoff distance (Å)
N_GLY_CUTOFF = 2.0  # Å

# O-glycosidic bond cutoff distance (Å)
O_GLY_CUTOFF = 2.0  # Å

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


# =========================
# DETECT GLYCAN ROOTS (N and O-linked)
# =========================

def detect_glycan_roots(residues, all_lines, glycan_resnames, residue_start_carn):
    """
    Identify N-glycosylation and O-glycosylation sites by distance.
    
    residue_start_carn format: "O:A2G,N:NDG" (default)
    """
    # Parse residue_start_carn
    n_initiation_residues = set()
    o_initiation_residues = set()
    
    for item in residue_start_carn.split(','):
        if ':' in item:
            linkage_type, resname = item.split(':')
            if linkage_type.upper() == 'N':
                n_initiation_residues.add(resname.strip().upper())
            elif linkage_type.upper() == 'O':
                o_initiation_residues.add(resname.strip().upper())
    
    print(f"  N-linked initiation residues: {n_initiation_residues}")
    print(f"  O-linked initiation residues: {o_initiation_residues}")
    
    # Collect ASN atoms for N-glycosylation
    asn_atoms = []
    for (chain, resnum), lines in residues.items():
        for line in lines:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if resname == "ASN" and atomname in ["ND2", "ND1"]:
                    coords = extract_coordinates_from_line(line)
                    if coords:
                        asn_atoms.append({
                            "coords": coords,
                            "chain": chain,
                            "resnum": resnum,
                            "resname": resname,
                            "atomname": atomname,
                            "line": line
                        })
    
    print(f"  Found {len(asn_atoms)} ASN ND2/ND1 atoms")
    
    # Collect SER/THR atoms for O-glycosylation
    ser_thr_atoms = []
    for (chain, resnum), lines in residues.items():
        for line in lines:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if resname in ["SER", "THR"]:
                    if (resname == "SER" and atomname == "OG") or \
                       (resname == "THR" and atomname == "OG1"):
                        coords = extract_coordinates_from_line(line)
                        if coords:
                            ser_thr_atoms.append({
                                "coords": coords,
                                "chain": chain,
                                "resnum": resnum,
                                "resname": resname,
                                "atomname": atomname,
                                "line": line
                            })
    
    print(f"  Found {len(ser_thr_atoms)} SER OG / THR OG1 atoms")
    
    # Collect initiation atoms from glycan residues (C1 atoms)
    n_initiation_atoms = []  # For N-linked
    o_initiation_atoms = []  # For O-linked
    
    for (chain, resnum), lines in residues.items():
        for line in lines:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if atomname == "C1":
                    coords = extract_coordinates_from_line(line)
                    if coords:
                        atom_info = {
                            "coords": coords,
                            "chain": chain,
                            "resnum": resnum,
                            "resname": resname,
                            "line": line
                        }
                        if resname in n_initiation_residues:
                            n_initiation_atoms.append(atom_info)
                        if resname in o_initiation_residues:
                            o_initiation_atoms.append(atom_info)
    
    print(f"  Found {len(n_initiation_atoms)} N-linked initiation C1 atoms")
    print(f"  Found {len(o_initiation_atoms)} O-linked initiation C1 atoms")
    
    # Detect N-glycosidic bonds
    n_roots = []
    for glycan in n_initiation_atoms:
        for protein in asn_atoms:
            dist = distance(glycan["coords"], protein["coords"])
            if dist <= N_GLY_CUTOFF:
                n_roots.append({
                    "type": "N-linked",
                    "glycan_chain": glycan["chain"],
                    "glycan_resnum": glycan["resnum"],
                    "glycan_resname": glycan["resname"],
                    "protein_chain": protein["chain"],
                    "protein_resnum": protein["resnum"],
                    "protein_resname": protein["resname"],
                    "protein_atom": protein["atomname"],
                    "glycan_atom": "C1",
                    "distance": round(dist, 3)
                })
                print(f"  N-glycosylation: {protein['resname']}{protein['resnum']}{protein['chain']} "
                      f"({protein['atomname']}) -> {glycan['resname']}{glycan['resnum']}{glycan['chain']} "
                      f"(dist: {dist:.3f} Å)")
    
    # Detect O-glycosidic bonds
    o_roots = []
    for glycan in o_initiation_atoms:
        for protein in ser_thr_atoms:
            dist = distance(glycan["coords"], protein["coords"])
            if dist <= O_GLY_CUTOFF:
                o_roots.append({
                    "type": "O-linked",
                    "glycan_chain": glycan["chain"],
                    "glycan_resnum": glycan["resnum"],
                    "glycan_resname": glycan["resname"],
                    "protein_chain": protein["chain"],
                    "protein_resnum": protein["resnum"],
                    "protein_resname": protein["resname"],
                    "protein_atom": protein["atomname"],
                    "glycan_atom": "C1",
                    "distance": round(dist, 3)
                })
                print(f"  O-glycosylation: {protein['resname']}{protein['resnum']}{protein['chain']} "
                      f"({protein['atomname']}) -> {glycan['resname']}{glycan['resnum']}{glycan['chain']} "
                      f"(dist: {dist:.3f} Å)")
    
    return n_roots, o_roots


# =========================
# SAVE TSV WITH GLYCOSYLATION SITES
# =========================

def save_linkages_tsv(n_roots, o_roots, output_dir):
    out = Path(output_dir) / "glycosylation_sites.tsv"
    
    with open(out, "w") as f:
        f.write(
            "site_protein_residue\tglycan_residue\t"
            "protein_residue_number\tprotein_chain\tprotein_atom\tglycan_atom\t"
            "glycan_resname\tglycan_chain\tglycan_resnum\tdistance\tlinking_type\n"
        )
        
        # Write N-linked sites
        for root in n_roots:
            site = f"{root['protein_resname']}{root['protein_resnum']}"
            glycan_site = f"{root['glycan_resname']}{root['glycan_resnum']}"
            f.write(f"{site}\t")
            f.write(f"{glycan_site}\t")
            f.write(f"{root['protein_resnum']}\t")
            f.write(f"{root['protein_chain']}\t")
            f.write(f"{root['protein_atom']}\t")
            f.write(f"{root['glycan_atom']}\t")
            f.write(f"{root['glycan_resname']}\t")
            f.write(f"{root['glycan_chain']}\t")
            f.write(f"{root['glycan_resnum']}\t")
            f.write(f"{root['distance']}\t")
            f.write(f"{root['type']}\n")
        
        # Write O-linked sites
        for root in o_roots:
            site = f"{root['protein_resname']}{root['protein_resnum']}"
            glycan_site = f"{root['glycan_resname']}{root['glycan_resnum']}"
            f.write(f"{site}\t")
            f.write(f"{glycan_site}\t")
            f.write(f"{root['protein_resnum']}\t")
            f.write(f"{root['protein_chain']}\t")
            f.write(f"{root['protein_atom']}\t")
            f.write(f"{root['glycan_atom']}\t")
            f.write(f"{root['glycan_resname']}\t")
            f.write(f"{root['glycan_chain']}\t")
            f.write(f"{root['glycan_resnum']}\t")
            f.write(f"{root['distance']}\t")
            f.write(f"{root['type']}\n")
    
    print(f"  Saved {len(n_roots) + len(o_roots)} glycosylation sites to {out}")


# =========================
# MAIN
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="Detect N and O-glycosylation sites in PDB structure",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input_pdb", required=True,
                        help="Input PDB file (glycosylated protein)")
    parser.add_argument("--output_noH", required=True,
                        help="Output PDB file without hydrogens")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for TSV file")
    parser.add_argument("--residue_start_carn", type=str, default="O:A2G,N:NDG",
                        help="Residue names for initiating N and O glycans (default: O:A2G,N:NDG)")
    parser.add_argument("--n_gly_cutoff", type=float, default=2.0,
                        help="Distance cutoff for N-glycosidic bonds (default: 2.0 Å)")
    parser.add_argument("--o_gly_cutoff", type=float, default=2.0,
                        help="Distance cutoff for O-glycosidic bonds (default: 2.0 Å)")
    
    args = parser.parse_args()
    
    # Set global cutoffs
    global N_GLY_CUTOFF, O_GLY_CUTOFF
    N_GLY_CUTOFF = args.n_gly_cutoff
    O_GLY_CUTOFF = args.o_gly_cutoff
    
    # Convert to absolute paths
    input_pdb = Path(args.input_pdb).resolve()
    output_noH = Path(args.output_noH).resolve()
    output_dir = Path(args.output_dir).resolve()
    
    if not input_pdb.exists():
        print(f"Error: Input PDB not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("DETECTING GLYCOSYLATION SITES IN PROTEIN")
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
    
    # Detect glycosylation roots (N and O-linked)
    print("\n" + "-" * 70)
    print("Step 4: Detecting N and O-glycosylation sites...")
    print("-" * 70)
    n_roots, o_roots = detect_glycan_roots(residues, all_lines, glycan_resnames, args.residue_start_carn)
    
    all_roots = n_roots + o_roots
    
    if not all_roots:
        print("  No glycosylation sites found!")
        sys.exit(0)
    
    print(f"\n  Total glycosylation sites found: {len(all_roots)}")
    print(f"    N-linked: {len(n_roots)}")
    print(f"    O-linked: {len(o_roots)}")
    
    # Save TSV
    print("\n" + "-" * 70)
    print("Step 5: Saving glycosylation sites TSV...")
    print("-" * 70)
    save_linkages_tsv(n_roots, o_roots, str(output_dir))
    
    print("\n" + "=" * 70)
    print("DETECTION COMPLETED SUCCESSFULLY!")
    print("=" * 70)
    print(f"  No-H PDB:      {output_noH}")
    print(f"  Sites TSV:     {output_dir}/glycosylation_sites.tsv")
    print("=" * 70)


if __name__ == "__main__":
    main()
