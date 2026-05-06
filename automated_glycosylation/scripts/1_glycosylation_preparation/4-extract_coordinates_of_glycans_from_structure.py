#!/usr/bin/env python3
"""
Script to process glycosylated PDB files (distance-based, robust)

- DOES NOT change sugar labels
- Detects N-glycosylation via ASN ND2 distance
- Does NOT rely on C1 atom
"""

import os
import numpy as np
from pathlib import Path
from collections import defaultdict
import argparse

# =========================
# CONFIG
# =========================

GLYCAN_RESNAMES = {
    "AMAN", "BMAN", "BGLC", "BGAL", "AGAL", "AFUC"
}

N_GLY_CUTOFF = 2.0  # Å

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
    return np.linalg.norm(np.array(a) - np.array(b))


# =========================
# REMOVE HYDROGENS
# =========================

def remove_hydrogens(input_pdb, output_pdb):
    with open(input_pdb) as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                atom = line[12:16].strip()
                element = line[76:78].strip()
                if element == "H" or atom.startswith("H"):
                    continue
            fout.write(line)


# =========================
# PARSE PDB
# =========================

def parse_pdb(pdb):
    protein, glycans = [], []
    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                protein.append(line)
            elif line.startswith("HETATM"):
                glycans.append(line)
            elif line.startswith("TER"):
                protein.append(line)
    return protein, glycans


# =========================
# DETECT GLYCAN ROOTS
# =========================

def detect_glycan_roots(protein_lines, glycan_lines):
    """
    Identify sugar residues linked to ASN ND2 by distance
    """
    asn_nd2_atoms = []
    sugar_atoms = []

    for line in protein_lines:
        if line.startswith("ATOM"):
            if line[17:20].strip() == "ASN" and line[12:16].strip() == "ND2":
                asn_nd2_atoms.append({
                    "coords": extract_coordinates_from_line(line),
                    "chain": line[21],
                    "resnum": int(line[22:26])
                })

    for line in glycan_lines:
        if line[17:20].strip() in GLYCAN_RESNAMES:
            sugar_atoms.append({
                "coords": extract_coordinates_from_line(line),
                "chain": line[21],
                "resnum": int(line[22:26])
            })

    roots = set()
    for s in sugar_atoms:
        for a in asn_nd2_atoms:
            if distance(s["coords"], a["coords"]) <= N_GLY_CUTOFF:
                roots.add((s["chain"], s["resnum"]))

    return roots


# =========================
# EXTRACT GLYCANS
# =========================

def extract_glycans(glycan_lines, roots):
    glycans = {}
    grouped = defaultdict(lambda: defaultdict(list))

    for line in glycan_lines:
        chain = line[21]
        resnum = int(line[22:26])
        grouped[chain][resnum].append(line)

    gid = 1
    for chain, residues in grouped.items():
        sorted_res = sorted(residues.keys())

        for root_chain, root_res in roots:
            if root_chain != chain:
                continue

            connected = [r for r in sorted_res if r >= root_res]
            lines = []
            for r in connected:
                lines.extend(residues[r])

            glycan_id = f"{chain}_{gid}"
            gid += 1

            glycans[glycan_id] = {
                "lines": lines,
                "chain": chain,
                "residue_numbers": connected,
                "start_residue": min(connected),
                "end_residue": max(connected)
            }

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

        residue_info = {}
        for l in data["lines"]:
            resn = l[17:20].strip()
            resi = int(l[22:26])
            residue_info.setdefault(resi, resn)

        seq = [residue_info[r] for r in sorted(residue_info)]
        data["residues"] = seq
        data["residue_sequence"] = "_".join(f"{r}{i+1}" for i, r in enumerate(seq))
        data["simple_residue_sequence"] = "_".join(seq)
        data["unit_number"] = len(seq)

    return pdb_dir


# =========================
# FIND LINKAGES
# =========================

def find_glycan_linkages(protein_lines, glycans):
    protein_atoms = []

    for l in protein_lines:
        if l.startswith("ATOM"):
            protein_atoms.append({
                "coords": extract_coordinates_from_line(l),
                "res": l[17:20].strip(),
                "resnum": int(l[22:26]),
                "chain": l[21],
                "atom": l[12:16].strip()
            })

    linkages = []

    for gid, g in glycans.items():
        first_line = g["lines"][0]
        gcoords = extract_coordinates_from_line(first_line)

        best = None
        dmin = 999

        for p in protein_atoms:
            if p["res"] == "ASN" and p["atom"] == "ND2":
                d = distance(gcoords, p["coords"])
                if d < dmin:
                    dmin = d
                    best = p

        if best and dmin < 5.0:
            linkages.append({
                "site_protein_residue": f"ASN{best['resnum']}",
                "glycan_binding": gid,
                "sequence_poly": g["residue_sequence"],
                "unit_number": g["unit_number"],
                "protein_residue_number": best["resnum"],
                "protein_chain": best["chain"],
                "protein_atom": "ND2",
                "glycan_atom": g["residues"][0],
                "distance": round(dmin, 3),
                "linking_type": "N-linked"
            })

    return linkages


# =========================
# SAVE TSV
# =========================

def save_linkages_tsv(linkages, output_dir):
    out = Path(output_dir) / "topol_carb.tsv"
    with open(out, "w") as f:
        f.write(
            "site_protein_residue\tglycan_binding\tsequence_poly\tunit_number\t"
            "protein_residue_number\tprotein_chain\tprotein_atom\tglycan_atom\t"
            "distance\tlinking_type\n"
        )
        for l in linkages:
            f.write("\t".join(str(l[k]) for k in l) + "\n")


# =========================
# MAIN
# =========================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb", required=True)
    parser.add_argument("--output_noH", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    remove_hydrogens(args.input_pdb, args.output_noH)
    protein, glycan_lines = parse_pdb(args.input_pdb)
    roots = detect_glycan_roots(protein, glycan_lines)
    glycans = extract_glycans(glycan_lines, roots)
    save_glycan_pdbs(glycans, args.output_dir)
    linkages = find_glycan_linkages(protein, glycans)
    save_linkages_tsv(linkages, args.output_dir)


if __name__ == "__main__":
    main()

