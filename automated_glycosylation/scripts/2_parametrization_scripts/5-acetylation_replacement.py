#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

if len(sys.argv) != 3:
    print("Usage: python3 5-acetylation_replacement.py <input.pdb> <output.pdb>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Replacements: (atom_name, residue_name) -> new_atom_name

replacements = {
    ("C1",  "NDG"): "C1",
    ("C2",  "NDG"): "C2",
    ("C3",  "NDG"): "C3",
    ("C4",  "NDG"): "C4",
    ("C5",  "NDG"): "C5",
    ("C6",  "NDG"): "C6",
    ("C7",  "NDG"): "C",
    ("C8",  "NDG"): "CT",
    ("O5",  "NDG"): "O5",
    ("O3",  "NDG"): "O3",
    ("O4",  "NDG"): "O4",
    ("O6",  "NDG"): "O6",
    ("O7",  "NDG"): "O",
    ("N2",  "NDG"): "N",
    ("HN2",  "NDG"): "HN",
    ("H1",  "NDG"): "H1",
    ("H2",  "NDG"): "H2",
    ("H3",  "NDG"): "H3",
    ("H4",  "NDG"): "H4",
    ("H5",  "NDG"): "H5",
    ("H61", "NDG"): "H61",
    ("H62", "NDG"): "H62",
    ("H81", "NDG"): "HT1",
    ("H82", "NDG"): "HT2",
    ("H83", "NDG"): "HT3",
    ("HO3", "NDG"): "HO3",
    ("C1",  "NAG"): "C1",
    ("C2",  "NAG"): "C2",
    ("C3",  "NAG"): "C3",
    ("C4",  "NAG"): "C4",
    ("C5",  "NAG"): "C5",
    ("C6",  "NAG"): "C6",
    ("C7",  "NAG"): "C",
    ("C8",  "NAG"): "CT",
    ("O5",  "NAG"): "O5",
    ("O3",  "NAG"): "O3",
    ("O4",  "NAG"): "O4",
    ("O6",  "NAG"): "O6",
    ("O7",  "NAG"): "O",
    ("N2",  "NAG"): "N",
    ("HN2",  "NAG"): "HN",
    ("H1",  "NAG"): "H1",
    ("H2",  "NAG"): "H2",
    ("H3",  "NAG"): "H3",
    ("H4",  "NAG"): "H4",
    ("H5",  "NAG"): "H5",
    ("H61", "NAG"): "H61",
    ("H62", "NAG"): "H62",
    ("H81", "NAG"): "HT1",
    ("H82", "NAG"): "HT2",
    ("H83", "NAG"): "HT3",
    ("HO3", "NAG"): "HO3",
    ("C1",   "SIA"): "C1",
    ("C2",   "SIA"): "C2",
    ("C3",   "SIA"): "C3",
    ("C4",   "SIA"): "C4",
    ("C5",   "SIA"): "C5",
    ("C6",   "SIA"): "C6",
    ("C7",   "SIA"): "C7",
    ("C8",   "SIA"): "C8",
    ("C9",   "SIA"): "C9",
    ("C10",  "SIA"): "C",
    ("C11",  "SIA"): "CT",
    ("N5",   "SIA"): "N",
    ("O1A",  "SIA"): "O11",
    ("O1B",  "SIA"): "O12",
    ("O4",   "SIA"): "O4",
    ("O6",   "SIA"): "O6",
    ("O7",   "SIA"): "O7",
    ("O8",   "SIA"): "O8",
    ("O9",   "SIA"): "O9",
    ("O10",  "SIA"): "O",
    ("H31",  "SIA"): "H31",
    ("H32",  "SIA"): "H32",
    ("H4",   "SIA"): "H4",
    ("H5",   "SIA"): "H5",
    ("H6",   "SIA"): "H6",
    ("H7",   "SIA"): "H7",
    ("H8",   "SIA"): "H8",
    ("H91",  "SIA"): "H91",
    ("H92",  "SIA"): "H92",
    ("H111", "SIA"): "HT1",
    ("H112", "SIA"): "HT2",
    ("H113", "SIA"): "HT3",
    ("HN5",  "SIA"): "HN",
    ("HO1B", "SIA"): "HO2",
    ("HO4",  "SIA"): "HO4",
    ("HO7",  "SIA"): "HO7",
    ("HO8",  "SIA"): "HO8",
    ("HO9",  "SIA"): "HO9",
    ("O5",  "A2G"): "O5",
    ("C1",  "A2G"): "C1",
    ("C2",  "A2G"): "C2",
    ("N2",  "A2G"): "N",
    ("C3",  "A2G"): "C3",
    ("O3",  "A2G"): "O3",
    ("C4",  "A2G"): "C4",
    ("O4",  "A2G"): "O4",
    ("C5",  "A2G"): "C5",
    ("C6",  "A2G"): "C6",
    ("O6",  "A2G"): "O6",
    ("C7",  "A2G"): "C",
    ("O7",  "A2G"): "O",
    ("C8",  "A2G"): "CT",
    ("H1",  "A2G"): "H1",
    ("H2",  "A2G"): "H2",
    ("HN2", "A2G"): "HN",
    ("H3",  "A2G"): "H3",
    ("H4",  "A2G"): "H4",
    ("HO4", "A2G"): "HO4",
    ("H5",  "A2G"): "H5",
    ("H61", "A2G"): "H61",
    ("H62", "A2G"): "H62",
    ("HO6", "A2G"): "HO6",
    ("H81", "A2G"): "HT1",
    ("H82", "A2G"): "HT2",
    ("H83", "A2G"): "HT3"
}

def pdb_atom_format(atom_name):
    """Format atom name according to PDB conventions (columns 13-16)."""
    if len(atom_name) == 1:
        return f" {atom_name}  "
    elif len(atom_name) == 2:
        return f" {atom_name} "
    elif len(atom_name) == 3:
        return f" {atom_name}"
    elif len(atom_name) == 4:
        return atom_name
    else:
        raise ValueError(f"Invalid atom name length: {atom_name}")

with open(input_file, "r") as f:
    lines = f.readlines()

modified_lines = []
for line in lines:
    if line.startswith(("HETATM", "ATOM")):
        atom_name = line[12:16].strip()
        res_name  = line[17:20].strip()
        key = (atom_name, res_name)
        if key in replacements:
            new_atom = pdb_atom_format(replacements[key])
            line = line[:12] + new_atom + line[16:]
    modified_lines.append(line)

with open(output_file, "w") as f:
    f.writelines(modified_lines)

print(f"Replacements completed. Modified file saved as '{output_file}'.")

