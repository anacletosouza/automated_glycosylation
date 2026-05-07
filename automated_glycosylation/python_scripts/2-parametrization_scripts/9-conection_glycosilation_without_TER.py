# -*- coding: utf-8 -*-
"""
Script to remove TER lines from a PDB file, add CONECT lines from a reference PDB,
and preserve the original END lines.

Usage:
    python 9-conection_glycosilation.py \
        --glycosylated input_pdb \
        --conect connect_pdb \
        --output output_pdb
"""

import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Remove TER lines from a glycosylated PDB file and add CONECT lines."
    )
    parser.add_argument(
        "--glycosylated",
        required=True,
        help="Path to the glycosylated PDB file (input)"
    )
    parser.add_argument(
        "--conect",
        required=True,
        help="Path to the PDB file containing CONECT lines (reference)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output PDB file"
    )

    args = parser.parse_args()

    # Step 1: Read original glycosylated PDB and remove TER lines
    with open(args.glycosylated, "r") as f:
        pdb_lines = f.readlines()

    pdb_lines_no_ter = [
        line for line in pdb_lines if line.strip() != "TER"
    ]

    # Step 2: Read CONECT lines from reference PDB
    with open(args.conect, "r") as f:
        conect_lines = [
            line for line in f if line.startswith("CONECT")
        ]

    # Step 3: Preserve END lines from the original PDB
    end_lines = []
    while pdb_lines_no_ter and pdb_lines_no_ter[-1].strip() == "END":
        end_lines.insert(0, pdb_lines_no_ter.pop())

    # If no END line was found, add one
    if not end_lines:
        end_lines = ["END\n"]

    # Step 4: Write the new PDB file with CONECT before END
    with open(args.output, "w") as f:
        f.writelines(pdb_lines_no_ter)
        f.writelines(conect_lines)
        f.writelines(end_lines)

    print(f"File created successfully: {args.output}")


if __name__ == "__main__":
    main()

