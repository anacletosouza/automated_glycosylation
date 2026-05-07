#!/usr/bin/env python3
import os
import argparse

def extract_protein_atoms(pdb_path, remove_hydrogens=True):
    """
    Keep only protein ATOM records.
    Optionally remove hydrogens.
    All HETATM entries from the input PDB are discarded.
    """
    protein_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                if remove_hydrogens:
                    atom_name = line[12:16].strip()
                    if atom_name.startswith("H"):
                        continue
                protein_lines.append(line)
    return protein_lines


def extract_carbohydrate_atoms(pdb_path, remove_hydrogens=True):
    """
    Extract carbohydrate atoms from a processed glycan PDB.
    Only HETATM entries are kept.
    Hydrogens are optionally removed.
    """
    carb_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):
                if remove_hydrogens:
                    atom_name = line[12:16].strip()
                    if atom_name.startswith("H"):
                        continue
                carb_lines.append(line)
    return carb_lines


def main(protein_pdb, carbs_dir, output_pdb,
         keep_hydrogens_protein, keep_hydrogens_carb):

    final_pdb = []

    remove_hydrogens_prot = not keep_hydrogens_protein
    remove_hydrogens_carb = not keep_hydrogens_carb

    print(f"Reading protein PDB: {protein_pdb}")
    protein_atoms = extract_protein_atoms(
        protein_pdb, remove_hydrogens=remove_hydrogens_prot
    )
    final_pdb.extend(protein_atoms)
    final_pdb.append("TER\n")

    print(f"Scanning carbohydrate directories in: {carbs_dir}")
    carb_dirs = sorted(
        d for d in os.listdir(carbs_dir)
        if os.path.isdir(os.path.join(carbs_dir, d))
    )

    for cdir in carb_dirs:
        modified_pdb = os.path.join(carbs_dir, cdir, f"{cdir}_modified.pdb")
        if os.path.isfile(modified_pdb):
            print(f"Adding glycan from: {modified_pdb}")
            carb_atoms = extract_carbohydrate_atoms(
                modified_pdb, remove_hydrogens=remove_hydrogens_carb
            )
            if carb_atoms:
                final_pdb.extend(carb_atoms)
                final_pdb.append("TER\n")
        else:
            print(f"Warning: {modified_pdb} not found, skipping.")

    final_pdb.append("END\n")

    with open(output_pdb, "w") as f:
        f.writelines(final_pdb)

    print(f"Final glycoprotein written to: {output_pdb}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build a glycoprotein PDB by combining a protein and processed glycans."
    )
    parser.add_argument("--protein", required=True, help="Protein PDB file (input)")
    parser.add_argument("--carbs_dir", required=True, help="Directory containing glycan folders")
    parser.add_argument("--output", required=True, help="Final glycoprotein PDB file")

    parser.add_argument(
        "--keep_hydrogens_protein",
        action="store_true",
        help="Keep hydrogen atoms in protein (default: removed)"
    )
    parser.add_argument(
        "--keep_hydrogens_carb",
        action="store_true",
        help="Keep hydrogen atoms in carbohydrates (default: removed)"
    )

    args = parser.parse_args()

    main(
        args.protein,
        args.carbs_dir,
        args.output,
        args.keep_hydrogens_protein,
        args.keep_hydrogens_carb
    )

