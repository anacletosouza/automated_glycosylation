import argparse
import re
import os
import shutil
import subprocess

def find_residuetypes():
    """Automatically locate residuetypes.dat based on GROMACS executable."""
    try:
        gmx_path = subprocess.check_output(["which", "gmx"], text=True).strip()
    except subprocess.CalledProcessError:
        print("GROMACS executable 'gmx' not found in PATH.")
        return None

    bin_dir = os.path.dirname(gmx_path)
    root_dir = os.path.dirname(bin_dir)  # .../gromacs-xxxx
    residuetypes_path = os.path.join(root_dir, "share", "gromacs", "top", "residuetypes.dat")

    if os.path.isfile(residuetypes_path):
        print(f"Found residuetypes.dat at: {residuetypes_path}")
        return residuetypes_path
    else:
        print(f"residuetypes.dat not found at expected location: {residuetypes_path}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Add all residue names from a carb_unique.rtp file to residuetypes.dat."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input carb_unique RPT file (carb_unique_redundance_removed.rtp)"
    )
    parser.add_argument(
        "--output", "-o", required=False,
        help="residuetypes.dat file to update (optional; will auto-detect if not provided)"
    )
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    # If no output is provided, try to find residuetypes.dat automatically
    if not output_file:
        output_file = find_residuetypes()
        if not output_file:
            raise FileNotFoundError("Could not locate residuetypes.dat automatically. Please provide --output.")

    # Create a backup of residuetypes.dat
    backup_file = output_file + ".backup"
    if os.path.isfile(output_file):
        shutil.copy2(output_file, backup_file)
        print(f"Backup created: {backup_file}")

    # Regex to find residue names: [ RESNAME ]
    resname_pattern = re.compile(r"^\[\s*(\w+)\s*\]", re.MULTILINE)

    # Read residue names from the RPT file
    with open(input_file, "r") as f:
        content = f.read()
    resnames = resname_pattern.findall(content)

    # Read existing residue names from residuetypes.dat
    existing_resnames = set()
    if os.path.isfile(output_file):
        with open(output_file, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    existing_resnames.add(line.split()[0])

    # Filter out residues already in the file
    new_resnames = [r for r in resnames if r not in existing_resnames]

    # Append new residues to the file with type "Carbohydrate"
    if new_resnames:
        with open(output_file, "a") as f:
            for r in new_resnames:
                f.write(f"{r}\tCarbohydrate\n")

    print(f"Added {len(new_resnames)} new residue(s) to {output_file}.")

if __name__ == "__main__":
    main()

