import json

def reconstruct_pdb_from_json(json_file, output_pdb):
    """
    Reconstruct a PDB file from JSON data with proper TER and residue numbering.
    
    Parameters:
    -----------
    json_file : str
        Path to input JSON file
    output_pdb : str
        Path to output PDB file
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    with open(output_pdb, 'w') as f:
        atom_counter = 1  # Optional: can renumber atoms sequentially if needed

        # -----------------------
        # 1. Write protein atoms
        # -----------------------
        for atom in data['protein']:
            line = f"{atom['record']:6s}{atom['atom_number']:5d} {atom['atom_name']:4s}" \
                   f"{atom['alt_loc']:1s}{atom['residue_name']:3s} {atom['chain_id']:1s}" \
                   f"{atom['residue_number']:4d}{atom['icode']:1s}   " \
                   f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}" \
                   f"{atom['occupancy']:6.2f}{atom['temp_factor']:6.2f}" \
                   f"          {atom['element']:2s}{atom['charge']:2s}\n"
            f.write(line)
            atom_counter += 1

        # Write TER after protein block
        f.write("TER\n")

        # -----------------------
        # 2. Write glycans
        # -----------------------
        for glycan_id, glycan_data in data['glycans'].items():
            residue_mapping = {}
            new_residue_number = 1  # restart numbering for each glycan

            # Determine mapping from old residue numbers to new ones
            for old_res_num in sorted(glycan_data['residue_numbers']):
                residue_mapping[old_res_num] = new_residue_number
                new_residue_number += 1

            # Write atoms for this glycan
            for atom in glycan_data['atoms']:
                new_res_num = residue_mapping[atom['residue_number']]
                line = f"{atom['record']:6s}{atom['atom_number']:5d} {atom['atom_name']:4s}" \
                       f"{atom['alt_loc']:1s}{atom['residue_name']:3s} {atom['chain_id']:1s}" \
                       f"{new_res_num:4d}{atom['icode']:1s}   " \
                       f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}" \
                       f"{atom['occupancy']:6.2f}{atom['temp_factor']:6.2f}" \
                       f"          {atom['element']:2s}{atom['charge']:2s}\n"
                f.write(line)
                atom_counter += 1

            # Write TER at the end of this glycan block
            f.write("TER\n")

        # -----------------------
        # 3. Write remaining other lines (CONECT, END)
        # -----------------------
        for line in data['other_lines']:
            if line.startswith('CONECT') or line.startswith('END'):
                f.write(line + '\n')

        # Ensure file ends with END if not already present
        if not data['other_lines'] or not data['other_lines'][-1].startswith('END'):
            f.write("END\n")

    print(f"PDB reconstruction complete: {output_pdb}")


# -----------------------
# Example usage
# -----------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Reconstruct PDB from JSON")
    parser.add_argument("--json_file", required=True, help="Input JSON file")
    parser.add_argument("--output_pdb", required=True, help="Output PDB file")

    args = parser.parse_args()

    reconstruct_pdb_from_json(args.json_file, args.output_pdb)

