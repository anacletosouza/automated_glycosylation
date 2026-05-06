import glycosylator as gl
import matplotlib.pyplot as plt
import pandas as pd
import math
import argparse

# -----------------------------
# Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Glycosylate a protein based on a glycosylation table.")
parser.add_argument("--input_tsv_glycosylator", required=True, help="Path to the glycosylation table (TSV format)")
parser.add_argument("--input_pdb_protein", required=True, help="Path to the protonated protein (PDB format)")
parser.add_argument("--protein_residue_start", type=int, default=10, help="Residue number where the protein starts")
parser.add_argument("--output", default="glycosylated_protein.pdb", help="Output PDB file path")
args = parser.parse_args()

PROTEIN_START_RESIDUE = args.protein_residue_start

# -----------------------------
# Helper: Find residue by number
# -----------------------------
def find_residue(protein, chain_id, residue_num):
    """Find a specific residue in the protein using residue number."""
    for chain in protein.chains:
        if chain.id == chain_id:
            for res in chain.residues:
                if isinstance(res.id, tuple):
                    res_number = res.id[1]
                elif isinstance(res.id, (int, float)):
                    res_number = int(res.id)
                else:
                    try:
                        res_number = int(str(res.id))
                    except:
                        continue
                if res_number == residue_num:
                    return res
    return None

# -----------------------------
# Helper: Rotate glycan to reduce clashes
# -----------------------------
def rotate_glycan_to_avoid_clashes(glycan, protein, pivot_atom_name='ND2', axis='y', step_deg=10, max_deg=360):
    """
    Try rotating the glycan around pivot_atom_name to minimize clashes with protein.
    axis: 'x', 'y', 'z'
    step_deg: rotation step in degrees
    max_deg: total rotation range
    """
    pivot = glycan.get_atom(pivot_atom_name)
    if pivot is None:
        return  # cannot rotate if pivot not found

    for angle in range(0, max_deg, step_deg):
        glycan.rotate_around_atom(pivot_atom_name, angle, axis=axis)
        clashes = glycan.find_clashes_with(protein)
        if not clashes:
            break

# -----------------------------
# Load glycosylation table
# -----------------------------
glyco_table = pd.read_csv(args.input_tsv_glycosylator, sep="\t")

# Load protein
protein = gl.protein(args.input_pdb_protein)
protein.reindex()

print(f"Structure loaded: {len(protein.chains)} chains")

# -----------------------------
# Calculate chain offsets
# -----------------------------
chain_offsets = {}
for chain in protein.chains:
    residues = list(chain.residues)
    res_numbers = []
    for res in residues:
        if isinstance(res.id, tuple):
            res_numbers.append(res.id[1])
        elif isinstance(res.id, (int, float)):
            res_numbers.append(int(res.id))
        else:
            try:
                res_numbers.append(int(str(res.id)))
            except:
                pass
    if res_numbers:
        chain_offsets[chain.id] = min(res_numbers) - 1
        print(f"Chain {chain.id}: residues {len(residues)}, numbering {min(res_numbers)}-{max(res_numbers)}, offset={chain_offsets[chain.id]}")

# -----------------------------
# Separate N- and O-glycosylations
# -----------------------------
n_glyco_data = glyco_table[glyco_table['site'].str.startswith('N')].copy()
o_glyco_data = glyco_table[glyco_table['site'].str.startswith(('S','T'))].copy()
residues_to_glycosylate = []
ignored_residues = []

# -----------------------------
# Process N-glycosylations
# -----------------------------
print("\nProcessing N-glycosylations...")
for _, row in n_glyco_data.iterrows():
    site = row['site']
    iupac_seq = row['iupac_glycosylator'].strip()
    table_res_num = int(row['residue_number'])
    chain_id = row['protein_chain']
    if chain_id not in chain_offsets:
        ignored_residues.append((f"{chain_id}:{table_res_num}", "Chain not found"))
        continue
    real_res_num = table_res_num + chain_offsets[chain_id] - (PROTEIN_START_RESIDUE - 1)
    res = find_residue(protein, chain_id, real_res_num)
    if res and res.name == 'ASN' and 'ND2' in [a.name for a in res.atoms]:
        residues_to_glycosylate.append((res, iupac_seq, 'N', f"{chain_id}:{table_res_num}"))
    else:
        ignored_residues.append((f"{chain_id}:{table_res_num}", f"Residue not suitable"))

# -----------------------------
# Process O-glycosylations
# -----------------------------
print("\nProcessing O-glycosylations...")
for _, row in o_glyco_data.iterrows():
    site = row['site']
    iupac_seq = row['iupac_glycosylator'].strip()
    table_res_num = int(row['residue_number'])
    chain_id = row['protein_chain']
    if chain_id not in chain_offsets:
        ignored_residues.append((f"{chain_id}:{table_res_num}", "Chain not found"))
        continue
    real_res_num = table_res_num + chain_offsets[chain_id] - (PROTEIN_START_RESIDUE - 1)
    res = find_residue(protein, chain_id, real_res_num)
    if res and ((res.name=='SER' and 'OG' in [a.name for a in res.atoms]) or (res.name=='THR' and 'OG1' in [a.name for a in res.atoms])):
        residues_to_glycosylate.append((res, iupac_seq, 'O', f"{chain_id}:{table_res_num}"))
    else:
        ignored_residues.append((f"{chain_id}:{table_res_num}", "Residue not suitable"))

print(f"\n{len(residues_to_glycosylate)} residues ready for glycosylation.")

# -----------------------------
# Apply standard bonds
# -----------------------------
all_residues = [res for res, _, _, _ in residues_to_glycosylate]
print(f"\nApplying standard bonds for {len(all_residues)} residues...")
try:
    protein.apply_standard_bonds_for(*all_residues)
    print("Standard bonds applied.")
except Exception as e:
    print(f"Error applying standard bonds: {e}")

# -----------------------------
# Glycosylation with rotation to avoid clashes
# -----------------------------
print("\nStarting glycosylation with clash avoidance...")
successful = 0
failed = 0

for res, iupac_seq, glyco_type, location in residues_to_glycosylate:
    try:
        glycan = gl.glycan(iupac_seq)
        glycan.infer_glycan_tree()
        glycan.remove_hydrogens()
        glycan.apply_standard_bonds()
        glycan.add_hydrogens()

        # 1. Optimize glycan geometry
        glycan.optimize()

        # 2. Rotate to avoid clashes with protein
        pivot_atom = 'ND2' if glyco_type=='N' else 'OG' if res.name=='SER' else 'OG1'
        rotate_glycan_to_avoid_clashes(glycan, protein, pivot_atom_name=pivot_atom)

        # 3. Attach to protein
        protein.glycosylate(glycan, residues=[res])

        print(f"✓ Glycosylated {location} ({res.name}) with {glyco_type}-glycan")
        successful += 1

    except Exception as e:
        print(f"✗ Failed glycosylation at {location}: {e}")
        failed += 1

print(f"\nGlycosylation completed: {successful} successful, {failed} failed.")

# -----------------------------
# SNFG Visualization
# -----------------------------
try:
    protein.snfg()
    plt.show()
except Exception as e:
    print(f"Error visualizing SNFG: {e}")

# -----------------------------
# Save final PDB
# -----------------------------
output_file = args.output
try:
    protein.to_pdb(output_file)
    print(f"Structure saved as {output_file}")
except Exception as e:
    print(f"Error saving PDB: {e}")

# -----------------------------
# Summary
# -----------------------------
print(f"\nTotal residues in table: {len(glyco_table)}")
print(f"Residues successfully glycosylated: {successful}")
print(f"Residues ignored/failed: {len(ignored_residues) + failed}")

