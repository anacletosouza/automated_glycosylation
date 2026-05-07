import re
import glycosylator as gl
import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse

# ============================================================
# PARSE COMMAND-LINE ARGUMENTS
# ============================================================
parser = argparse.ArgumentParser(description="Convert CHARMM sequences to glycosylator-compatible IUPAC and generate glycans")
parser.add_argument("--input_tsv", required=True, help="Path to input TSV file")
parser.add_argument("--output_tsv", required=True, help="Path to output TSV file")
parser.add_argument("--output_dir", required=True, help="Path to output directory")
args = parser.parse_args()

INPUT_TSV = args.input_tsv
OUTPUT_TSV = args.output_tsv
OUTPUT_DIR = args.output_dir

# ============================================================
# LOAD INPUT
# ============================================================
df = pd.read_csv(INPUT_TSV, sep="\t")

if "sequence" not in df.columns:
    raise ValueError("Column 'sequence' not found in input TSV")

# ============================================================
# EXTRACT PROTEIN CHAIN AND RESIDUE NUMBER (PROA-123)
# ============================================================
def extract_protein_info(seq: str):
    """
    Extract protein chain and residue number from PRO{chain}-{resnum}
    Example: PROA-123 → chain=A, residue_number=123
    """
    if pd.isna(seq):
        return pd.Series([None, None])

    match = re.search(r"PRO([A-Z])-(\d+)", seq)
    if match:
        return pd.Series([match.group(1), int(match.group(2))])

    return pd.Series([None, None])

df[["protein_chain", "residue_number"]] = df["sequence"].apply(extract_protein_info)

# ============================================================
# CHARMM → GLYCAM-like
# ============================================================
def charmm_to_glycam(seq: str) -> str:
    if pd.isna(seq):
        return seq

    # Remove PRO{chain}-{residue}
    seq = re.sub(r"PRO[A-Z]-\d+", "", seq)

    replacements = {
        "aDMan": "DManpa",
        "bDMan": "DManpb",
        "aDGlcNAc": "DGlcpNAca",
        "bDGlcNAc": "DGlcpNAcb",
        "aDGal": "DGalpa",
        "bDGal": "DGalpb",
        "aDGalNAc": "DGalpNAca",
        "bDGalNAc": "DGalpNAcb",
        "aLFuc": "LFucpa",
        "bLFuc": "LFucpb",
        "aDNeu5Ac": "DNeu5Aca",
        "bDNeu5Ac": "DNeu5Acb",
    }

    for k, v in replacements.items():
        seq = seq.replace(k, v)

    seq = (
        seq.replace("(1→6)", "1-6")
           .replace("(1→4)", "1-4")
           .replace("(1→3)", "1-3")
           .replace("(1→2)", "1-2")
           .replace("(2→6)", "2-6")
           .replace("(2→3)", "2-3")
           .replace("(1→)", "1-OH")
    )

    return seq.replace(" ", "").strip()

# ============================================================
# GLYCAM-like → IUPAC condensed
# ============================================================
def glycam_to_iupac(s: str) -> str:
    if pd.isna(s):
        return s

    s = re.sub(r"-OH$", "", s)
    s = re.sub(r"([A-Za-z]+)[pf]", r"\1", s)
    s = re.sub(r"([ab]\d-\d)", r"(\1)", s)

    return s

# ============================================================
# IUPAC → GLYCOSYLATOR-COMPATIBLE IUPAC
# ============================================================
def to_glycosylator_iupac(s: str) -> str:
    if pd.isna(s):
        return s

    # Remove D/L attached to residues
    s = re.sub(
        r"([ab]-)?[DL](Man|GlcNAc|GalNAc|Gal|Fuc|Neu5Ac|Neu5Gc)",
        r"\1\2",
        s,
    )

    # Remove any remaining D/L tokens
    s = re.sub(r"\b[DL]\b", "", s)

    # Remove reducing-end linkage
    s = re.sub(r"(a|b)1$", "", s)

    # Remove trailing numbers
    s = re.sub(r"\d+$", "", s)

    # Fix malformed names
    s = s.replace("GalaNAc", "GalNAc")
    s = s.replace("GlcaNAc", "GlcNAc")

    return s

# ============================================================
# APPLY PIPELINE
# ============================================================
df["glycam"] = df["sequence"].apply(charmm_to_glycam)
df["iupac_condensed"] = df["glycam"].apply(glycam_to_iupac)
df["iupac_glycosylator"] = df["iupac_condensed"].apply(to_glycosylator_iupac)

# ============================================================
# SAVE OUTPUT TSV
# ============================================================
df.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("Conversion completed successfully")
print(f"Output saved to: {OUTPUT_TSV}")

# ============================================================
# OUTPUT DIRECTORIES
# ============================================================
BASE_DIR = OUTPUT_DIR
SVG_DIR = os.path.join(BASE_DIR, "SVG")
PNG_DIR = os.path.join(BASE_DIR, "PNG")

os.makedirs(SVG_DIR, exist_ok=True)
os.makedirs(PNG_DIR, exist_ok=True)

# ============================================================
# SAFETY CHECK
# ============================================================
required_cols = {"name", "iupac_glycosylator"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# ============================================================
# DRAW AND SAVE GLYCANS
# ============================================================
for _, row in df.iterrows():

    name = str(row["name"])
    iupac = row["iupac_glycosylator"]

    if pd.isna(iupac) or not iupac.strip():
        continue

    try:
        g = gl.glycan(iupac)
        ax = g.draw2d()
        fig = ax.figure
        fig.tight_layout()

        fig.savefig(os.path.join(SVG_DIR, f"{name}.svg"), bbox_inches="tight")
        fig.savefig(os.path.join(PNG_DIR, f"{name}.png"), dpi=300, bbox_inches="tight")

        plt.close(fig)
        print(f"Saved: {name}")

    except Exception as e:
        print(f"Failed for {name}: {e}")

# ============================================================
# GENERATE PDB FILES
# ============================================================
PDB_DIR = os.path.join(OUTPUT_DIR, "PDB_files_glycosylator")
os.makedirs(PDB_DIR, exist_ok=True)

for _, row in df.iterrows():

    name = str(row["name"])
    iupac = row["iupac_glycosylator"]

    if pd.isna(iupac) or not iupac.strip():
        continue

    try:
        g = gl.glycan(iupac)
        g.optimize()
        g.to_pdb(os.path.join(PDB_DIR, f"{name}.pdb"))
        print(f"Saved PDB: {name}")

    except Exception as e:
        print(f"Failed for {name}: {e}")

# ============================================================
# PREPROCESS PDBs
# ============================================================
PDB_PREPROCESS_DIR = os.path.join(OUTPUT_DIR, "PDB_preprocessed")
os.makedirs(PDB_PREPROCESS_DIR, exist_ok=True)

def preprocess_glycan(name: str, iupac: str, output_dir: str):

    if pd.isna(iupac) or not iupac.strip():
        return

    try:
        g = gl.glycan(iupac)
        g.remove_hydrogens()
        g.apply_standard_bonds()
        g.add_hydrogens()

        for a1, a2 in g.find_clashes():
            g.adjust_bond_length(a1, a2, 1.4)

        g.optimize()
        g.to_pdb(os.path.join(output_dir, f"{name}.pdb"))

        print(f"Preprocessed PDB saved: {name}")

    except Exception as e:
        print(f"Failed for {name}: {e}")

for _, row in df.iterrows():
    preprocess_glycan(row["name"], row["iupac_glycosylator"], PDB_PREPROCESS_DIR)

