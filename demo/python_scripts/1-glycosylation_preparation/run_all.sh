#!/bin/bash
# -----------------------------------------------------------------------------
# run_all.sh
# Shell script to run the full glycosylation workflow
# Each step corresponds to a separate Python script execution
# -----------------------------------------------------------------------------

set -e  # Stop script if any command fails

# -------------------------------
# Define paths
# -------------------------------
INPUT_PATH="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/DELTA/INPUT_DIR"
OUTPUT_PATH="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/DELTA/1-GLYCOPROTEIN_PREPARATION"
SCRIPTS_PYTHON="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/DELTA/python_scripts/1-glycosylation_preparation"

# -------------------------------
# Create output directories
# -------------------------------
echo "Creating output directories..."

mkdir -p "$OUTPUT_PATH/TSV"
mkdir -p "$OUTPUT_PATH/PDB_PROTEIN_GLYCOSYLATED"
mkdir -p "$OUTPUT_PATH/EXTRACTED_CARBOHYDRATES"
mkdir -p "$OUTPUT_PATH/TO_TOP"

# -------------------------------
# Step 0: Correct Caselino table
# -------------------------------
INPUT_TSV="$INPUT_PATH/caselino_2020_table_Delta.tsv"
OUTPUT_TSV="$OUTPUT_PATH/TSV/caselino_2020_tables_corrected.tsv"
OUTPUT_RESULTS="$OUTPUT_PATH/EXTRACTED_CARBOHYDRATES"

echo "Step 0: Correcting Caselino table..."
if [[ ! -f "$INPUT_TSV" ]]; then
    echo "ERROR: Input TSV not found: $INPUT_TSV"
    exit 1
fi

python3 "$SCRIPTS_PYTHON"/0-correcting_caselino_table_for_variants.py \
    --input "$INPUT_TSV" \
    --output "$OUTPUT_TSV"

# ---------------------------------------
# Step 1: Convert to IUPAC notation
# ---------------------------------------
OUTPUT_GLYCOSYLATOR="$OUTPUT_PATH/TSV/caselino_2020_tables_glycosylator.tsv"

echo "Step 1: Converting to IUPAC notation..."
python3 "$SCRIPTS_PYTHON"/1-iupac_converted.py \
    --input_tsv "$OUTPUT_TSV" \
    --output_tsv "$OUTPUT_GLYCOSYLATOR" \
    --output_dir "$OUTPUT_RESULTS"

# ---------------------------------------
# Step 2: Glycosylation script
# ---------------------------------------

PROTEIN_RESIDUE_START=10
PDB_INPUT_PROTEIN_PREPARED="$INPUT_PATH/spike_Delta.pdb"
PDB_OUTPUT_PROTEIN_PREPARED_2="$OUTPUT_PATH/PDB_PROTEIN_GLYCOSYLATED/spike_Delta_asn_orientation.pdb"
GLYCOSYLATED_OUTPUT_PDB="$OUTPUT_PATH/PDB_PROTEIN_GLYCOSYLATED/spike_glycosylated.pdb"

# ---------------------------------------
# Preparation of the asparagine orientations
# ---------------------------------------

echo "Step 2a: Optimizing asparagine orientations..."
python3 "$SCRIPTS_PYTHON"/asn_orientation.py "$PDB_INPUT_PROTEIN_PREPARED" \
    --rotate-atoms "OD1,CG,ND2,HD22,HD21,HB2,HB3" \
    --fixed-atom CB \
    --center-atom CA \
    --radius 30.0 \
    --rotation-step 1 \
    -o "$PDB_OUTPUT_PROTEIN_PREPARED_2"

echo "Step 2b: Running glycosylation script..."
python3 "$SCRIPTS_PYTHON"/2-glycosylation_script.py \
    --input_tsv_glycosylator "$OUTPUT_GLYCOSYLATOR" \
    --input_pdb_protein "$PDB_OUTPUT_PROTEIN_PREPARED_2" \
    --protein_residue_start "$PROTEIN_RESIDUE_START" \
    --output "$GLYCOSYLATED_OUTPUT_PDB"

# ---------------------------------------------------
# Step 3: Correct chain labels and residue numbers
# ---------------------------------------------------
RENAMED_PDB="$OUTPUT_PATH/PDB_PROTEIN_GLYCOSYLATED/spike_glycosylated_renumbered.pdb"

echo "Step 3: Correcting chain labels and residue numbers..."
python3 "$SCRIPTS_PYTHON"/3-correction_chain_labels_residue_numbers.py \
    "$GLYCOSYLATED_OUTPUT_PDB" \
    "$RENAMED_PDB"

# ---------------------------------------------------
# Step 4: Extract coordinates of glycans from PDB
# ---------------------------------------------------
INPUT_PDB_FINAL="$RENAMED_PDB"
OUTPUT_NOH_PDB="$OUTPUT_PATH/PDB_PROTEIN_GLYCOSYLATED/spike_glycosylated_renumbered_without_H.pdb"
OUTPUT_GLYCAN_DIR="$OUTPUT_PATH/TO_TOP"

echo "Step 4: Extracting glycans coordinates..."
python3 "$SCRIPTS_PYTHON"/4-extract_coordinates_of_glycans_from_structure.py \
    --input_pdb "$INPUT_PDB_FINAL" \
    --output_noH "$OUTPUT_NOH_PDB" \
    --output_dir "$OUTPUT_GLYCAN_DIR"

echo "Workflow completed successfully!"
