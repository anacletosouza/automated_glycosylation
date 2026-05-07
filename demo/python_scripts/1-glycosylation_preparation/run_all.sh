#!/bin/bash
# -----------------------------------------------------------------------------
# run_all.sh
# Shell script to run the full glycosylation workflow
# Each step corresponds to a separate Python script execution
# -----------------------------------------------------------------------------

set -e  # Stop script if any command fails

# -------------------------------
# Parse command line arguments
# -------------------------------
usage() {
    echo "Usage: $0 --pdb <input.pdb> --tsv <input.tsv> --output <output_dir>"
    echo "  --pdb     : Input PDB file path"
    echo "  --tsv     : Input TSV file path"
    echo "  --output  : Output directory path"
    exit 1
}

# Initialize variables
INPUT_PDB=""
INPUT_TSV=""
OUTPUT_BASE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --pdb)
            INPUT_PDB="$2"
            shift 2
            ;;
        --tsv)
            INPUT_TSV="$2"
            shift 2
            ;;
        --output)
            OUTPUT_BASE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_PDB" ]] || [[ -z "$INPUT_TSV" ]] || [[ -z "$OUTPUT_BASE" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Convert to absolute paths
INPUT_PDB=$(realpath "$INPUT_PDB")
INPUT_TSV=$(realpath "$INPUT_TSV")
OUTPUT_BASE=$(realpath "$OUTPUT_BASE")

# Check if input files exist
if [[ ! -f "$INPUT_PDB" ]]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

if [[ ! -f "$INPUT_TSV" ]]; then
    echo "ERROR: Input TSV not found: $INPUT_TSV"
    exit 1
fi

# -------------------------------
# Define paths
# -------------------------------
# Get the absolute path where the Python scripts are located
# This assumes the run_all.sh script is in the same directory as the python_scripts folder
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_PYTHON="${SCRIPT_DIR}"

# If the python scripts are in a subdirectory, uncomment and adjust the line below
# SCRIPTS_PYTHON="${SCRIPT_DIR}/python_scripts"

# Output directories
OUTPUT_PATH="$OUTPUT_BASE"
TSV_DIR="$OUTPUT_PATH/TSV"
PDB_DIR="$OUTPUT_PATH/PDB_PROTEIN_GLYCOSYLATED"
GLYCAN_DIR="$OUTPUT_PATH/EXTRACTED_CARBOHYDRATES"
TO_TOP_DIR="$OUTPUT_PATH/TO_TOP"

# -------------------------------
# Create output directories
# -------------------------------
echo "Creating output directories..."
mkdir -p "$TSV_DIR"
mkdir -p "$PDB_DIR"
mkdir -p "$GLYCAN_DIR"
mkdir -p "$TO_TOP_DIR"

# Extract base filename without extension
BASENAME=$(basename "$INPUT_TSV" .tsv)
BASENAME_PDB=$(basename "$INPUT_PDB" .pdb)

# -------------------------------
# Step 0: Correct input table
# -------------------------------
INPUT_TSV_FILE="$INPUT_TSV"
OUTPUT_TSV="$TSV_DIR/${BASENAME}_corrected.tsv"
OUTPUT_RESULTS="$GLYCAN_DIR"

echo "Step 0: Correcting input table..."
echo "Python script path: ${SCRIPTS_PYTHON}/0-correcting_caselino_table_for_variants.py"

if [[ ! -f "${SCRIPTS_PYTHON}/0-correcting_caselino_table_for_variants.py" ]]; then
    echo "ERROR: Python script not found: ${SCRIPTS_PYTHON}/0-correcting_caselino_table_for_variants.py"
    exit 1
fi

python3 "${SCRIPTS_PYTHON}/0-correcting_caselino_table_for_variants.py" \
    --input "$INPUT_TSV_FILE" \
    --output "$OUTPUT_TSV"

# ---------------------------------------
# Step 1: Convert to IUPAC notation
# ---------------------------------------
OUTPUT_GLYCOSYLATOR="$TSV_DIR/${BASENAME}_glycosylator.tsv"

echo "Step 1: Converting to IUPAC notation..."
if [[ ! -f "${SCRIPTS_PYTHON}/1-iupac_converted.py" ]]; then
    echo "ERROR: Python script not found: ${SCRIPTS_PYTHON}/1-iupac_converted.py"
    exit 1
fi

python3 "${SCRIPTS_PYTHON}/1-iupac_converted.py" \
    --input_tsv "$OUTPUT_TSV" \
    --output_tsv "$OUTPUT_GLYCOSYLATOR" \
    --output_dir "$OUTPUT_RESULTS"

# ---------------------------------------
# Step 2: Glycosylation script
# ---------------------------------------
PROTEIN_RESIDUE_START=10
PDB_INPUT_PROTEIN_PREPARED="$INPUT_PDB"
PDB_OUTPUT_PROTEIN_PREPARED_2="$PDB_DIR/${BASENAME_PDB}_asn_orientation.pdb"
GLYCOSYLATED_OUTPUT_PDB="$PDB_DIR/${BASENAME_PDB}_glycosylated.pdb"

# ---------------------------------------
# Preparation of the asparagine orientations
# ---------------------------------------
echo "Step 2a: Optimizing asparagine orientations..."
if [[ ! -f "${SCRIPTS_PYTHON}/asn_orientation.py" ]]; then
    echo "ERROR: Python script not found: ${SCRIPTS_PYTHON}/asn_orientation.py"
    exit 1
fi

python3 "${SCRIPTS_PYTHON}/asn_orientation.py" "$PDB_INPUT_PROTEIN_PREPARED" \
    --rotate-atoms "OD1,CG,ND2,HD22,HD21,HB2,HB3" \
    --fixed-atom CB \
    --center-atom CA \
    --radius 30.0 \
    --rotation-step 1 \
    -o "$PDB_OUTPUT_PROTEIN_PREPARED_2"

echo "Step 2b: Running glycosylation script..."
if [[ ! -f "${SCRIPTS_PYTHON}/2-glycosylation_script.py" ]]; then
    echo "ERROR: Python script not found: ${SCRIPTS_PYTHON}/2-glycosylation_script.py"
    exit 1
fi

python3 "${SCRIPTS_PYTHON}/2-glycosylation_script.py" \
    --input_tsv_glycosylator "$OUTPUT_GLYCOSYLATOR" \
    --input_pdb_protein "$PDB_OUTPUT_PROTEIN_PREPARED_2" \
    --protein_residue_start "$PROTEIN_RESIDUE_START" \
    --output "$GLYCOSYLATED_OUTPUT_PDB"

# ---------------------------------------------------
# Step 3: Correct chain labels and residue numbers
# ---------------------------------------------------
RENAMED_PDB="$PDB_DIR/${BASENAME_PDB}_glycosylated_renumbered.pdb"

echo "Step 3: Correcting chain labels and residue numbers..."
if [[ ! -f "${SCRIPTS_PYTHON}/3-correction_chain_labels_residue_numbers.py" ]]; then
    echo "ERROR: Python script not found: ${SCRIPTS_PYTHON}/3-correction_chain_labels_residue_numbers.py"
    exit 1
fi

python3 "${SCRIPTS_PYTHON}/3-correction_chain_labels_residue_numbers.py" \
    "$GLYCOSYLATED_OUTPUT_PDB" \
    "$RENAMED_PDB"

# ---------------------------------------------------
# Step 4: Extract coordinates of glycans from PDB
# ---------------------------------------------------
INPUT_PDB_FINAL="$RENAMED_PDB"
OUTPUT_NOH_PDB="$PDB_DIR/${BASENAME_PDB}_glycosylated_renumbered_without_H.pdb"
OUTPUT_GLYCAN_DIR="$TO_TOP_DIR"

echo "Step 4: Extracting glycans coordinates..."
if [[ ! -f "${SCRIPTS_PYTHON}/4-extract_coordinates_of_glycans_from_structure.py" ]]; then
    echo "ERROR: Python script not found: ${SCRIPTS_PYTHON}/4-extract_coordinates_of_glycans_from_structure.py"
    exit 1
fi

python3 "${SCRIPTS_PYTHON}/4-extract_coordinates_of_glycans_from_structure.py" \
    --input_pdb "$INPUT_PDB_FINAL" \
    --output_noH "$OUTPUT_NOH_PDB" \
    --output_dir "$OUTPUT_GLYCAN_DIR"

echo "Workflow completed successfully!"
echo "Output files are in: $OUTPUT_PATH"
