#!/bin/bash
set -euo pipefail

############################################
# USAGE FUNCTION
############################################
usage() {
    echo "Usage: $0 --input_dir <step2_directory> --output_dir <step3_results_directory>"
    echo "  --input_dir  : Directory containing STEP2 output (topology files)"
    echo "  --output_dir : Directory for STEP3 results (minimization results)"
    exit 1
}

############################################
# PARSE COMMAND LINE ARGUMENTS
############################################
STEP2_DIR=""
STEP3_RESULTS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            STEP2_DIR="$2"
            shift 2
            ;;
        --output_dir)
            STEP3_RESULTS="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$STEP2_DIR" ]] || [[ -z "$STEP3_RESULTS" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Convert to absolute paths
STEP2_DIR=$(realpath "$STEP2_DIR")
STEP3_RESULTS=$(realpath "$STEP3_RESULTS")

# Check if input directory exists
if [[ ! -d "$STEP2_DIR" ]]; then
    echo "ERROR: Input directory not found: $STEP2_DIR"
    exit 1
fi

############################################
# DEFINE PATHS (Generic)
############################################

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$SCRIPT_DIR"

# CHARMM36 directory (inside STEP2)
CHARMM36="$STEP2_DIR/charmm36.ff"
if [[ ! -d "$CHARMM36" ]]; then
    echo "ERROR: CHARMM36 directory not found: $CHARMM36"
    exit 1
fi

# JSON directory (inside STEP2)
JSON="$STEP2_DIR/JSON"
if [[ ! -d "$JSON" ]]; then
    echo "ERROR: JSON directory not found: $JSON"
    exit 1
fi

# Find the glycoprotein variants directory (generic)
GLYCOPROTEIN_STEP2_DIR="$STEP2_DIR/VALENCE_GLYCAN_VARIANTS"
if [[ ! -d "$GLYCOPROTEIN_STEP2_DIR" ]]; then
    echo "WARNING: VALENCE_GLYCAN_VARIANTS directory not found, looking for variants..."
    # Try to find any directory with variants
    GLYCOPROTEIN_STEP2_DIR=$(find "$STEP2_DIR" -type d -name "*VARIANTS*" | head -n 1)
    if [[ -z "$GLYCOPROTEIN_STEP2_DIR" ]]; then
        # Use JSON directory as fallback
        GLYCOPROTEIN_STEP2_DIR="$JSON"
    fi
fi

# Find the protein PDB file (generic - look for *_variants.pdb or *_final_valence_corrected.pdb)
PROTEIN_STEP2=$(find "$GLYCOPROTEIN_STEP2_DIR" -maxdepth 1 -name "*_variants.pdb" | head -n 1)
if [[ -z "$PROTEIN_STEP2" ]]; then
    PROTEIN_STEP2=$(find "$STEP2_DIR/PDB_GLYCOPROTEIN" -maxdepth 1 -name "*_final_valence_corrected.pdb" | head -n 1)
fi
if [[ -z "$PROTEIN_STEP2" ]]; then
    PROTEIN_STEP2=$(find "$STEP2_DIR/PDB_GLYCOPROTEIN" -maxdepth 1 -name "*.pdb" | grep -v "without_H" | head -n 1)
fi

if [[ -z "$PROTEIN_STEP2" ]] || [[ ! -f "$PROTEIN_STEP2" ]]; then
    echo "ERROR: Could not find protein PDB file in $STEP2_DIR"
    echo "Looking for: *_variants.pdb or *_final_valence_corrected.pdb"
    exit 1
fi

echo "Found input protein: $PROTEIN_STEP2"

# Extract basename for output files
BASENAME=$(basename "$PROTEIN_STEP2" .pdb)
# Remove suffixes if present
BASENAME=${BASENAME%_variants}
BASENAME=${BASENAME%_final_valence_corrected}
BASENAME=${BASENAME%_glycosylated_final_valence_corrected}

echo "Using basename: $BASENAME"

# Create output directories
mkdir -p "$STEP3_RESULTS/JSON_FILES"
mkdir -p "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED"
mkdir -p "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PDB_CARB_ONLY"

# Define output files with generic names
PDB_TO_JSON_OUTPUT="$STEP3_RESULTS/JSON_FILES/pdb_to_json.json"
CHARMM36_JSON_OUTPUT="$STEP3_RESULTS/JSON_FILES/glycan_data_charmm36.json"
OPTIMIZED_JSON_OUTPUT="$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/glycan_optimized.json"
OPTIMIZED_PDB_OUTPUT="$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/${BASENAME}_optimized.pdb"
REPORT_FILE="$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/report.txt"
GLYCANS_OUTPUT_DIR="$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PDB_CARB_ONLY"

############################################
# CHECK PYTHON SCRIPTS
############################################

# List of required Python scripts
REQUIRED_SCRIPTS=(
    "1-pdb_to_json.py"
    "3-adding_chamm36_parameters.py"
    "4-optimize_glycans_mcmc.py"
)

for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [[ ! -f "$SCRIPTS/$script" ]]; then
        echo "ERROR: Python script not found: $SCRIPTS/$script"
        exit 1
    fi
done

############################################
# STEP 1: Convert PDB to JSON
############################################
echo "========================================"
echo "Step 1: Converting PDB to JSON..."
echo "========================================"
echo "Input PDB: $PROTEIN_STEP2"
echo "Output JSON: $PDB_TO_JSON_OUTPUT"

python3 "$SCRIPTS/1-pdb_to_json.py" \
    --input_pdb "$PROTEIN_STEP2" \
    --output_json "$PDB_TO_JSON_OUTPUT"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 1 failed"
    exit 1
fi

############################################
# STEP 2: Add CHARMM36 parameters
############################################
echo "========================================"
echo "Step 2: Adding CHARMM36 parameters..."
echo "========================================"
echo "Input JSON: $PDB_TO_JSON_OUTPUT"
echo "CHARMM directory: $CHARMM36"
echo "Output JSON: $CHARMM36_JSON_OUTPUT"

python3 "$SCRIPTS/3-adding_chamm36_parameters.py" \
    --input_json "$PDB_TO_JSON_OUTPUT" \
    --charmm_dir "$CHARMM36" \
    --output_json "$CHARMM36_JSON_OUTPUT"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 2 failed"
    exit 1
fi

############################################
# STEP 3: Optimize glycans using MCMC
############################################
echo "========================================"
echo "Step 3: Optimizing glycans using MCMC..."
echo "========================================"
echo "Input JSON: $CHARMM36_JSON_OUTPUT"
echo "Output JSON: $OPTIMIZED_JSON_OUTPUT"
echo "Output PDB: $OPTIMIZED_PDB_OUTPUT"
echo "Glycans output dir: $GLYCANS_OUTPUT_DIR"
echo "Report file: $REPORT_FILE"

python3 "$SCRIPTS/4-optimize_glycans_mcmc.py" \
    --input_json "$CHARMM36_JSON_OUTPUT" \
    --output_json "$OPTIMIZED_JSON_OUTPUT" \
    --output_pdb "$OPTIMIZED_PDB_OUTPUT" \
    --glycans_output_dir "$GLYCANS_OUTPUT_DIR" \
    --theta_step 10 \
    --n_steps 10 \
    --max_cycles 5 \
    --radius 300 \
    --use_coulomb no \
    --n_workers 12 \
    --report_file "$REPORT_FILE" \
    --save_individual_glycans \
    --save_before_after \
    --verbose

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 3 failed"
    exit 1
fi

############################################
# COMPLETION
############################################
echo "========================================"
echo "All steps completed successfully!"
echo "========================================"
echo "Results saved in: $STEP3_RESULTS"
echo "  - JSON_FILES/: JSON configuration files"
echo "  - PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/: Optimized structures"
echo "  - PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PDB_CARB_ONLY/: Individual glycans"
echo ""
echo "Main output files:"
echo "  - Optimized PDB: $OPTIMIZED_PDB_OUTPUT"
echo "  - Report: $REPORT_FILE"
echo "========================================"
