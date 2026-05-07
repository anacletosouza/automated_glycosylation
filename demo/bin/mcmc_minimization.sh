#!/bin/bash
set -euo pipefail

# -----------------------------------------------------------------------------
# mcmc_minimization.sh
# Unified script to run the complete glycosylation workflow:
#   Step 1: Glycosylation preparation
#   Step 2: Parametrization scripts
#   Step 3: Carbohydrate orientation optimization
# -----------------------------------------------------------------------------

############################################
# USAGE FUNCTION
############################################
usage() {
    echo "Usage: $0 --pdb <input.pdb> --tsv <input.tsv> [OPTIONS]"
    echo ""
    echo "REQUIRED:"
    echo "  --pdb     : Input PDB file path"
    echo "  --tsv     : Input TSV file path"
    echo ""
    echo "OPTIONAL:"
    echo "  --output_dir <dir>               : Output directory (default: ./output)"
    echo "  --input_python_scripts_dir <dir> : Directory containing Python scripts"
    echo "                                     (default: directory of this script)"
    echo "  --url_charmm36 <url>             : URL for CHARMM36 force field download"
    echo "                                     (default: https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz)"
    echo ""
    echo "  --theta_step <deg>                : Rotation step in degrees (default: 10)"
    echo "  --n_steps <int>                   : Number of rotation steps (default: 10)"
    echo "  --max_cycles <int>                : Maximum MCMC cycles (default: 5)"
    echo "  --radius <float>                  : Radius for clash detection (default: 300)"
    echo "  --use_coulomb <yes/no>            : Use Coulomb interactions (default: no)"
    echo "  --n_workers <int>                 : Number of parallel workers (default: 12)"
    echo "  --report_file <path>              : Report file path (default: <output_dir>/report.txt)"
    echo "  --save_individual_glycans         : Save individual glycan PDB files"
    echo "  --save_before_after               : Save before/after comparison files"
    echo "  --verbose                         : Enable verbose output"
    echo ""
    echo "  --help                            : Show this help message"
    exit 1
}

############################################
# HELPER FUNCTION TO FIND PYTHON SCRIPT
############################################
find_python_script() {
    local script_name="$1"
    local search_dir="$2"
    
    # First check in the main directory
    if [[ -f "${search_dir}/${script_name}" ]]; then
        echo "${search_dir}/${script_name}"
        return 0
    fi
    
    # Search in subdirectories
    local found=$(find "$search_dir" -name "$script_name" -type f | head -n 1)
    if [[ -n "$found" ]]; then
        echo "$found"
        return 0
    fi
    
    # Search in conventional subdirectory names
    for subdir in "1-glycosylation_preparation" "2-parametrization_scripts" "3-carbohydrate_orientation"; do
        if [[ -f "${search_dir}/${subdir}/${script_name}" ]]; then
            echo "${search_dir}/${subdir}/${script_name}"
            return 0
        fi
    done
    
    echo ""
    return 1
}

############################################
# PARSE COMMAND LINE ARGUMENTS
############################################

# Required arguments
INPUT_PDB=""
INPUT_TSV=""

# Optional arguments with defaults
OUTPUT_BASE=""
PYTHON_SCRIPTS_DIR=""
URL_CHARMM36="https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz"

# Step 3 parameters
THETA_STEP=10
N_STEPS=10
MAX_CYCLES=5
RADIUS=300
USE_COULOMB="no"
N_WORKERS=12
REPORT_FILE=""
SAVE_INDIVIDUAL_GLYCANS=""
SAVE_BEFORE_AFTER=""
VERBOSE=""

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
        --output_dir)
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --input_python_scripts_dir)
            PYTHON_SCRIPTS_DIR="$2"
            shift 2
            ;;
        --url_charmm36)
            URL_CHARMM36="$2"
            shift 2
            ;;
        --theta_step)
            THETA_STEP="$2"
            shift 2
            ;;
        --n_steps)
            N_STEPS="$2"
            shift 2
            ;;
        --max_cycles)
            MAX_CYCLES="$2"
            shift 2
            ;;
        --radius)
            RADIUS="$2"
            shift 2
            ;;
        --use_coulomb)
            USE_COULOMB="$2"
            shift 2
            ;;
        --n_workers)
            N_WORKERS="$2"
            shift 2
            ;;
        --report_file)
            REPORT_FILE="$2"
            shift 2
            ;;
        --save_individual_glycans)
            SAVE_INDIVIDUAL_GLYCANS="--save_individual_glycans"
            shift
            ;;
        --save_before_after)
            SAVE_BEFORE_AFTER="--save_before_after"
            shift
            ;;
        --verbose)
            VERBOSE="--verbose"
            shift
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_PDB" ]] || [[ -z "$INPUT_TSV" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Set default output directory if not provided
if [[ -z "$OUTPUT_BASE" ]]; then
    OUTPUT_BASE="./output"
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

# Set Python scripts directory
if [[ -z "$PYTHON_SCRIPTS_DIR" ]]; then
    # Get the absolute path where this script is located, then go up to find python_scripts
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    # Assuming the script is in ../bin/ and python_scripts is in ../
    PYTHON_SCRIPTS_DIR="$(dirname "$SCRIPT_DIR")/python_scripts"
fi

if [[ ! -d "$PYTHON_SCRIPTS_DIR" ]]; then
    echo "ERROR: Python scripts directory not found: $PYTHON_SCRIPTS_DIR"
    exit 1
fi

PYTHON_SCRIPTS_DIR=$(realpath "$PYTHON_SCRIPTS_DIR")

# Set default report file if not provided
if [[ -z "$REPORT_FILE" ]]; then
    REPORT_FILE="$OUTPUT_BASE/STEP3/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/report.txt"
else
    REPORT_FILE=$(realpath "$REPORT_FILE")
fi

echo "========================================"
echo "Configuration:"
echo "  Input PDB: $INPUT_PDB"
echo "  Input TSV: $INPUT_TSV"
echo "  Output directory: $OUTPUT_BASE"
echo "  Python scripts dir: $PYTHON_SCRIPTS_DIR"
echo "  CHARMM36 URL: $URL_CHARMM36"
echo "========================================"

############################################
# CREATE OUTPUT DIRECTORY STRUCTURE
############################################

# Step 1 directories
STEP1_OUTPUT="$OUTPUT_BASE/STEP1"
TSV_DIR="$STEP1_OUTPUT/TSV"
PDB_DIR="$STEP1_OUTPUT/PDB_PROTEIN_GLYCOSYLATED"
GLYCAN_DIR="$STEP1_OUTPUT/EXTRACTED_CARBOHYDRATES"
TO_TOP_DIR="$STEP1_OUTPUT/TO_TOP"

# Step 2 directories
STEP2_OUTPUT="$OUTPUT_BASE/STEP2"
TOPO_DIR="$STEP2_OUTPUT"
JSON_DIR="$TOPO_DIR/JSON"
PDB_GLYCO_DIR="$TOPO_DIR/PDB_GLYCOPROTEIN"
VARIANTS_DIR="$TOPO_DIR/VALENCE_GLYCAN_VARIANTS"

# Step 3 directories
STEP3_OUTPUT="$OUTPUT_BASE/STEP3"
JSON_FILES_DIR="$STEP3_OUTPUT/JSON_FILES"
OPTIMIZED_DIR="$STEP3_OUTPUT/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED"
GLYCANS_ONLY_DIR="$OPTIMIZED_DIR/PDB_CARB_ONLY"

# Create all directories
mkdir -p "$TSV_DIR"
mkdir -p "$PDB_DIR"
mkdir -p "$GLYCAN_DIR"
mkdir -p "$TO_TOP_DIR"
mkdir -p "$JSON_DIR"
mkdir -p "$PDB_GLYCO_DIR"
mkdir -p "$VARIANTS_DIR"
mkdir -p "$JSON_FILES_DIR"
mkdir -p "$OPTIMIZED_DIR"
mkdir -p "$GLYCANS_ONLY_DIR"

# Extract base filename without extension
BASENAME=$(basename "$INPUT_TSV" .tsv)
BASENAME_PDB=$(basename "$INPUT_PDB" .pdb)

echo "========================================"
echo "STEP 1: GLYCOSYLATION PREPARATION"
echo "========================================"

# ---------------------------------------
# Find required scripts for Step 1
# ---------------------------------------
SCRIPT_CORRECT=$(find_python_script "0-correcting_caselino_table_for_variants.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_CORRECT" ]]; then
    echo "ERROR: Could not find 0-correcting_caselino_table_for_variants.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_CORRECT"

SCRIPT_IUPAC=$(find_python_script "1-iupac_converted.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_IUPAC" ]]; then
    echo "ERROR: Could not find 1-iupac_converted.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_IUPAC"

SCRIPT_ASN=$(find_python_script "asn_orientation.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_ASN" ]]; then
    echo "ERROR: Could not find asn_orientation.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_ASN"

SCRIPT_GLYCOSYLATION=$(find_python_script "2-glycosylation_script.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_GLYCOSYLATION" ]]; then
    echo "ERROR: Could not find 2-glycosylation_script.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_GLYCOSYLATION"

SCRIPT_CORRECTION=$(find_python_script "3-correction_chain_labels_residue_numbers.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_CORRECTION" ]]; then
    echo "ERROR: Could not find 3-correction_chain_labels_residue_numbers.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_CORRECTION"

SCRIPT_EXTRACT=$(find_python_script "4-extract_coordinates_of_glycans_from_structure.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_EXTRACT" ]]; then
    echo "ERROR: Could not find 4-extract_coordinates_of_glycans_from_structure.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_EXTRACT"

# ---------------------------------------
# Step 1.0: Correct input table
# ---------------------------------------
INPUT_TSV_FILE="$INPUT_TSV"
OUTPUT_TSV="$TSV_DIR/${BASENAME}_corrected.tsv"
OUTPUT_RESULTS="$GLYCAN_DIR"

echo "Step 1.0: Correcting input table..."

python3 "$SCRIPT_CORRECT" \
    --input "$INPUT_TSV_FILE" \
    --output "$OUTPUT_TSV"

# ---------------------------------------
# Step 1.1: Convert to IUPAC notation
# ---------------------------------------
OUTPUT_GLYCOSYLATOR="$TSV_DIR/${BASENAME}_glycosylator.tsv"

echo "Step 1.1: Converting to IUPAC notation..."

python3 "$SCRIPT_IUPAC" \
    --input_tsv "$OUTPUT_TSV" \
    --output_tsv "$OUTPUT_GLYCOSYLATOR" \
    --output_dir "$OUTPUT_RESULTS"

# ---------------------------------------
# Step 1.2a: Preparation of asparagine orientations
# ---------------------------------------
PROTEIN_RESIDUE_START=10
PDB_INPUT_PROTEIN_PREPARED="$INPUT_PDB"
PDB_OUTPUT_PROTEIN_PREPARED_2="$PDB_DIR/${BASENAME_PDB}_asn_orientation.pdb"
GLYCOSYLATED_OUTPUT_PDB="$PDB_DIR/${BASENAME_PDB}_glycosylated.pdb"

echo "Step 1.2a: Optimizing asparagine orientations..."

python3 "$SCRIPT_ASN" "$PDB_INPUT_PROTEIN_PREPARED" \
    --rotate-atoms "OD1,CG,ND2,HD22,HD21,HB2,HB3" \
    --fixed-atom CB \
    --center-atom CA \
    --radius 30.0 \
    --rotation-step 1 \
    -o "$PDB_OUTPUT_PROTEIN_PREPARED_2"

# ---------------------------------------
# Step 1.2b: Run glycosylation script
# ---------------------------------------
echo "Step 1.2b: Running glycosylation script..."

python3 "$SCRIPT_GLYCOSYLATION" \
    --input_tsv_glycosylator "$OUTPUT_GLYCOSYLATOR" \
    --input_pdb_protein "$PDB_OUTPUT_PROTEIN_PREPARED_2" \
    --protein_residue_start "$PROTEIN_RESIDUE_START" \
    --output "$GLYCOSYLATED_OUTPUT_PDB"

# ---------------------------------------
# Step 1.3: Correct chain labels and residue numbers
# ---------------------------------------
RENAMED_PDB="$PDB_DIR/${BASENAME_PDB}_glycosylated_renumbered.pdb"

echo "Step 1.3: Correcting chain labels and residue numbers..."

python3 "$SCRIPT_CORRECTION" \
    "$GLYCOSYLATED_OUTPUT_PDB" \
    "$RENAMED_PDB"

# ---------------------------------------
# Step 1.4: Extract coordinates of glycans from PDB
# ---------------------------------------
INPUT_PDB_FINAL="$RENAMED_PDB"
OUTPUT_NOH_PDB="$PDB_DIR/${BASENAME_PDB}_glycosylated_renumbered_without_H.pdb"
OUTPUT_GLYCAN_DIR="$TO_TOP_DIR"

echo "Step 1.4: Extracting glycans coordinates..."

python3 "$SCRIPT_EXTRACT" \
    --input_pdb "$INPUT_PDB_FINAL" \
    --output_noH "$OUTPUT_NOH_PDB" \
    --output_dir "$OUTPUT_GLYCAN_DIR"

echo "STEP 1 completed successfully!"

echo "========================================"
echo "STEP 2: PARAMETRIZATION SCRIPTS"
echo "========================================"

# ---------------------------------------
# Find required scripts for Step 2
# ---------------------------------------
SCRIPT_JSON_GEN=$(find_python_script "0-JSON_generator.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_JSON_GEN" ]]; then
    echo "ERROR: Could not find 0-JSON_generator.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_JSON_GEN"

SCRIPT_PARSER=$(find_python_script "1-parser_pdb.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_CARB_RTP=$(find_python_script "2-parser_carb_rtp.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_COMPARE=$(find_python_script "3-comparison_pdb_rtp.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_RTP1=$(find_python_script "4-rtp_generator_part1.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_RTP2=$(find_python_script "4-rtp_generator_part2.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_RTP3=$(find_python_script "4-rtp_generator_part3.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_ACETYL=$(find_python_script "5-acetylation_replacement.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_CLEAN=$(find_python_script "6-clean_rtp.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_TOGETHER1=$(find_python_script "7-together_part_1.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_TOGETHER2=$(find_python_script "7-together_part_2.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_TOGETHER3=$(find_python_script "7-together_part_3.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_TOGETHER4=$(find_python_script "7-together_part_4.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_TOGETHER5=$(find_python_script "7-together_part_5.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_GLYCOPROT=$(find_python_script "8-glycoprotein.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_CONNECT=$(find_python_script "9-conection_glycosilation_without_TER.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_IDENTIFY=$(find_python_script "glycosylation_identifying.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_HDB=$(find_python_script "10-generation_hdb.py" "$PYTHON_SCRIPTS_DIR")
SCRIPT_VARIANTS=$(find_python_script "glycosylation_variants.py" "$PYTHON_SCRIPTS_DIR")

# CHARMM Force Field Setup
CHARMM_DIR="$TOPO_DIR/charmm36.ff"
CHARMM_RTP="$CHARMM_DIR/carb.rtp"
CHARMM_HDB="$CHARMM_DIR/carb.hdb"
CHARMM_HDB_BACKUP="$CHARMM_DIR/carb.hdb.backup"
CHARMM_RTP_BACKUP="$CHARMM_DIR/carb.rtp.backup"

cd "$TOPO_DIR"

# (1) Download or restore CHARMM force field
if [ ! -f "$CHARMM_HDB_BACKUP" ] || [ ! -f "$CHARMM_RTP_BACKUP" ]; then
    echo "Backup files not found. Downloading CHARMM force field..."
    
    rm -rf "$CHARMM_DIR"
    
    wget "$URL_CHARMM36" -O charmm36.ff.tgz
    
    tar -xzf charmm36.ff.tgz
    
    # Find the extracted directory name (may vary)
    EXTRACTED_DIR=$(find . -maxdepth 1 -type d -name "charmm36-*" | head -n 1)
    if [[ -z "$EXTRACTED_DIR" ]]; then
        EXTRACTED_DIR=$(find . -maxdepth 1 -type d -name "charmm36*" | head -n 1)
    fi
    
    if [[ -n "$EXTRACTED_DIR" ]]; then
        mv "$EXTRACTED_DIR" "$CHARMM_DIR"
    else
        echo "ERROR: Could not find extracted CHARMM directory"
        exit 1
    fi
    
    rm -f charmm36.ff.tgz
    
    # Create initial backups
    cp "$CHARMM_HDB" "$CHARMM_HDB_BACKUP"
    cp "$CHARMM_RTP" "$CHARMM_RTP_BACKUP"
    
    echo "CHARMM force field downloaded and backups created at: $CHARMM_DIR"
else
    echo "Backup files found. Restoring from backups..."
    
    mkdir -p "$CHARMM_DIR"
    
    cp "$CHARMM_HDB_BACKUP" "$CHARMM_HDB"
    cp "$CHARMM_RTP_BACKUP" "$CHARMM_RTP"
    
    echo "Force field restored from backups."
fi

# Restore residuetypes.dat if backup exists
if [ -f "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat.backup" ]; then
    cp "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat.backup" \
       "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat" 2>/dev/null || true
fi

# Check required files exist
for FILE in "$CHARMM_RTP" "$CHARMM_HDB"; do
    if [ ! -f "$FILE" ]; then
        echo "ERROR: File not found: $FILE"
        exit 1
    fi
done

# ---------------------------------------
# Generation of JSONs for each glycan
# ---------------------------------------
echo "Generating JSON for each glycan..."

python3 "$SCRIPT_JSON_GEN" \
    --base_dir "$STEP1_OUTPUT" \
    --output_dir "$JSON_DIR"

# ---------------------------------------
# Process each glycan directory
# ---------------------------------------
echo "Starting processing for each glycan directory..."

for DIR in "$JSON_DIR"/*; do
    [ -d "$DIR" ] || continue

    DIR_BASENAME=$(basename "$DIR")

    PDB_FILE="$DIR/$DIR_BASENAME.pdb"
    JSON_FILE="$DIR/$DIR_BASENAME.json"
    PARSER_FILE="$DIR/${DIR_BASENAME}_parser.pkl"
    RTP_PICKLE="$DIR/carb_residues.pkl"
    RTP_MODIFIED="$DIR/carb_modified.rtp"
    RTP_UNIQUE="$DIR/carb_unique.rtp"

    if [[ ! -f "$PDB_FILE" ]]; then
        echo "WARNING: PDB file not found: $PDB_FILE, skipping..."
        continue
    fi

    echo "----------------------------------------"
    echo "Processing: $DIR_BASENAME"

    cd "$DIR"

    python3 "$SCRIPT_PARSER" "$PDB_FILE" -o "$PARSER_FILE" || {
        echo "ERROR in 1-parser_pdb.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPT_CARB_RTP" "$CHARMM_RTP" -o "$RTP_PICKLE" || {
        echo "ERROR in 2-parser_carb_rtp.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPT_COMPARE" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" || {
        echo "ERROR in 3-comparison_pdb_rtp.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPT_RTP1" \
        -p "$PARSER_FILE" \
        -r "$RTP_PICKLE" || {
        echo "ERROR in 4-rtp_generator_part1.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPT_RTP2" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" \
        --json "$JSON_FILE" || {
        echo "ERROR in 4-rtp_generator_part2.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPT_RTP3" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" \
        --json "$JSON_FILE" \
        --output "$RTP_MODIFIED" || {
        echo "ERROR in 4-rtp_generator_part3.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPT_ACETYL" \
        "$PDB_FILE" \
        "${DIR}/${DIR_BASENAME}_modified.pdb" || {
        echo "ERROR in 5-acetylation_replacement.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    DIR_LETTER="${DIR_BASENAME:0:1}"

    python3 "$SCRIPT_CLEAN" \
        "$RTP_MODIFIED" \
        "$RTP_UNIQUE" \
        "$DIR_LETTER" || {
        echo "ERROR in 6-clean_rtp.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    cd - > /dev/null

    echo "Finished $DIR_BASENAME"
done

# ---------------------------------------
# Unification of RTP / HDB files
# ---------------------------------------
echo "Unifying RTP files..."

if [ -n "$(ls -A "$JSON_DIR"/*/carb_unique.rtp 2>/dev/null)" ]; then
    python3 "$SCRIPT_TOGETHER1" \
        --input "$JSON_DIR" \
        --output "$JSON_DIR/carb_unique_total.rtp" || {
        echo "ERROR in 7-together_part_1.py"
        exit 1
    }

    python3 "$SCRIPT_TOGETHER2" \
        --input "$JSON_DIR/carb_unique_total.rtp" \
        --output "$JSON_DIR/carb_redundance_removed.rtp" || {
        echo "ERROR in 7-together_part_2.py"
        exit 1
    }

    python3 "$SCRIPT_TOGETHER3" \
        --input "$JSON_DIR/carb_redundance_removed.rtp" \
        --output "$CHARMM_RTP" || {
        echo "ERROR in 7-together_part_3.py"
        exit 1
    }

    python3 "$SCRIPT_TOGETHER4" \
        --input "$JSON_DIR/carb_redundance_removed.rtp" || {
        echo "ERROR in 7-together_part_4.py"
        exit 1
    }

    python3 "$SCRIPT_TOGETHER5" \
        --input "$JSON_DIR/carb_redundance_removed.rtp" \
        --output "$JSON_DIR/carb_redundance_removed.hdb" || {
        echo "ERROR in 7-together_part_5.py"
        exit 1
    }
else
    echo "WARNING: No carb_unique.rtp files found. Skipping unification steps."
fi

# ---------------------------------------
# Final Glycoprotein Construction
# ---------------------------------------
echo "Building final glycoprotein..."

# Find the input protein from STEP1
INPUT_PROTEIN=$(find "$PDB_DIR" -name "*_glycosylated_renumbered.pdb" | head -n 1)
if [[ -z "$INPUT_PROTEIN" ]]; then
    echo "ERROR: Could not find *_glycosylated_renumbered.pdb in $PDB_DIR"
    exit 1
fi

echo "Found input protein: $INPUT_PROTEIN"

OUTPUT_PROTEIN="$PDB_GLYCO_DIR/${BASENAME_PDB}_glycosylated_corrected.pdb"
FINAL_STRUCTURE="$PDB_GLYCO_DIR/${BASENAME_PDB}_glycosylated_final_connected.pdb"
FINAL_STRUCTURE_2="$PDB_GLYCO_DIR/${BASENAME_PDB}_glycosylated_final_valence_corrected.pdb"
FINAL_STRUCTURE_2_NOH="$PDB_GLYCO_DIR/${BASENAME_PDB}_glycosylated_final_valence_corrected_noh.pdb"

if [ -n "$(ls -d "$JSON_DIR"/*/ 2>/dev/null)" ]; then
    python3 "$SCRIPT_GLYCOPROT" \
        --protein "$INPUT_PROTEIN" \
        --carbs_dir "$JSON_DIR" \
        --output "$OUTPUT_PROTEIN" \
        --keep_hydrogens_carb \
        --keep_hydrogens_prot || {
        echo "ERROR in 8-glycoprotein.py"
        exit 1
    }

    python3 "$SCRIPT_CONNECT" \
        --glycosylated "$OUTPUT_PROTEIN" \
        --conect "$INPUT_PROTEIN" \
        --output "$FINAL_STRUCTURE" || {
        echo "ERROR in 9-conection_glycosilation_without_TER.py"
        exit 1
    }
    
    python3 "$SCRIPT_IDENTIFY" "$FINAL_STRUCTURE" "$FINAL_STRUCTURE_2" || {
        echo "ERROR in glycosylation_identifying.py"
        exit 1
    }

    if [[ ! -f "$FINAL_STRUCTURE_2" ]]; then
        echo "WARNING: $FINAL_STRUCTURE_2 not found"
        if [[ -f "$FINAL_STRUCTURE_2_NOH" ]]; then
            echo "Found alternative file: $FINAL_STRUCTURE_2_NOH"
            FINAL_STRUCTURE_2="$FINAL_STRUCTURE_2_NOH"
        else
            echo "ERROR: Neither $FINAL_STRUCTURE_2 nor $FINAL_STRUCTURE_2_NOH found"
            exit 1
        fi
    fi

    echo "Final structure file: $FINAL_STRUCTURE_2"
else
    echo "WARNING: No glycan directories found in $JSON_DIR"
    echo "Skipping glycoprotein construction steps."
fi

# ---------------------------------------
# HDB Update
# ---------------------------------------
echo "Updating HDB files..."

if [[ -f "$JSON_DIR/carb_redundance_removed.rtp" ]]; then
    python3 "$SCRIPT_HDB" \
        "$JSON_DIR/carb_redundance_removed.rtp" \
        "$CHARMM_HDB" \
        -o "$JSON_DIR/carb_modified.hdb" || {
        echo "ERROR in 10-generation_hdb.py"
        exit 1
    }
else
    echo "WARNING: carb_redundance_removed.rtp not found, skipping HDB generation"
fi

# ---------------------------------------
# Variants Generation
# ---------------------------------------
echo "Generating glycan variants..."

if [[ -f "$FINAL_STRUCTURE_2" ]] && [[ -f "$JSON_DIR/carb_redundance_removed.rtp" ]] && [[ -f "$JSON_DIR/carb_redundance_removed.hdb" ]]; then
    python3 "$SCRIPT_VARIANTS" \
        -p "$FINAL_STRUCTURE_2" \
        -r "$JSON_DIR/carb_redundance_removed.rtp" \
        -d "$JSON_DIR/carb_redundance_removed.hdb" \
        -o "$VARIANTS_DIR" || {
        echo "WARNING: glycosylation_variants.py failed, but continuing..."
    }
else
    echo "WARNING: Required files for variant generation not found"
fi

# ---------------------------------------
# Include variants in CHARMM files
# ---------------------------------------
echo "Including generated variants in CHARMM force field files..."

VARIANT_RTP=$(find "$VARIANTS_DIR" -name "*_variants.rtp" 2>/dev/null | head -n 1)
if [[ -z "$VARIANT_RTP" ]]; then
    VARIANT_RTP=$(find "$VARIANTS_DIR" -name "*.rtp" 2>/dev/null | grep -v "carb_" | head -n 1)
fi

VARIANT_HDB=$(find "$VARIANTS_DIR" -name "*_variants.hdb" 2>/dev/null | head -n 1)
if [[ -z "$VARIANT_HDB" ]]; then
    VARIANT_HDB=$(find "$VARIANTS_DIR" -name "*.hdb" 2>/dev/null | grep -v "carb_" | head -n 1)
fi

if [[ -f "$VARIANT_RTP" ]] && [[ -f "$VARIANT_HDB" ]]; then
    echo "Found variant RTP: $VARIANT_RTP"
    echo "Found variant HDB: $VARIANT_HDB"
    
    LINES_BEFORE_RTP=$(wc -l < "$CHARMM_RTP")
    LINES_BEFORE_HDB=$(wc -l < "$CHARMM_HDB")
    
    echo "Appending variant RTP to $CHARMM_RTP..."
    cat "$VARIANT_RTP" >> "$CHARMM_RTP"
    
    echo "Appending variant HDB to $CHARMM_HDB..."
    cat "$VARIANT_HDB" >> "$CHARMM_HDB"
    
    LINES_AFTER_RTP=$(wc -l < "$CHARMM_RTP")
    LINES_AFTER_HDB=$(wc -l < "$CHARMM_HDB")
    
    echo "RTP: Added $((LINES_AFTER_RTP - LINES_BEFORE_RTP)) lines"
    echo "HDB: Added $((LINES_AFTER_HDB - LINES_BEFORE_HDB)) lines"
    echo "Variants successfully included in CHARMM force field files."
    
    cp "$CHARMM_HDB" "$CHARMM_HDB_BACKUP"
    cp "$CHARMM_RTP" "$CHARMM_RTP_BACKUP"
    echo "Backup files updated with new variants."
else
    echo "WARNING: Variant files not found in $VARIANTS_DIR"
fi

echo "STEP 2 completed successfully!"

echo "========================================"
echo "STEP 3: CARBOHYDRATE ORIENTATION"
echo "========================================"

# ---------------------------------------
# Find required scripts for Step 3
# ---------------------------------------
SCRIPT_PDB_TO_JSON=$(find_python_script "1-pdb_to_json.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_PDB_TO_JSON" ]]; then
    echo "ERROR: Could not find 1-pdb_to_json.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_PDB_TO_JSON"

SCRIPT_ADD_CHARMM=$(find_python_script "3-adding_chamm36_parameters.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_ADD_CHARMM" ]]; then
    echo "ERROR: Could not find 3-adding_chamm36_parameters.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_ADD_CHARMM"

SCRIPT_OPTIMIZE=$(find_python_script "4-optimize_glycans_mcmc.py" "$PYTHON_SCRIPTS_DIR")
if [[ -z "$SCRIPT_OPTIMIZE" ]]; then
    echo "ERROR: Could not find 4-optimize_glycans_mcmc.py in $PYTHON_SCRIPTS_DIR"
    exit 1
fi
echo "Found: $SCRIPT_OPTIMIZE"

# Find the glycoprotein variants directory
GLYCOPROTEIN_STEP2_DIR="$VARIANTS_DIR"
if [[ ! -d "$GLYCOPROTEIN_STEP2_DIR" ]]; then
    echo "WARNING: VALENCE_GLYCAN_VARIANTS directory not found, looking for variants..."
    GLYCOPROTEIN_STEP2_DIR=$(find "$STEP2_OUTPUT" -type d -name "*VARIANTS*" | head -n 1)
    if [[ -z "$GLYCOPROTEIN_STEP2_DIR" ]]; then
        GLYCOPROTEIN_STEP2_DIR="$JSON_DIR"
    fi
fi

# Find the protein PDB file
PROTEIN_STEP2=$(find "$GLYCOPROTEIN_STEP2_DIR" -maxdepth 1 -name "*_variants.pdb" 2>/dev/null | head -n 1)
if [[ -z "$PROTEIN_STEP2" ]]; then
    PROTEIN_STEP2=$(find "$PDB_GLYCO_DIR" -maxdepth 1 -name "*_final_valence_corrected.pdb" 2>/dev/null | head -n 1)
fi
if [[ -z "$PROTEIN_STEP2" ]]; then
    PROTEIN_STEP2=$(find "$PDB_GLYCO_DIR" -maxdepth 1 -name "*.pdb" 2>/dev/null | grep -v "without_H" | head -n 1)
fi

if [[ -z "$PROTEIN_STEP2" ]] || [[ ! -f "$PROTEIN_STEP2" ]]; then
    echo "ERROR: Could not find protein PDB file in $STEP2_OUTPUT"
    exit 1
fi

echo "Found input protein for Step 3: $PROTEIN_STEP2"

# Extract basename for output files
STEP3_BASENAME=$(basename "$PROTEIN_STEP2" .pdb)
STEP3_BASENAME=${STEP3_BASENAME%_variants}
STEP3_BASENAME=${STEP3_BASENAME%_final_valence_corrected}
STEP3_BASENAME=${STEP3_BASENAME%_glycosylated_final_valence_corrected}

echo "Using basename: $STEP3_BASENAME"

# Define output files for Step 3
PDB_TO_JSON_OUTPUT="$JSON_FILES_DIR/pdb_to_json.json"
CHARMM36_JSON_OUTPUT="$JSON_FILES_DIR/glycan_data_charmm36.json"
OPTIMIZED_JSON_OUTPUT="$OPTIMIZED_DIR/glycan_optimized.json"
OPTIMIZED_PDB_OUTPUT="$OPTIMIZED_DIR/${STEP3_BASENAME}_optimized.pdb"
FINAL_REPORT_FILE="$REPORT_FILE"

# ---------------------------------------
# Step 3.1: Convert PDB to JSON
# ---------------------------------------
echo "Step 3.1: Converting PDB to JSON..."

python3 "$SCRIPT_PDB_TO_JSON" \
    --input_pdb "$PROTEIN_STEP2" \
    --output_json "$PDB_TO_JSON_OUTPUT"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 3.1 failed"
    exit 1
fi

# ---------------------------------------
# Step 3.2: Add CHARMM36 parameters
# ---------------------------------------
echo "Step 3.2: Adding CHARMM36 parameters..."

python3 "$SCRIPT_ADD_CHARMM" \
    --input_json "$PDB_TO_JSON_OUTPUT" \
    --charmm_dir "$CHARMM_DIR" \
    --output_json "$CHARMM36_JSON_OUTPUT"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 3.2 failed"
    exit 1
fi

# ---------------------------------------
# Step 3.3: Optimize glycans using MCMC
# ---------------------------------------
echo "Step 3.3: Optimizing glycans using MCMC..."

python3 "$SCRIPT_OPTIMIZE" \
    --input_json "$CHARMM36_JSON_OUTPUT" \
    --output_json "$OPTIMIZED_JSON_OUTPUT" \
    --output_pdb "$OPTIMIZED_PDB_OUTPUT" \
    --glycans_output_dir "$GLYCANS_ONLY_DIR" \
    --theta_step "$THETA_STEP" \
    --n_steps "$N_STEPS" \
    --max_cycles "$MAX_CYCLES" \
    --radius "$RADIUS" \
    --use_coulomb "$USE_COULOMB" \
    --n_workers "$N_WORKERS" \
    --report_file "$FINAL_REPORT_FILE" \
    $SAVE_INDIVIDUAL_GLYCANS \
    $SAVE_BEFORE_AFTER \
    $VERBOSE

if [[ $? -ne 0 ]]; then
    echo "ERROR: Step 3.3 failed"
    exit 1
fi

echo "STEP 3 completed successfully!"

############################################
# COMPLETION
############################################
echo "========================================"
echo "All steps completed successfully!"
echo "========================================"
echo "Results saved in: $OUTPUT_BASE"
echo ""
echo "STEP 1 (Glycosylation Preparation):"
echo "  - TSV files: $TSV_DIR"
echo "  - PDB files: $PDB_DIR"
echo "  - Glycans: $GLYCAN_DIR"
echo ""
echo "STEP 2 (Parametrization):"
echo "  - JSON files: $JSON_DIR"
echo "  - Glycoprotein PDB: $PDB_GLYCO_DIR"
echo "  - Variants: $VARIANTS_DIR"
echo "  - CHARMM36: $CHARMM_DIR"
echo ""
echo "STEP 3 (Carbohydrate Orientation):"
echo "  - JSON files: $JSON_FILES_DIR"
echo "  - Optimized structures: $OPTIMIZED_DIR"
echo "  - Individual glycans: $GLYCANS_ONLY_DIR"
echo ""
echo "Main output files:"
echo "  - Optimized PDB: $OPTIMIZED_PDB_OUTPUT"
echo "  - Report: $FINAL_REPORT_FILE"
echo "========================================"
