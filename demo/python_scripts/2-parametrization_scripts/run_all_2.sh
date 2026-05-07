#!/bin/bash
set -euo pipefail

############################################
# USAGE FUNCTION
############################################
usage() {
    echo "Usage: $0 --base_dir <input_base_directory> --output_dir <output_directory>"
    echo "  --base_dir   : Directory containing output from previous pipeline (e.g., STEP_1/)"
    echo "  --output_dir : Output directory for topology files (e.g., STEP_2/)"
    exit 1
}

############################################
# PARSE COMMAND LINE ARGUMENTS
############################################
BASE_DIR=""
OUTPUT_BASE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --base_dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --output_dir)
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
if [[ -z "$BASE_DIR" ]] || [[ -z "$OUTPUT_BASE" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Convert to absolute paths
BASE_DIR=$(realpath "$BASE_DIR")
OUTPUT_BASE=$(realpath "$OUTPUT_BASE")

# Check if base directory exists
if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: Base directory not found: $BASE_DIR"
    exit 1
fi

############################################
# DEFINITION OF PATHS
############################################

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR"

# Output directories
TOPO_DIR="$OUTPUT_BASE"
JSON_DIR="$TOPO_DIR/JSON"
PDB_DIR="$TOPO_DIR/PDB_GLYCOPROTEIN"
VARIANTS_DIR="$TOPO_DIR/VALENCE_GLYCAN_VARIANTS"

mkdir -p "$JSON_DIR"
mkdir -p "$PDB_DIR"
mkdir -p "$VARIANTS_DIR"

# Find input files from previous pipeline
PREP_DIR="$BASE_DIR"
PDB_PREP_DIR="$PREP_DIR/PDB_PROTEIN_GLYCOSYLATED"

# Find the glycosylated renumbered PDB file (generic)
INPUT_PROTEIN=$(find "$PDB_PREP_DIR" -name "*_glycosylated_renumbered.pdb" | head -n 1)
if [[ -z "$INPUT_PROTEIN" ]]; then
    echo "ERROR: Could not find *_glycosylated_renumbered.pdb in $PDB_PREP_DIR"
    exit 1
fi

echo "Found input protein: $INPUT_PROTEIN"

# Extract basename for output files
BASENAME=$(basename "$INPUT_PROTEIN" _glycosylated_renumbered.pdb)
echo "Base filename: $BASENAME"

OUTPUT_PROTEIN="$PDB_DIR/${BASENAME}_glycosylated_corrected.pdb"
FINAL_STRUCTURE="$PDB_DIR/${BASENAME}_glycosylated_final_connected.pdb"
FINAL_STRUCTURE_2="$PDB_DIR/${BASENAME}_glycosylated_final_valence_corrected.pdb"
FINAL_STRUCTURE_2_NOH="$PDB_DIR/${BASENAME}_glycosylated_final_valence_corrected_noh.pdb"

############################################
# CHARMM FILES
############################################

CHARMM_DIR="$TOPO_DIR/charmm36.ff"
CHARMM_RTP="$CHARMM_DIR/carb.rtp"
CHARMM_HDB="$CHARMM_DIR/carb.hdb"
CHARMM_HDB_BACKUP="$CHARMM_DIR/carb.hdb.backup"
CHARMM_RTP_BACKUP="$CHARMM_DIR/carb.rtp.backup"

############################################
# PREPARATION OF FORCE FIELD
############################################

cd "$TOPO_DIR"

# (1) MODIFICAÇÃO: Verificar se já existe o backup e pular download se existir
if [ ! -f "$CHARMM_HDB_BACKUP" ] || [ ! -f "$CHARMM_RTP_BACKUP" ]; then
    echo "Backup files not found. Downloading CHARMM force field..."
    
    rm -rf "$CHARMM_DIR"
    
    wget https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz
    
    tar -xzf download.php\?filename\=CHARMM_ff_params_files%2Fcharmm36-jul2022.ff.tgz
    
    mv charmm36-jul2022.ff "$CHARMM_DIR"
    
    # Create initial backups
    cp "$CHARMM_HDB" "$CHARMM_HDB_BACKUP"
    cp "$CHARMM_RTP" "$CHARMM_RTP_BACKUP"
    
    echo "CHARMM force field downloaded and backups created."
else
    echo "Backup files found. Restoring from backups..."
    
    # Create directory if it doesn't exist
    mkdir -p "$CHARMM_DIR"
    
    # Restore from backups
    cp "$CHARMM_HDB_BACKUP" "$CHARMM_HDB"
    cp "$CHARMM_RTP_BACKUP" "$CHARMM_RTP"
    
    echo "Force field restored from backups."
fi

############################################
# BACKUPS RESTORATION
############################################

if [ -f "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat.backup" ]; then
    cp "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat.backup" \
       "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat" || true
fi

############################################
# CHECKING
############################################

for FILE in "$CHARMM_RTP" "$CHARMM_HDB"; do
    if [ ! -f "$FILE" ]; then
        echo "ERROR: File not found: $FILE"
        exit 1
    fi
done

############################################
# GENERATION OF THE JSONs
############################################

echo "Generating JSON for each glycan..."

python3 "$SCRIPTS_DIR/0-JSON_generator.py" \
    --base_dir "$BASE_DIR" \
    --output_dir "$JSON_DIR"

############################################
# PROCESSING
############################################

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

    # Skip if PDB file doesn't exist
    if [[ ! -f "$PDB_FILE" ]]; then
        echo "WARNING: PDB file not found: $PDB_FILE, skipping..."
        continue
    fi

    echo "----------------------------------------"
    echo "Processing: $DIR_BASENAME"

    cd "$DIR"

    python3 "$SCRIPTS_DIR/1-parser_pdb.py" "$PDB_FILE" -o "$PARSER_FILE" || {
        echo "ERROR in 1-parser_pdb.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPTS_DIR/2-parser_carb_rtp.py" "$CHARMM_RTP" -o "$RTP_PICKLE" || {
        echo "ERROR in 2-parser_carb_rtp.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPTS_DIR/3-comparison_pdb_rtp.py" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" || {
        echo "ERROR in 3-comparison_pdb_rtp.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPTS_DIR/4-rtp_generator_part1.py" \
        -p "$PARSER_FILE" \
        -r "$RTP_PICKLE" || {
        echo "ERROR in 4-rtp_generator_part1.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPTS_DIR/4-rtp_generator_part2.py" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" \
        --json "$JSON_FILE" || {
        echo "ERROR in 4-rtp_generator_part2.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPTS_DIR/4-rtp_generator_part3.py" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" \
        --json "$JSON_FILE" \
        --output "$RTP_MODIFIED" || {
        echo "ERROR in 4-rtp_generator_part3.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    python3 "$SCRIPTS_DIR/5-acetylation_replacement.py" \
        "$PDB_FILE" \
        "${DIR}/${DIR_BASENAME}_modified.pdb" || {
        echo "ERROR in 5-acetylation_replacement.py for $DIR_BASENAME, skipping..."
        cd - > /dev/null
        continue
    }

    DIR_LETTER="${DIR_BASENAME:0:1}"

    python3 "$SCRIPTS_DIR/6-clean_rtp.py" \
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

############################################
# UNIFICATION OF RTP / HDB
############################################

echo "Unifying RTP files..."

# Check if any RTP files were generated
if [ -z "$(ls -A "$JSON_DIR"/*/carb_unique.rtp 2>/dev/null)" ]; then
    echo "WARNING: No carb_unique.rtp files found. Skipping unification steps."
else
    python3 "$SCRIPTS_DIR/7-together_part_1.py" \
        --input "$JSON_DIR" \
        --output "$JSON_DIR/carb_unique_total.rtp" || {
        echo "ERROR in 7-together_part_1.py"
        exit 1
    }

    python3 "$SCRIPTS_DIR/7-together_part_2.py" \
        --input "$JSON_DIR/carb_unique_total.rtp" \
        --output "$JSON_DIR/carb_redundance_removed.rtp" || {
        echo "ERROR in 7-together_part_2.py"
        exit 1
    }

    python3 "$SCRIPTS_DIR/7-together_part_3.py" \
        --input "$JSON_DIR/carb_redundance_removed.rtp" \
        --output "$CHARMM_RTP" || {
        echo "ERROR in 7-together_part_3.py"
        exit 1
    }

    python3 "$SCRIPTS_DIR/7-together_part_4.py" \
        --input "$JSON_DIR/carb_redundance_removed.rtp" || {
        echo "ERROR in 7-together_part_4.py"
        exit 1
    }

    python3 "$SCRIPTS_DIR/7-together_part_5.py" \
        --input "$JSON_DIR/carb_redundance_removed.rtp" \
        --output "$JSON_DIR/carb_redundance_removed.hdb" || {
        echo "ERROR in 7-together_part_5.py"
        exit 1
    }
fi

############################################
# FINAL GLYCOPROTEIN
############################################

echo "Building final glycoprotein..."

# Check if JSON_DIR has any subdirectories with necessary files
if [ -z "$(ls -d "$JSON_DIR"/*/ 2>/dev/null)" ]; then
    echo "WARNING: No glycan directories found in $JSON_DIR"
    echo "Skipping glycoprotein construction steps."
else
    python3 "$SCRIPTS_DIR/8-glycoprotein.py" \
        --protein "$INPUT_PROTEIN" \
        --carbs_dir "$JSON_DIR" \
        --output "$OUTPUT_PROTEIN" \
        --keep_hydrogens_carb \
        --keep_hydrogens_prot || {
        echo "ERROR in 8-glycoprotein.py"
        exit 1
    }

    python3 "$SCRIPTS_DIR/9-conection_glycosilation_without_TER.py" \
        --glycosylated "$OUTPUT_PROTEIN" \
        --conect "$INPUT_PROTEIN" \
        --output "$FINAL_STRUCTURE" || {
        echo "ERROR in 9-conection_glycosilation_without_TER.py"
        exit 1
    }
    
    python3 "$SCRIPTS_DIR/glycosylation_identifying.py" "$FINAL_STRUCTURE" "$FINAL_STRUCTURE_2" || {
        echo "ERROR in glycosylation_identifying.py"
        exit 1
    }

    # Check if the expected output file exists, otherwise look for alternative
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
fi

############################################
# HDB UPDATE
############################################

echo "Updating HDB files..."

if [[ -f "$JSON_DIR/carb_redundance_removed.rtp" ]]; then
    python3 "$SCRIPTS_DIR/10-generation_hdb.py" \
        "$JSON_DIR/carb_redundance_removed.rtp" \
        "$CHARMM_HDB" \
        -o "$JSON_DIR/carb_modified.hdb" || {
        echo "ERROR in 10-generation_hdb.py"
        exit 1
    }
else
    echo "WARNING: carb_redundance_removed.rtp not found, skipping HDB generation"
fi

############################################
# VARIANTS GENERATION
############################################

echo "Generating glycan variants..."

if [[ -f "$FINAL_STRUCTURE_2" ]] && [[ -f "$JSON_DIR/carb_redundance_removed.rtp" ]] && [[ -f "$JSON_DIR/carb_redundance_removed.hdb" ]]; then
    python3 "$SCRIPTS_DIR/glycosylation_variants.py" \
        -p "$FINAL_STRUCTURE_2" \
        -r "$JSON_DIR/carb_redundance_removed.rtp" \
        -d "$JSON_DIR/carb_redundance_removed.hdb" \
        -o "$VARIANTS_DIR" || {
        echo "WARNING: glycosylation_variants.py failed, but continuing..."
    }
else
    echo "WARNING: Required files for variant generation not found"
    echo "  - FINAL_STRUCTURE_2: $FINAL_STRUCTURE_2"
    echo "  - RTP file: $JSON_DIR/carb_redundance_removed.rtp"
    echo "  - HDB file: $JSON_DIR/carb_redundance_removed.hdb"
fi

############################################
# INCLUDE VARIANTS IN CHARMM FILES
############################################

echo "Including generated variants in CHARMM force field files..."

# Find variant files generically with multiple patterns
VARIANT_RTP=$(find "$VARIANTS_DIR" -name "*_variants.rtp" 2>/dev/null | head -n 1)
if [[ -z "$VARIANT_RTP" ]]; then
    VARIANT_RTP=$(find "$VARIANTS_DIR" -name "*.rtp" 2>/dev/null | grep -v "carb_" | head -n 1)
fi

VARIANT_HDB=$(find "$VARIANTS_DIR" -name "*_variants.hdb" 2>/dev/null | head -n 1)
if [[ -z "$VARIANT_HDB" ]]; then
    VARIANT_HDB=$(find "$VARIANTS_DIR" -name "*.hdb" 2>/dev/null | grep -v "carb_" | head -n 1)
fi

# Verificar se os arquivos de variantes existem
if [[ -f "$VARIANT_RTP" ]] && [[ -f "$VARIANT_HDB" ]]; then
    echo "Found variant RTP: $VARIANT_RTP"
    echo "Found variant HDB: $VARIANT_HDB"
    
    # Count lines before appending
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
    
    # Atualizar backups com os novos arquivos
    cp "$CHARMM_HDB" "$CHARMM_HDB_BACKUP"
    cp "$CHARMM_RTP" "$CHARMM_RTP_BACKUP"
    echo "Backup files updated with new variants."
else
    echo "WARNING: Variant files not found in $VARIANTS_DIR"
    echo "  RTP found: ${VARIANT_RTP:-none}"
    echo "  HDB found: ${VARIANT_HDB:-none}"
    echo "Skipping inclusion of variants in CHARMM files."
fi

############################################
# CONCLUSION
############################################

echo "========================================"
echo "The process was finished successfully!"
echo "Output directory: $TOPO_DIR"
echo "JSON directory: $JSON_DIR"
echo "PDB directory: $PDB_DIR"
echo "Variants directory: $VARIANTS_DIR"
echo "========================================"
