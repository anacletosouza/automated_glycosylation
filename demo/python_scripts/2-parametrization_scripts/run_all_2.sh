#!/bin/bash
set -euo pipefail

############################################
# DEFINITION
############################################

BASE_DIR="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/DELTA"

SCRIPTS_DIR="$BASE_DIR/python_scripts/2-parametrization_scripts"

TOPO_DIR="$BASE_DIR/2-GLYCOPROTEIN_TOPOLOGY"
JSON_DIR="$TOPO_DIR/JSON"
PDB_DIR="$TOPO_DIR/PDB_GLYCOPROTEIN"

mkdir -p "$JSON_DIR"
mkdir -p "$PDB_DIR"
mkdir -p "$BASE_DIR/2-GLYCOPROTEIN_TOPOLOGY/VALENCE_GLYCAN_VARIANTS"

INPUT_PROTEIN="$BASE_DIR/1-GLYCOPROTEIN_PREPARATION/PDB_PROTEIN_GLYCOSYLATED/spike_glycosylated_renumbered.pdb"
OUTPUT_PROTEIN="$PDB_DIR/spike_glycosylated_corrected.pdb"
FINAL_STRUCTURE="$PDB_DIR/spike_glycosylated_final_connected.pdb"
FINAL_STRUCTURE_2="$PDB_DIR/spike_glycosylated_final_valence_corrected.pdb"

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
       "$HOME/programs/GROMACS/share/gromacs/top/residuetypes.dat"
fi

# Restauração já feita acima, comentada para evitar duplicação
# if [ -f "$CHARMM_HDB_BACKUP" ]; then
#     cp "$CHARMM_HDB_BACKUP" "$CHARMM_HDB"
# fi

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
    --base_dir "$BASE_DIR/1-GLYCOPROTEIN_PREPARATION" \
    --output_dir "$JSON_DIR"

############################################
# PROCESSING
############################################

echo "Starting processing for each glycan directory..."

for DIR in "$JSON_DIR"/*; do
    [ -d "$DIR" ] || continue

    BASENAME=$(basename "$DIR")

    PDB_FILE="$DIR/$BASENAME.pdb"
    JSON_FILE="$DIR/$BASENAME.json"
    PARSER_FILE="$DIR/${BASENAME}_parser.pkl"
    RTP_PICKLE="$DIR/carb_residues.pkl"
    RTP_MODIFIED="$DIR/carb_modified.rtp"
    RTP_UNIQUE="$DIR/carb_unique.rtp"

    echo "----------------------------------------"
    echo "Processing: $BASENAME"

    cd "$DIR"

    python3 "$SCRIPTS_DIR/1-parser_pdb.py" "$PDB_FILE" -o "$PARSER_FILE"

    python3 "$SCRIPTS_DIR/2-parser_carb_rtp.py" "$CHARMM_RTP" -o "$RTP_PICKLE"

    python3 "$SCRIPTS_DIR/3-comparison_pdb_rtp.py" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE"

    python3 "$SCRIPTS_DIR/4-rtp_generator_part1.py" \
        -p "$PARSER_FILE" \
        -r "$RTP_PICKLE"

    python3 "$SCRIPTS_DIR/4-rtp_generator_part2.py" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" \
        --json "$JSON_FILE"

    python3 "$SCRIPTS_DIR/4-rtp_generator_part3.py" \
        --pdb "$PARSER_FILE" \
        --rtp "$RTP_PICKLE" \
        --json "$JSON_FILE" \
        --output "$RTP_MODIFIED"

    python3 "$SCRIPTS_DIR/5-acetylation_replacement.py" \
        "$PDB_FILE" \
        "${DIR}/${BASENAME}_modified.pdb"

    DIR_LETTER="${BASENAME:0:1}"

    python3 "$SCRIPTS_DIR/6-clean_rtp.py" \
        "$RTP_MODIFIED" \
        "$RTP_UNIQUE" \
        "$DIR_LETTER"

    cd - > /dev/null

    echo "Finished $BASENAME"
done

############################################
# UNIFICATION OF RTP / HDB
############################################

python3 "$SCRIPTS_DIR/7-together_part_1.py" \
    --input "$JSON_DIR" \
    --output "$JSON_DIR/carb_unique_total.rtp"

python3 "$SCRIPTS_DIR/7-together_part_2.py" \
    --input "$JSON_DIR/carb_unique_total.rtp" \
    --output "$JSON_DIR/carb_redundance_removed.rtp"

python3 "$SCRIPTS_DIR/7-together_part_3.py" \
    --input "$JSON_DIR/carb_redundance_removed.rtp" \
    --output "$CHARMM_RTP"

python3 "$SCRIPTS_DIR/7-together_part_4.py" \
    --input "$JSON_DIR/carb_redundance_removed.rtp"

python3 "$SCRIPTS_DIR/7-together_part_5.py" \
    --input "$JSON_DIR/carb_redundance_removed.rtp" \
    --output "$JSON_DIR/carb_redundance_removed.hdb"

############################################
# FINAL GLYCOPROTEIN
############################################

python3 "$SCRIPTS_DIR/8-glycoprotein.py" \
    --protein "$INPUT_PROTEIN" \
    --carbs_dir "$JSON_DIR" \
    --output "$OUTPUT_PROTEIN" \
    --keep_hydrogens_carb \
    --keep_hydrogens_prot

python3 "$SCRIPTS_DIR/9-conection_glycosilation_without_TER.py" \
    --glycosylated "$OUTPUT_PROTEIN" \
    --conect "$INPUT_PROTEIN" \
    --output "$FINAL_STRUCTURE"
    
python3 "$SCRIPTS_DIR/glycosylation_identifying.py" "$FINAL_STRUCTURE" "$FINAL_STRUCTURE_2"

############################################
# HDB UPDATE
############################################

python3 "$SCRIPTS_DIR/10-generation_hdb.py" \
    "$JSON_DIR/carb_redundance_removed.rtp" \
    "$CHARMM_HDB" \
    -o "$JSON_DIR/carb_modified.hdb"

python3 "$SCRIPTS_DIR/glycosylation_variants.py" \
        -p "$FINAL_STRUCTURE_2" -r "$JSON_DIR/carb_redundance_removed.rtp" \
        -d "$JSON_DIR/carb_redundance_removed.hdb" \
        -o "$BASE_DIR/2-GLYCOPROTEIN_TOPOLOGY/VALENCE_GLYCAN_VARIANTS"

# (2) MODIFICAÇÃO: Incluir as variantes geradas nos arquivos CHARMM RTP e HDB
echo "Including generated variants in CHARMM force field files..."

VARIANT_RTP="$BASE_DIR/2-GLYCOPROTEIN_TOPOLOGY/VALENCE_GLYCAN_VARIANTS/spike_glycosylated_final_valence_corrected_variants.rtp"
VARIANT_HDB="$BASE_DIR/2-GLYCOPROTEIN_TOPOLOGY/VALENCE_GLYCAN_VARIANTS/spike_glycosylated_final_valence_corrected_variants.hdb"

# Verificar se os arquivos de variantes existem
if [ -f "$VARIANT_RTP" ] && [ -f "$VARIANT_HDB" ]; then
    echo "Appending variant RTP to $CHARMM_RTP..."
    cat "$VARIANT_RTP" >> "$CHARMM_RTP"
    
    echo "Appending variant HDB to $CHARMM_HDB..."
    cat "$VARIANT_HDB" >> "$CHARMM_HDB"
    
    echo "Variants successfully included in CHARMM force field files."
else
    echo "WARNING: Variant files not found:"
    echo "  RTP: $VARIANT_RTP"
    echo "  HDB: $VARIANT_HDB"
    echo "Skipping inclusion of variants in CHARMM files."
fi

# Atualizar backups com os novos arquivos
cp "$CHARMM_HDB" "$CHARMM_HDB_BACKUP"
cp "$CHARMM_RTP" "$CHARMM_RTP_BACKUP"

echo "Backup files updated with new variants."

############################################
# CONCLUSION
############################################

echo "The process was finished successfully!"
