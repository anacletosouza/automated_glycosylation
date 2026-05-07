#!/bin/bash

# Define base directories
STEP2_DIR="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/OMICRON/2-GLYCOPROTEIN_TOPOLOGY"
STEP3_RESULTS="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/OMICRON/3-MINIMIZATION_CARBOHYDRATE"

# CHARMM36 directory
CHARMM36="$STEP2_DIR/charmm36.ff"

# JSON directory
JSON="$STEP2_DIR/JSON"

# Glycoprotein directories
GLYCOPROTEIN_STEP2_DIR="$STEP2_DIR/VALENCE_GLYCAN_VARIANTS"
PROTEIN_STEP2="$GLYCOPROTEIN_STEP2_DIR/spike_glycosylated_final_valence_corrected_variants.pdb"

# Scripts directory
SCRIPTS="/grain/anacleto/projects/project_2_automatized_glycosylation_in_glycoproteins/S_GLYCOSYLATION/OMICRON/python_scripts/3-carbohydrate_orientation"

# Create output directories
mkdir -p "$STEP3_RESULTS/JSON_FILES"
mkdir -p "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED"
mkdir -p "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PBD_CARB_ONLY"

# Step 1: Convert PDB to JSON
echo "Step 1: Converting PDB to JSON..."
python3 "$SCRIPTS/1-pdb_to_json.py" \
    --input_pdb "$PROTEIN_STEP2" \
    --output_json "$STEP3_RESULTS/JSON_FILES/pdb_to_json.json"

# Step 3: Add CHARMM36 parameters
echo "Step 3: Adding CHARMM36 parameters..."
python3 "$SCRIPTS/3-adding_chamm36_parameters.py" \
    --input_json "$STEP3_RESULTS/JSON_FILES/pdb_to_json.json" \
    --charmm_dir "$CHARMM36" \
    --output_json "$STEP3_RESULTS/JSON_FILES/glycan_data_charmm36.json"

# Step 4: Optimize glycans using MCMC
echo "Step 4: Optimizing glycans using MCMC..."
python3 "$SCRIPTS/4-optimize_glycans_mcmc.py" \
    --input_json "$STEP3_RESULTS/JSON_FILES/glycan_data_charmm36.json" \
    --output_json "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/glycan_optimized.json" \
    --output_pdb "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/spike_glycosylated_final_optimized.pdb" \
    --glycans_output_dir "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/PBD_CARB_ONLY" \
    --theta_step 10 \
    --n_steps 10 \
    --max_cycles 5 \
    --radius 300 \
    --use_coulomb no \
    --n_workers 12 \
    --report_file "$STEP3_RESULTS/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/report.txt" \
    --save_individual_glycans \
    --save_before_after \
    --verbose

echo "All steps completed successfully!"
