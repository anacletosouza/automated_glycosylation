#!/bin/bash
# -----------------------------------------------------------------------------
# run_all.sh - Complete pipeline runner for glycosylation workflow
# 
# This script runs the full glycosylation pipeline with all options
# It can run individual steps or the complete workflow
# -----------------------------------------------------------------------------

set -euo pipefail

# Default values
INPUT_PDB=""
INPUT_TSV=""
INPUT_GLYCOSYLATOR_TSV=""
PREP_OUTPUT_DIR=""
PARAM_OUTPUT_DIR=""
ORIENT_OUTPUT_DIR=""
CHARMM_DIR=""

# Prep parameters
PROTEIN_RESIDUE_START=1
ROTATE_ATOMS="OD1,CG,ND2,HD22,HD21,HB2,HB3"
FIXED_ATOM="CB"
CENTER_ATOM="CA"
RADIUS_PREP=30.0
ROTATION_STEP=1

# Param parameters
SKIP_CHARMM_DOWNLOAD=false

# Orient parameters
THETA_STEP=10
N_STEPS=10
MAX_CYCLES=5
RADIUS_ORIENT=300
USE_COULOMB="no"
N_WORKERS=1
SAVE_INDIVIDUAL_GLYCANS=false
SAVE_BEFORE_AFTER=false
VERBOSE=false

# Step selection
RUN_PREP=true
RUN_PARAM=true
RUN_ORIENT=true

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Complete pipeline runner for glycosylation workflow

REQUIRED OPTIONS:
    --input-pdb PATH            Input PDB file
    --prep-output-dir PATH      Output directory for preparation step
    --param-output-dir PATH     Output directory for parametrization step
    --orient-output-dir PATH    Output directory for orientation step

INPUT OPTIONS:
    --input-tsv PATH            Input TSV file (Caselino format)
    --input-glycosylator-tsv PATH  Pre-processed glycosylator TSV

PREPARATION OPTIONS:
    --protein-residue-start NUM Protein residue start number (default: 1)
    --rotate-atoms LIST         Atoms to rotate (default: OD1,CG,ND2,HD22,HD21,HB2,HB3)
    --fixed-atom ATOM           Fixed atom (default: CB)
    --center-atom ATOM          Center atom (default: CA)
    --radius-prep NUM           Radius for orientation (default: 30.0)
    --rotation-step NUM         Rotation step in degrees (default: 1)

PARAMETRIZATION OPTIONS:
    --skip-charmm-download      Skip CHARMM force field download

ORIENTATION OPTIONS:
    --charmm-dir PATH           CHARMM36 directory (optional)
    --theta-step NUM            Theta step for MCMC (default: 10)
    --n-steps NUM               Number of steps for MCMC (default: 10)
    --max-cycles NUM            Maximum cycles for MCMC (default: 5)
    --radius-orient NUM         Radius for orientation (default: 300)
    --use-coulomb yes|no        Use Coulomb potential (default: no)
    --n-workers NUM             Number of workers (default: 1)
    --save-individual-glycans   Save individual glycans
    --save-before-after         Save before/after structures
    --verbose                   Verbose output

STEP SELECTION:
    --prep-only                 Run only preparation step
    --param-only                Run only parametrization step
    --orient-only               Run only orientation step

OTHER:
    -h, --help                  Show this help message

EXAMPLES:
    # Run complete pipeline
    $0 --input-pdb protein.pdb \\
        --prep-output-dir ./prep \\
        --param-output-dir ./param \\
        --orient-output-dir ./orient

    # Run with TSV file and custom parameters
    $0 --input-pdb protein.pdb \\
        --input-tsv table.tsv \\
        --prep-output-dir ./prep \\
        --param-output-dir ./param \\
        --orient-output-dir ./orient \\
        --protein-residue-start 10 \\
        --theta-step 20 \\
        --n-steps 50

    # Run only preparation
    $0 --prep-only --input-pdb protein.pdb --prep-output-dir ./prep

    # Run only orientation with custom settings
    $0 --orient-only \\
        --input-pdb protein_param.pdb \\
        --param-output-dir ./param \\
        --orient-output-dir ./orient \\
        --verbose
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input-pdb)
            INPUT_PDB="$2"
            shift 2
            ;;
        --input-tsv)
            INPUT_TSV="$2"
            shift 2
            ;;
        --input-glycosylator-tsv)
            INPUT_GLYCOSYLATOR_TSV="$2"
            shift 2
            ;;
        --prep-output-dir)
            PREP_OUTPUT_DIR="$2"
            shift 2
            ;;
        --param-output-dir)
            PARAM_OUTPUT_DIR="$2"
            shift 2
            ;;
        --orient-output-dir)
            ORIENT_OUTPUT_DIR="$2"
            shift 2
            ;;
        --charmm-dir)
            CHARMM_DIR="$2"
            shift 2
            ;;
        --protein-residue-start)
            PROTEIN_RESIDUE_START="$2"
            shift 2
            ;;
        --rotate-atoms)
            ROTATE_ATOMS="$2"
            shift 2
            ;;
        --fixed-atom)
            FIXED_ATOM="$2"
            shift 2
            ;;
        --center-atom)
            CENTER_ATOM="$2"
            shift 2
            ;;
        --radius-prep)
            RADIUS_PREP="$2"
            shift 2
            ;;
        --rotation-step)
            ROTATION_STEP="$2"
            shift 2
            ;;
        --skip-charmm-download)
            SKIP_CHARMM_DOWNLOAD=true
            shift
            ;;
        --theta-step)
            THETA_STEP="$2"
            shift 2
            ;;
        --n-steps)
            N_STEPS="$2"
            shift 2
            ;;
        --max-cycles)
            MAX_CYCLES="$2"
            shift 2
            ;;
        --radius-orient)
            RADIUS_ORIENT="$2"
            shift 2
            ;;
        --use-coulomb)
            USE_COULOMB="$2"
            shift 2
            ;;
        --n-workers)
            N_WORKERS="$2"
            shift 2
            ;;
        --save-individual-glycans)
            SAVE_INDIVIDUAL_GLYCANS=true
            shift
            ;;
        --save-before-after)
            SAVE_BEFORE_AFTER=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --prep-only)
            RUN_PREP=true
            RUN_PARAM=false
            RUN_ORIENT=false
            shift
            ;;
        --param-only)
            RUN_PREP=false
            RUN_PARAM=true
            RUN_ORIENT=false
            shift
            ;;
        --orient-only)
            RUN_PREP=false
            RUN_PARAM=false
            RUN_ORIENT=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate required arguments
validate_requirements() {
    if [ "$RUN_PREP" = true ] || [ "$RUN_PARAM" = true ] || [ "$RUN_ORIENT" = true ]; then
        if [ -z "$INPUT_PDB" ]; then
            echo "ERROR: --input-pdb is required"
            exit 1
        fi
    fi
    
    if [ "$RUN_PREP" = true ] && [ -z "$PREP_OUTPUT_DIR" ]; then
        echo "ERROR: --prep-output-dir is required for preparation step"
        exit 1
    fi
    
    if [ "$RUN_PARAM" = true ] && [ -z "$PARAM_OUTPUT_DIR" ]; then
        echo "ERROR: --param-output-dir is required for parametrization step"
        exit 1
    fi
    
    if [ "$RUN_ORIENT" = true ] && [ -z "$ORIENT_OUTPUT_DIR" ]; then
        echo "ERROR: --orient-output-dir is required for orientation step"
        exit 1
    fi
}

# Build and run preparation command
run_preparation() {
    echo ""
    echo "============================================================"
    echo "RUNNING PREPARATION STEP"
    echo "============================================================"
    
    local cmd="glyco-prep"
    cmd="$cmd --input-pdb $INPUT_PDB"
    cmd="$cmd --output-dir $PREP_OUTPUT_DIR"
    
    [ -n "$INPUT_TSV" ] && cmd="$cmd --input-tsv $INPUT_TSV"
    [ -n "$INPUT_GLYCOSYLATOR_TSV" ] && cmd="$cmd --input-glycosylator-tsv $INPUT_GLYCOSYLATOR_TSV"
    
    cmd="$cmd --protein-residue-start $PROTEIN_RESIDUE_START"
    cmd="$cmd --rotate-atoms \"$ROTATE_ATOMS\""
    cmd="$cmd --fixed-atom $FIXED_ATOM"
    cmd="$cmd --center-atom $CENTER_ATOM"
    cmd="$cmd --radius $RADIUS_PREP"
    cmd="$cmd --rotation-step $ROTATION_STEP"
    
    echo "Running: $cmd"
    eval $cmd
    
    # Get the final PDB from preparation
    FINAL_PREP_PDB="$PREP_OUTPUT_DIR/PDB_PROTEIN_GLYCOSYLATED/protein_renumbered.pdb"
    if [ ! -f "$FINAL_PREP_PDB" ]; then
        FINAL_PREP_PDB="$PREP_OUTPUT_DIR/PDB_PROTEIN_GLYCOSYLATED/protein_asn_orientation.pdb"
    fi
    
    echo "Preparation completed. Final PDB: $FINAL_PREP_PDB"
}

# Build and run parametrization command
run_parametrization() {
    echo ""
    echo "============================================================"
    echo "RUNNING PARAMETRIZATION STEP"
    echo "============================================================"
    
    local input_pdb="${1:-$FINAL_PREP_PDB}"
    
    local cmd="glyco-param"
    cmd="$cmd --prep-output-dir $PREP_OUTPUT_DIR"
    cmd="$cmd --input-pdb $input_pdb"
    cmd="$cmd --output-dir $PARAM_OUTPUT_DIR"
    
    if [ "$SKIP_CHARMM_DOWNLOAD" = true ]; then
        cmd="$cmd --skip-charmm-download"
    fi
    
    echo "Running: $cmd"
    eval $cmd
    
    # Get the final PDB from parametrization
    FINAL_PARAM_PDB="$PARAM_OUTPUT_DIR/PDB_GLYCOPROTEIN/protein_final_valence_corrected.pdb"
    
    echo "Parametrization completed. Final PDB: $FINAL_PARAM_PDB"
}

# Build and run orientation command
run_orientation() {
    echo ""
    echo "============================================================"
    echo "RUNNING ORIENTATION STEP"
    echo "============================================================"
    
    local input_pdb="${1:-$FINAL_PARAM_PDB}"
    
    local cmd="glyco-orient"
    cmd="$cmd --input-pdb $input_pdb"
    cmd="$cmd --param-output-dir $PARAM_OUTPUT_DIR"
    cmd="$cmd --output-dir $ORIENT_OUTPUT_DIR"
    
    [ -n "$CHARMM_DIR" ] && cmd="$cmd --charmm-dir $CHARMM_DIR"
    
    cmd="$cmd --theta-step $THETA_STEP"
    cmd="$cmd --n-steps $N_STEPS"
    cmd="$cmd --max-cycles $MAX_CYCLES"
    cmd="$cmd --radius $RADIUS_ORIENT"
    cmd="$cmd --use-coulomb $USE_COULOMB"
    cmd="$cmd --n-workers $N_WORKERS"
    
    if [ "$SAVE_INDIVIDUAL_GLYCANS" = true ]; then
        cmd="$cmd --save-individual-glycans"
    fi
    
    if [ "$SAVE_BEFORE_AFTER" = true ]; then
        cmd="$cmd --save-before-after"
    fi
    
    if [ "$VERBOSE" = true ]; then
        cmd="$cmd --verbose"
    fi
    
    echo "Running: $cmd"
    eval $cmd
    
    # Get the final PDB from orientation
    FINAL_ORIENT_PDB="$ORIENT_OUTPUT_DIR/PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/protein_optimized.pdb"
    
    echo "Orientation completed. Final PDB: $FINAL_ORIENT_PDB"
}

# Main execution
main() {
    echo ""
    echo "============================================================"
    echo "AUTOMATED GLYCOSYLATION PIPELINE"
    echo "============================================================"
    echo "Start time: $(date)"
    echo ""
    
    validate_requirements
    
    # Track final PDBs
    local final_prep=""
    local final_param=""
    local final_orient=""
    
    # Run steps based on selection
    if [ "$RUN_PREP" = true ]; then
        run_preparation
        final_prep="$FINAL_PREP_PDB"
    fi
    
    if [ "$RUN_PARAM" = true ]; then
        local input_for_param="${final_prep:-$INPUT_PDB}"
        run_parametrization "$input_for_param"
        final_param="$FINAL_PARAM_PDB"
    fi
    
    if [ "$RUN_ORIENT" = true ]; then
        local input_for_orient="${final_param:-$INPUT_PDB}"
        run_orientation "$input_for_orient"
        final_orient="$FINAL_ORIENT_PDB"
    fi
    
    # Summary
    echo ""
    echo "============================================================"
    echo "PIPELINE COMPLETED SUCCESSFULLY!"
    echo "============================================================"
    echo "End time: $(date)"
    echo ""
    echo "Results:"
    [ -n "$final_prep" ] && echo "  Preparation:    $final_prep"
    [ -n "$final_param" ] && echo "  Parametrization: $final_param"
    [ -n "$final_orient" ] && echo "  Orientation:    $final_orient"
    echo ""
}

# Run main function
main
