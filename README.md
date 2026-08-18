# Auto_glyco: Automated glycosylation pipeline for glycoprotein modeling in all-atom model resolution.

## Observation: This tool is in test phase.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Development Status](https://img.shields.io/badge/status-beta-orange.svg)](https://github.com/anacletosouza/automated_glycosylation)

A comprehensive pipeline for automated glycosylation of proteins, including glycan attachment, parametrization, and carbohydrate orientation optimization using MCMC methods.

## Overview

The Automated Glycosylation Pipeline provides a complete workflow for:

1. **Glycosylation Preparation** - Attach glycans to protein structures based on TSV input
2. **Parametrization** - Generate CHARMM36 force field parameters for glycans
3. **Carbohydrate Orientation** - Optimize glycan orientations using Monte Carlo methods

## Features

- **Automated glycan attachment** from IUPAC sequences
- **CHARMM36 force field** parameter generation
- **MCMC-based orientation optimization** for realistic glycan positioning
- **Parallel processing** support for large glycoproteins
- **Comprehensive reporting** of clashes and optimization results
- **PDB output** with optimized glycoprotein structures

## Installation

### Prerequisites

- Python 3.8 or higher
- pip (Python package installer)
- wget (for CHARMM force field download)

### Install from GitHub

```bash
# Clone the repository
git clone https://github.com/anacletosouza/automated_glycosylation.git
cd automated_glycosylation

# Install in development mode (recommended for development)
pip install -e .

# Or install normally
pip install .
```

### Install with Conda (alternative)

```bash
# Create a new conda environment
conda create -n glycosylation python=3.9
conda activate glycosylation

# Install the package
git clone https://github.com/anacletosouza/automated_glycosylation.git
cd automated_glycosylation
pip install -e .
```

## Quick Start

### Basic Usage

In this example, it is necessary to determine the protonation state of the protein before running the code

```bash
auto_glyco --pdb protein.pdb --tsv input/glycans.tsv
```

### Complete Example

```bash
auto_glyco \
    --pdb examples/test.pdb \
    --tsv examples/test.tsv \
    --output_dir OUTPUT \
    --theta_step 5 \
    --n_steps 20 \
    --max_cycles 10 \
    --radius 200 \
    --n_workers 8 \
    --save_individual_glycans \
    --verbose
```

## Input Files

### PDB File
A standard PDB file containing your protein structure. The file should be properly formatted with correct atom coordinates.

Example:
```pdb
ATOM      1  N   GLY A   1      12.345  67.890  23.456  1.00  0.00           N
ATOM      2  CA  GLY A   1      13.456  68.901  24.567  1.00  0.00           C
...
```

### TSV File
A tab-separated values file defining glycosylation sites. Each row represents a glycan to be attached.

Format:
```
Residue_Number	Chain	Glycan_IUPAC	Residue_Name
45	A	GlcNAc	ASN
78	B	Man(GlcNAc)2	ASN
123	A	GlcNAc(Man)3	ASN
```

- **Residue_Number**: The residue number in the PDB file
- **Chain**: Chain identifier (A, B, C, etc.)
- **Glycan_IUPAC**: IUPAC notation of the glycan structure
- **Residue_Name**: Target residue name (typically ASN for N-linked glycosylation)

## Command Line Arguments

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--pdb` | Input PDB file path | `--pdb protein.pdb` |
| `--tsv` | Input TSV file path | `--tsv sites.tsv` |

### Optional Arguments

#### General Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--output_dir` | `./output` | Output directory for all results |
| `--url_charmm36` | CHARMM36 URL | URL for CHARMM36 force field download |
| `--protein_residue_start` | `10` | Starting residue number for protein |
| `--verbose` | `False` | Enable verbose output |

#### Orientation Optimization Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--theta_step` | `10` | Rotation step in degrees |
| `--n_steps` | `10` | Number of rotation steps per cycle |
| `--max_cycles` | `5` | Maximum MCMC cycles |
| `--radius` | `300` | Radius (in Å) for clash detection |
| `--use_coulomb` | `no` | Use Coulomb interactions (`yes`/`no`) |
| `--n_workers` | `12` | Number of parallel workers |
| `--report_file` | Auto-generated | Path to report file |

#### Output Options

| Argument | Description |
|----------|-------------|
| `--save_individual_glycans` | Save individual glycan PDB files |
| `--save_before_after` | Save before/after comparison files |

## Output Structure

```
output/
├── STEP1/                          # Glycosylation preparation
│   ├── TSV/                        # Processed TSV files
│   │   ├── input_corrected.tsv
│   │   └── input_glycosylator.tsv
│   ├── PDB_PROTEIN_GLYCOSYLATED/   # Glycosylated protein structures
│   │   ├── protein_asn_orientation.pdb
│   │   ├── protein_glycosylated.pdb
│   │   └── protein_glycosylated_renumbered.pdb
│   ├── EXTRACTED_CARBOHYDRATES/    # Extracted glycan structures
│   └── TO_TOP/                     # Temporary files
│
├── STEP2/                          # Parametrization
│   ├── JSON/                       # Glycan JSON representations
│   ├── PDB_GLYCOPROTEIN/           # Final glycoprotein structures
│   ├── VALENCE_GLYCAN_VARIANTS/    # Glycan variants
│   └── charmm36.ff/                # CHARMM36 force field
│       ├── carb.rtp                # Updated RTP file
│       └── carb.hdb                # Updated HDB file
│
├── STEP3/                          # Carbohydrate orientation
│   ├── JSON_FILES/                 # Intermediate JSON files
│   ├── PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/
│   │   ├── protein_optimized.pdb   # Final optimized structure
│   │   ├── report.txt              # Optimization report
│   │   └── PDB_CARB_ONLY/          # Individual optimized glycans
│
└── report.txt                      # Main report file
```

## Understanding the Report

The optimization report (`report.txt`) contains:

- **Initial Energy**: Starting conformation energy
- **Final Energy**: Optimized conformation energy
- **Acceptance Rate**: MCMC acceptance ratio
- **Clash Statistics**: Number of clashes before/after
- **Dihedral Angles**: Final glycan conformations

Example report snippet:
```
========================================
GLYCAN OPTIMIZATION REPORT
========================================
Glycan: G1 (Chain A, Residue 45)
  Initial clashes: 12
  Final clashes: 0
  Initial energy: 45.23 kcal/mol
  Final energy: 12.45 kcal/mol
  Acceptance rate: 0.78
  Cycles completed: 5
  
Optimized dihedrals:
  phi: -67.3°
  psi: 152.8°
  omega: 178.2°
========================================
```

## Examples

### Example 1: Basic glycoprotein modeling

```bash
auto_glyco \
    --pdb examples/antibody.pdb \
    --tsv examples/glycosylation_sites.tsv \
    --output_dir ./antibody_glycosylated
```

### Example 2: High-precision optimization

```bash
auto_glyco \
    --pdb complex.pdb \
    --tsv sites.tsv \
    --theta_step 2 \
    --n_steps 50 \
    --max_cycles 20 \
    --radius 150 \
    --use_coulomb yes \
    --n_workers 16 \
    --save_individual_glycans \
    --save_before_after \
    --verbose
```

### Example 3: Quick screening

```bash
auto_glyco \
    --pdb screening.pdb \
    --tsv candidates.tsv \
    --output_dir ./screening \
    --theta_step 30 \
    --n_steps 5 \
    --max_cycles 3 \
    --n_workers 4
```

## Citation

If you use this pipeline in your research, please cite:

Anacleto Silva de Souza, Cristiane Rodrigues Guzzo and Siewert-Jan Marrink, Auto_glyco: Automated Glycosylation Pipeline, 2026, https://github.com/anacletosouza/automated_glycosylation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. This tool was supported by São Paulo Research Foundation (FAPESP) 2023/18211-0 and 2025/05583-1.

## Authors

- Anacleto Silva de Souza, Cristiane Rodrigues Guzzo and Siewert-Jan Marrink,
- email: anacletosilvadesouza@usp.br

## Contact

For questions, issues, or contributions, please:
- Open an issue on [GitHub](https://github.com/anacletosouza/automated_glycosylation/issues)
- Email the author at anacletosilvadesouza@usp.br

