# Automated Glycosylation Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

A comprehensive computational pipeline for automated glycosylation of proteins, including structure preparation, force field parametrization, and carbohydrate orientation optimization using Markov Chain Monte Carlo (MCMC) methods.

## Table of Contents

- [Overview](#overview)
- [Theoretical Background](#theoretical-background)
  - [Glycosylation in Structural Biology](#glycosylation-in-structural-biology)
  - [Markov Chain Monte Carlo (MCMC) for Glycan Orientation](#markov-chain-monte-carlo-mcmc-for-glycan-orientation)
  - [Metropolis-Hastings Algorithm](#metropolis-hastings-algorithm)
  - [Energy Functions and Force Fields](#energy-functions-and-force-fields)
- [Installation](#installation)
- [Pipeline Structure](#pipeline-structure)
- [Command Line Interface](#command-line-interface)
  - [Step 1: Glycosylation Preparation](#step-1-glycosylation-preparation)
  - [Step 2: Parametrization](#step-2-parametrization)
  - [Step 3: Carbohydrate Orientation](#step-3-carbohydrate-orientation)
  - [Complete Pipeline](#complete-pipeline)
- [Examples](#examples)
- [Output Structure](#output-structure)
- [Mathematical Details](#mathematical-details)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

## Overview

This pipeline automates the process of adding glycans to protein structures, generating force field parameters, and optimizing carbohydrate orientations using statistical mechanics principles. It is designed for structural biologists, computational chemists, and glycobiologists working with glycoprotein modeling.

## Theoretical Background

### Glycosylation in Structural Biology

Glycosylation is a post-translational modification where carbohydrate moieties (glycans) are attached to specific amino acid residues (typically asparagine in N-linked glycosylation). The orientation of these glycans significantly affects:

- Protein stability and folding
- Molecular recognition events
- Immune response modulation
- Drug efficacy and bioavailability

### Markov Chain Monte Carlo (MCMC) for Glycan Orientation

The orientation of glycans attached to proteins is governed by complex energy landscapes with multiple local minima. Traditional deterministic optimization methods often get trapped in local minima. MCMC provides a robust framework for exploring these energy landscapes.

#### Theoretical Foundation

Consider a glycan with $n$ rotatable dihedral angles $\boldsymbol{\theta} = (\theta_1, \theta_2, ..., \theta_n)$. The probability distribution of conformations follows the Boltzmann distribution:

$$
P(\boldsymbol{\theta}) = \frac{1}{Z} \exp\left(-\frac{E(\boldsymbol{\theta})}{k_B T}\right)
$$

where:
- $E(\boldsymbol{\theta})$ is the potential energy of the conformation
- $k_B$ is Boltzmann's constant
- $T$ is the absolute temperature
- $Z = \int \exp(-E(\boldsymbol{\theta})/k_B T) d\boldsymbol{\theta}$ is the partition function

### Metropolis-Hastings Algorithm

The Metropolis-Hastings algorithm generates a Markov chain that samples from $P(\boldsymbol{\theta})$ without requiring knowledge of $Z$.

**Algorithm Steps:**

1. **Initialize** starting conformation $\boldsymbol{\theta}_0$

2. **For each iteration** $t = 1, 2, ..., N$:
   
   a. **Propose** a new conformation $\boldsymbol{\theta}'$ from proposal distribution $q(\boldsymbol{\theta}' | \boldsymbol{\theta}_{t-1})$
   
   b. **Calculate acceptance probability**:
   
   $$
   \alpha = \min\left(1, \frac{P(\boldsymbol{\theta}')}{P(\boldsymbol{\theta}_{t-1})} \cdot \frac{q(\boldsymbol{\theta}_{t-1} | \boldsymbol{\theta}')}{q(\boldsymbol{\theta}' | \boldsymbol{\theta}_{t-1})}\right)
   $$
   
   For symmetric proposal distributions (e.g., Gaussian random walk), $q(\boldsymbol{\theta}_{t-1} | \boldsymbol{\theta}') = q(\boldsymbol{\theta}' | \boldsymbol{\theta}_{t-1})$, simplifying to:
   
   $$
   \alpha = \min\left(1, \frac{P(\boldsymbol{\theta}')}{P(\boldsymbol{\theta}_{t-1})}\right) = \min\left(1, \exp\left(-\frac{\Delta E}{k_B T}\right)\right)
   $$
   
   c. **Accept or reject**:
   - Generate random number $u \sim \mathcal{U}(0,1)$
   - If $u \leq \alpha$, accept: $\boldsymbol{\theta}_t = \boldsymbol{\theta}'$
   - Else, reject: $\boldsymbol{\theta}_t = \boldsymbol{\theta}_{t-1}$

3. **After convergence**, the samples $\{\boldsymbol{\theta}_t\}$ approximate the Boltzmann distribution.

### Energy Functions and Force Fields

The total potential energy is calculated using the CHARMM36 force field:

$$
E_{\text{total}} = E_{\text{bond}} + E_{\text{angle}} + E_{\text{dihedral}} + E_{\text{improper}} + E_{\text{nonbonded}}
$$

#### Bond Stretching (Harmonic oscillator approximation)

$$
E_{\text{bond}} = \sum_{\text{bonds}} k_b (r - r_0)^2
$$

where $k_b$ is the bond force constant, $r$ is the current bond length, and $r_0$ is the equilibrium bond length.

#### Angle Bending

$$
E_{\text{angle}} = \sum_{\text{angles}} k_\theta (\theta - \theta_0)^2
$$

where $k_\theta$ is the angle force constant, $\theta$ is the current angle, and $\theta_0$ is the equilibrium angle.

#### Dihedral Torsions

$$
E_{\text{dihedral}} = \sum_{\text{dihedrals}} k_\phi [1 + \cos(n\phi - \delta)]
$$

where $k_\phi$ is the dihedral force constant, $n$ is the multiplicity, $\phi$ is the dihedral angle, and $\delta$ is the phase shift.

#### Non-bonded Interactions

**Lennard-Jones potential (van der Waals):**

$$
E_{\text{vdW}} = \sum_{i<j} 4\varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right]
$$

where $\varepsilon_{ij}$ is the well depth, $\sigma_{ij}$ is the distance at zero potential, and $r_{ij}$ is the distance between atoms $i$ and $j$.

**Coulomb potential (electrostatics):**

$$
E_{\text{elec}} = \sum_{i<j} \frac{q_i q_j}{4\pi\epsilon_0 \epsilon_r r_{ij}}
$$

where $q_i$ and $q_j$ are partial atomic charges, $\epsilon_0$ is the vacuum permittivity, and $\epsilon_r$ is the relative permittivity.

#### Total Non-bonded Energy

$$
E_{\text{nonbonded}} = E_{\text{vdW}} + E_{\text{elec}}
$$

## Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager
- GROMACS (optional, for molecular dynamics simulations)

### Install from GitHub

```bash
# Clone the repository
git clone https://github.com/anacletosouza/automated_glycosylation.git
cd automated_glycosylation

# Install in development mode
pip install -e .

# Or install normally
pip install .
```

### Verify Installation

```bash
# Check available commands
glyco-prep --help
glyco-param --help
glyco-orient --help
glyco-all --help
```

## Pipeline Structure

```
INPUT PDB File
    │
    ├──> Step 1: Glycosylation Preparation
    │    ├── Asparagine orientation optimization
    │    ├── Glycan attachment
    │    ├── Chain/residue renumbering
    │    └── Glycan coordinate extraction
    │
    ├──> Step 2: Parametrization
    │    ├── CHARMM36 force field setup
    │    ├── JSON generation for each glycan
    │    ├── RTP file generation
    │    ├── Topology unification
    │    └── Glycoprotein construction
    │
    └──> Step 3: Carbohydrate Orientation
         ├── PDB to JSON conversion
         ├── Force field parameter integration
         └── MCMC optimization of glycans
```

## Command Line Interface

### Step 1: Glycosylation Preparation

Adds glycans to protein structure and prepares initial geometries.

```bash
glyco-prep -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input`: Input PDB file (protein structure without glycans)
- `-o, --output-dir`: Output directory for results

**Optional arguments:**
- `--asn-tsv`: TSV file with glycosylation sites (if not provided, auto-detects)
- `--rotate-atoms`: Atoms to rotate for asparagine orientation (default: "OD1,CG,ND2,HD22,HD21,HB2,HB3")
- `--fixed-atom`: Fixed atom for rotation (default: "CB")
- `--center-atom`: Center atom for rotation (default: "CA")
- `--radius`: Radius for neighbor detection in Angstroms (default: 30.0)
- `--rotation-step`: Rotation step in degrees (default: 1)
- `--protein-residue-start`: Starting residue number for protein (default: 1)
- `--keep-temp`: Keep temporary files

### Step 2: Parametrization

Generates force field parameters and topology files.

```bash
glyco-param -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input-pdb`: Input PDB file from Step 1 (glycosylated_protein_renumbered.pdb)
- `-o, --output-dir`: Output directory for topology files

**Optional arguments:**
- `--download-charmm`: Download CHARMM36 force field automatically
- `--charmm-url`: Custom URL for CHARMM download
- `--force-download`: Force download even if backup exists
- `--n-cpus`: Number of CPUs for parallel processing (default: 1)
- `--keep-intermediate`: Keep intermediate files

### Step 3: Carbohydrate Orientation

Optimizes glycan orientations using MCMC.

```bash
glyco-orient -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input-pdb`: Input PDB file from Step 2
- `-o, --output-dir`: Output directory for optimized structures

**MCMC Parameters:**
- `--theta-step`: Dihedral angle step size in degrees (default: 10)
- `--n-steps`: Number of MCMC steps per cycle (default: 10)
- `--max-cycles`: Maximum number of optimization cycles (default: 5)
- `--radius`: Interaction radius in Angstroms (default: 300.0)
- `--use-coulomb`: Include Coulomb electrostatics (default: false)

**Parallel Processing:**
- `--n-workers`: Number of CPU workers for parallel MCMC (default: 1)

**Output Options:**
- `--save-individual-glycans`: Save individual glycan PDB files
- `--save-before-after`: Save before/after comparison files
- `--verbose`: Verbose output
- `--report-file`: Custom report file path

**Force Field Options:**
- `--charmm-dir`: Custom CHARMM36 force field directory

### Complete Pipeline

Runs all three steps sequentially.

```bash
glyco-all -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input`: Input PDB file
- `-o, --output-dir`: Output directory for all results

**General Options:**
- `--asn-tsv`: TSV file with glycosylation sites
- `--download-charmm`: Download CHARMM36 force field
- `--n-cpus`: CPUs for parametrization (default: 1)
- `--n-workers`: Workers for MCMC (default: 1)
- `--keep-temp`: Keep temporary files
- `--verbose`: Verbose output

**Asparagine Orientation Options:**
- `--rotate-atoms`: Atoms to rotate (default: "OD1,CG,ND2,HD22,HD21,HB2,HB3")
- `--fixed-atom`: Fixed atom (default: "CB")
- `--center-atom`: Center atom (default: "CA")
- `--radius`: Neighbor radius in Angstroms (default: 30.0)
- `--rotation-step`: Rotation step in degrees (default: 1)

**MCMC Options:**
- `--theta-step`: Dihedral step size (default: 10)
- `--n-steps`: MCMC steps per cycle (default: 10)
- `--max-cycles`: Maximum cycles (default: 5)
- `--mcmc-radius`: Interaction radius in Angstroms (default: 300.0)
- `--use-coulomb`: Include electrostatics (default: false)

## Examples

### Example 1: Basic Glycosylation of a Protein

```bash
# Add glycans to protein structure
glyco-prep -i protein.pdb -o glycosylation_results

# Check outputs
ls glycosylation_results/
# PDB_PROTEIN_GLYCOSYLATED/  TO_TOP/  TSV/
```

### Example 2: Full Pipeline with Custom Parameters

```bash
# Run complete pipeline with custom settings
glyco-all -i antibody.pdb -o antibody_glycosylation \
    --asn-tsv glycosylation_sites.tsv \
    --download-charmm \
    --n-cpus 8 \
    --n-workers 8 \
    --theta-step 5 \
    --n-steps 20 \
    --max-cycles 10 \
    --radius 25.0 \
    --rotation-step 2 \
    --verbose
```

### Example 3: Individual Step Execution

```bash
# Step 1: Prepare glycosylated protein
glyco-prep -i protein.pdb -o step1_output \
    --rotate-atoms "OD1,CG,ND2,HD22,HD21" \
    --radius 35.0 \
    --rotation-step 2

# Step 2: Generate parameters
glyco-param -i step1_output/PDB_PROTEIN_GLYCOSYLATED/glycosylated_protein_renumbered.pdb \
    -o step2_output \
    --download-charmm \
    --n-cpus 8

# Step 3: Optimize orientations
glyco-orient -i step2_output/VALENCE_GLYCAN_VARIANTS/glycosylated_protein_final_valence_corrected_variants.pdb \
    -o step3_output \
    --theta-step 10 \
    --n-steps 50 \
    --max-cycles 20 \
    --save-individual-glycans \
    --save-before-after \
    --n-workers 8
```

### Example 4: High-Throughput Processing

```bash
# Process multiple proteins
for protein in *.pdb; do
    name="${protein%.pdb}"
    glyco-all -i "$protein" -o "${name}_glycosylation" \
        --download-charmm \
        --n-cpus 4 \
        --n-workers 4 \
        --keep-temp
done
```

### Example 5: MCMC Convergence Analysis

```bash
# Run with detailed output for convergence analysis
glyco-orient -i glycoprotein.pdb -o mcmc_analysis \
    --theta-step 5 \
    --n-steps 1000 \
    --max-cycles 1 \
    --save-individual-glycans \
    --verbose \
    --report-file convergence.txt

# The report file contains energy vs. step data for analysis
```

## Output Structure

### Step 1 Output

```
OUTPUT_DIR/
├── PDB_PROTEIN_GLYCOSYLATED/
│   ├── asn_orientation.pdb
│   ├── glycosylated_protein.pdb
│   ├── glycosylated_protein_renumbered.pdb
│   └── glycosylated_protein_renumbered_without_H.pdb
├── TO_TOP/
│   └── [individual glycan PDB files]
├── TSV/
│   └── [glycosylation site data]
└── EXTRACTED_CARBOHYDRATES/
```

### Step 2 Output

```
OUTPUT_DIR/
├── charmm36.ff/              # CHARMM36 force field
│   ├── carb.rtp             # Modified RTP file
│   ├── carb.hdb             # Modified HDB file
│   └── *.backup             # Backup files
├── JSON/                     # Glycan JSON files
│   ├── [glycan_name]/
│   │   ├── [glycan].pdb
│   │   ├── [glycan].json
│   │   ├── *.pkl
│   │   ├── *.rtp
│   │   └── *.hdb
│   └── carb_redundance_removed.rtp
├── PDB_GLYCOPROTEIN/
│   ├── glycosylated_protein_corrected.pdb
│   ├── glycosylated_protein_final_connected.pdb
│   └── glycosylated_protein_final_valence_corrected.pdb
└── VALENCE_GLYCAN_VARIANTS/
    ├── *.rtp                # Variant RTP files
    ├── *.hdb                # Variant HDB files
    └── *.pdb                # Final glycoprotein structures
```

### Step 3 Output

```
OUTPUT_DIR/
├── JSON_FILES/
│   ├── pdb_to_json.json
│   └── glycan_data_charmm36.json
├── PDB_CARBOHYDRATE_ORIENTATION_OPTIMIZED/
│   ├── glycan_optimized.json
│   ├── glycosylated_protein_final_optimized.pdb
│   ├── report.txt           # MCMC convergence report
│   └── PDB_CARB_ONLY/       # Individual optimized glycans
│       ├── [glycan]_before.pdb
│       ├── [glycan]_after.pdb
│       └── [glycan]_trajectory/
│           ├── step_*.pdb
│           └── energy_*.dat
```

## Mathematical Details

### MCMC Acceptance Probability Derivation

The Metropolis-Hastings acceptance ratio ensures detailed balance:

$$
\pi(\boldsymbol{\theta}) P(\boldsymbol{\theta} \rightarrow \boldsymbol{\theta}') = \pi(\boldsymbol{\theta}') P(\boldsymbol{\theta}' \rightarrow \boldsymbol{\theta})
$$

where $\pi(\boldsymbol{\theta})$ is the target distribution and $P(\boldsymbol{\theta} \rightarrow \boldsymbol{\theta}')$ is the transition probability:

$$
P(\boldsymbol{\theta} \rightarrow \boldsymbol{\theta}') = q(\boldsymbol{\theta}' | \boldsymbol{\theta}) \alpha(\boldsymbol{\theta}, \boldsymbol{\theta}')
$$

Substituting:

$$
\pi(\boldsymbol{\theta}) q(\boldsymbol{\theta}' | \boldsymbol{\theta}) \alpha(\boldsymbol{\theta}, \boldsymbol{\theta}') = \pi(\boldsymbol{\theta}') q(\boldsymbol{\theta} | \boldsymbol{\theta}') \alpha(\boldsymbol{\theta}', \boldsymbol{\theta})
$$

Solving for $\alpha$:

$$
\frac{\alpha(\boldsymbol{\theta}, \boldsymbol{\theta}')}{\alpha(\boldsymbol{\theta}', \boldsymbol{\theta})} = \frac{\pi(\boldsymbol{\theta}') q(\boldsymbol{\theta} | \boldsymbol{\theta}')}{\pi(\boldsymbol{\theta}) q(\boldsymbol{\theta}' | \boldsymbol{\theta})}
$$

The Metropolis choice is:

$$
\alpha(\boldsymbol{\theta}, \boldsymbol{\theta}') = \min\left(1, \frac{\pi(\boldsymbol{\theta}') q(\boldsymbol{\theta} | \boldsymbol{\theta}')}{\pi(\boldsymbol{\theta}) q(\boldsymbol{\theta}' | \boldsymbol{\theta})}\right)
$$

### Energy Minimization in MCMC

The potential energy for a glycan conformation is computed as:

$$
E(\boldsymbol{\theta}) = E_{\text{intra}}(\boldsymbol{\theta}) + E_{\text{inter}}(\boldsymbol{\theta})
$$

where:
- $E_{\text{intra}}$ is the internal energy of the glycan (bond, angle, dihedral terms)
- $E_{\text{inter}}$ is the interaction energy with the protein (van der Waals + electrostatics)

### Convergence Criteria

The MCMC simulation is considered converged when:

1. **Gelman-Rubin statistic** $\hat{R} < 1.1$:

$$
\hat{R} = \sqrt{\frac{\text{Var}(\boldsymbol{\theta})}{W}}
$$

where $\text{Var}(\boldsymbol{\theta})$ is the between-chain variance and $W$ is the within-chain variance.

2. **Autocorrelation time** $\tau < 10$ steps:

$$
\tau = 1 + 2 \sum_{k=1}^{\infty} \rho_k
$$

where $\rho_k$ is the autocorrelation at lag $k$.

3. **Energy stabilization**: The running average of $\Delta E$ changes by < 1% over 100 steps.

### Temperature Schedule (Simulated Annealing)

The pipeline optionally implements a cooling schedule:

$$
T_k = T_0 \cdot \beta^k
$$

where:
- $T_0$ is the initial temperature (default: 300 K)
- $\beta$ is the cooling factor (default: 0.99)
- $k$ is the step number

This allows the system to escape local minima and find global energy minima.

## Troubleshooting

### Common Issues and Solutions

**Issue**: CHARMM force field download fails

**Solution**: Manually download using a mirror URL
```bash
glyco-param -i input.pdb -o output --charmm-url "your_mirror_url"
```

**Issue**: MCMC not converging

**Solution**: Increase steps and cycles
```bash
glyco-orient -i input.pdb -o output --n-steps 100 --max-cycles 20
```

**Issue**: Memory error during parametrization

**Solution**: Reduce parallel workers
```bash
glyco-param -i input.pdb -o output --n-cpus 2
```

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{automated_glycosylation_2024,
  author = {Silva de Souza, Anacleto},
  title = {Automated Glycosylation Pipeline for Glycoproteins},
  year = {2024},
  url = {https://github.com/anacletosouza/automated_glycosylation},
  doi = {10.5281/zenodo.xxxxxxx}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- CHARMM36 force field developers
- GROMACS community
- São Paulo Research Foundation (FAPESP) for funding

---

**Contact**: anacletosilvadesouza@usp.br

**GitHub**: [https://github.com/anacletosouza/automated_glycosylation](https://github.com/anacletosouza/automated_glycosylation)

