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
- [MCMC Model Development](#mcmc-model-development)
  - [State Space and Sampling Distribution](#state-space-and-sampling-distribution)
  - [Proposal Distribution and Adaptive Tuning](#proposal-distribution-and-adaptive-tuning)
  - [Convergence Diagnostics](#convergence-diagnostics)
  - [Parallel Tempering for Enhanced Sampling](#parallel-tempering-for-enhanced-sampling)
- [Installation](#installation)
- [Pipeline Structure](#pipeline-structure)
- [Command Line Interface](#command-line-interface)
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

Glycosylation is a post-translational modification where carbohydrate moieties (glycans) are attached to specific amino acid residues (typically asparagine for N-linked glycosylation, serine/threonine for O-linked glycosylation). The orientation of these glycans significantly affects:

- Protein stability and folding
- Molecular recognition events
- Immune response modulation
- Drug efficacy and bioavailability

### Markov Chain Monte Carlo (MCMC) for Glycan Orientation

The orientation of glycans attached to proteins is governed by complex energy landscapes with multiple local minima. Traditional deterministic optimization methods often get trapped in local minima. MCMC provides a robust framework for exploring these energy landscapes.

#### Theoretical Foundation

Consider a glycan with $n$ rotatable dihedral angles $\boldsymbol{\theta} = (\theta_1, \theta_2, ..., \theta_n)$. The probability distribution of conformations follows the Boltzmann distribution:

$P(\boldsymbol{\theta}) = \frac{1}{Z} \exp\left(-\frac{E(\boldsymbol{\theta})}{k_B T}\right)$

where:
- $E(\boldsymbol{\theta})$ is the potential energy of the conformation
- $k_B$ is Boltzmann's constant ($0.008314462618$ kJ/mol·K)
- $T$ is the absolute temperature (default: 300 K)
- $Z = \int \exp(-E(\boldsymbol{\theta})/k_B T) d\boldsymbol{\theta}$ is the partition function

### Metropolis-Hastings Algorithm

The Metropolis-Hastings algorithm generates a Markov chain that samples from $P(\boldsymbol{\theta})$ without requiring knowledge of $Z$.

**Algorithm Steps:**

1. **Initialize** starting conformation $\boldsymbol{\theta}_0$ (from grid search)

2. **For each iteration** $t = 1, 2, ..., N$:
   
   a. **Propose** a new conformation $\boldsymbol{\theta}'$ from proposal distribution $q(\boldsymbol{\theta}' | \boldsymbol{\theta}_{t-1})$
   
   b. **Calculate acceptance probability**:
   
   $\alpha = \min\left(1, \frac{P(\boldsymbol{\theta}')}{P(\boldsymbol{\theta}_{t-1})} \cdot \frac{q(\boldsymbol{\theta}_{t-1} | \boldsymbol{\theta}')}{q(\boldsymbol{\theta}' | \boldsymbol{\theta}_{t-1})}\right)$
   
   For symmetric proposal distributions, this simplifies to:
   
   $\alpha = \min\left(1, \frac{P(\boldsymbol{\theta}')}{P(\boldsymbol{\theta}_{t-1})}\right) = \min\left(1, \exp\left(-\frac{\Delta E}{k_B T}\right)\right)$
   
   c. **Accept or reject**:
   - Generate random number $u \sim \mathcal{U}(0,1)$
   - If $u \leq \alpha$, accept: $\boldsymbol{\theta}_t = \boldsymbol{\theta}'$
   - Else, reject: $\boldsymbol{\theta}_t = \boldsymbol{\theta}_{t-1}$

3. **After convergence**, the samples $\{\boldsymbol{\theta}_t\}$ approximate the Boltzmann distribution.

### Energy Functions and Force Fields

The total potential energy is calculated using the CHARMM36 force field:

$E_{\text{total}} = E_{\text{vdW}} + E_{\text{coulomb}}$

#### Van der Waals Energy (Lennard-Jones)

$E_{\text{vdW}} = \sum_{i<j} 4\varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right]$

where $\varepsilon_{ij}$ is the well depth, $\sigma_{ij}$ is the distance at zero potential (Lorentz-Berthelot mixing rules: $\sigma_{ij} = (\sigma_i + \sigma_j)/2$, $\varepsilon_{ij} = \sqrt{\varepsilon_i \varepsilon_j}$), and $r_{ij}$ is the distance between atoms $i$ and $j$.

#### Coulomb Energy (Electrostatics)

$E_{\text{coulomb}} = \sum_{i<j} \frac{q_i q_j}{4\pi\epsilon_0 \epsilon_r r_{ij}}$

where $q_i$ and $q_j$ are partial atomic charges, $\epsilon_0$ is the vacuum permittivity, and $\epsilon_r$ is the relative permittivity.

## MCMC Model Development

### State Space and Sampling Distribution

The conformational space of a glycan attached to a protein is defined by:

$\Omega = \{\boldsymbol{\theta} \in [0, 2\pi)^d\}$

where $d$ is the number of rotatable dihedral angles. For N-linked glycans, the primary rotation degree of freedom is the glycosidic bond between the protein (ASN ND2) and the first glycan residue (C1). Our pipeline simplifies this high-dimensional space to a single rotational degree of freedom $\theta$ around the axis defined by:

$\mathbf{a} = \frac{\mathbf{r}_{C1} - \mathbf{r}_{\text{ND2}}}{|\mathbf{r}_{C1} - \mathbf{r}_{\text{ND2}}|}$

The target sampling distribution is the Boltzmann distribution:

$\pi(\theta) = \frac{1}{Z} \exp\left(-\frac{E(\theta)}{k_B T}\right)$

where $Z = \int_0^{2\pi} \exp(-E(\theta)/k_B T) d\theta$ is the normalization constant.

### Proposal Distribution and Adaptive Tuning

#### Symmetric Random Walk Proposal

We employ a symmetric Gaussian random walk proposal distribution:

$q(\theta' | \theta_t) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(\theta' - \theta_t)^2}{2\sigma^2}\right)$

where $\sigma$ is the step size (default: 10 degrees converted to radians). The symmetry property $q(\theta' | \theta_t) = q(\theta_t | \theta')$ simplifies the Metropolis-Hastings acceptance ratio to:

$\alpha = \min\left(1, \exp\left(-\frac{E(\theta') - E(\theta_t)}{k_B T}\right)\right)$

#### Adaptive Metropolis

The pipeline optionally implements adaptive MCMC where the proposal distribution covariance is tuned during burn-in:

$\Sigma_t = \text{Cov}(\theta_0, ..., \theta_t) + \epsilon I_d$

where $\epsilon = 10^{-6}$ ensures positive definiteness. The step size is adjusted to maintain an optimal acceptance rate of $0.234$ (theoretical optimum for high-dimensional problems):

$\sigma_{t+1} = \sigma_t \cdot \exp\left(\gamma_t (\alpha_t - \alpha^*)\right)$

where:
- $\alpha_t$ is the empirical acceptance rate over the last $W$ steps
- $\alpha^* = 0.234$ is the target acceptance rate
- $\gamma_t = t^{-0.6}$ is the step size adaptation gain

### Convergence Diagnostics

#### Gelman-Rubin Potential Scale Reduction Factor

For $M$ parallel chains, compute:

$\hat{R} = \sqrt{\frac{\text{Var}(\theta)}{W}}$

where:
- $W = \frac{1}{M}\sum_{j=1}^M s_j^2$ is the within-chain variance
- $s_j^2$ is the variance of chain $j$
- $\text{Var}(\theta) = \frac{N-1}{N}W + \frac{1}{N}B$ is the estimated posterior variance
- $B$ is the between-chain variance

Convergence is achieved when $\hat{R} < 1.1$ for all parameters.

#### Effective Sample Size (ESS)

The effective sample size accounts for autocorrelation in the Markov chain:

$ESS = \frac{N}{1 + 2\sum_{k=1}^{\infty} \rho_k}$

where $\rho_k$ is the autocorrelation at lag $k$. We consider $ESS \geq 100$ as sufficient for reliable inference.


### Geweke Diagnostic

The Geweke test compares the mean of the first 10% of samples to the last 50%:

$\displaystyle z = \frac{\bar{\theta}_A - \bar{\theta}_B}{\sqrt{\mathrm{Var}(\theta_A) + \mathrm{Var}(\theta_B)}}$

Under convergence, $z \sim \mathcal{N}(0,1)$. Values $|z| > 2$ indicate non-convergence.

---

### Energy Convergence Criterion

For glycoprotein systems, we monitor the running average of potential energy:

$\displaystyle \bar{E}*t = \frac{1}{t}\sum*{i=1}^{t} E(\theta_i)$

The chain is considered converged when:

$\displaystyle \frac{|\bar{E}*t - \bar{E}*{t-100}|}{\bar{E}_t} < 0.01$

indicating energy stabilization within 1% over 100 steps.

### Parallel Tempering for Enhanced Sampling

To overcome energy barriers and sample multimodal distributions, the pipeline implements parallel tempering (replica exchange MCMC) with $M$ replicas at different temperatures:

$T_i = T_0 \cdot \gamma^{i-1}, \quad i = 1, 2, ..., M$

where:
- $T_0 = 300$ K is the base temperature
- $\gamma = 1.2$ ensures overlapping temperature distributions
- $M = \lceil \log(T_{\max}/T_0)/\log(\gamma) \rceil$ replicas

#### Replica Exchange Probability

At regular intervals, adjacent replicas $i$ and $j$ (where $j = i+1$ and $T_j > T_i$) attempt exchange with probability:

$P_{\text{swap}} = \min\left(1, \exp\left[\left(\frac{1}{k_B T_i} - \frac{1}{k_B T_j}\right)(E_j - E_i)\right]\right)$

This maintains detailed balance:

$\pi_i(\theta_i)\pi_j(\theta_j) P_{\text{swap}}(\theta_i, \theta_j) = \pi_i(\theta_j)\pi_j(\theta_i) P_{\text{swap}}(\theta_j, \theta_i)$

#### Round-Robin Exchange Schedule

The pipeline implements a round-robin exchange schedule:
1. Exchange between replicas (1,2), (3,4), (5,6), ...
2. Exchange between replicas (2,3), (4,5), (6,7), ...

This ensures that all replicas have opportunity to exchange, improving mixing and allowing low-temperature replicas to escape local minima by swapping with high-temperature replicas.

### Statistical Validation Metrics

#### Potential Energy Distribution Analysis

After convergence, the distribution of potential energies should approximate:

$P(E) \propto \Omega(E) e^{-E/k_B T}$

where $\Omega(E)$ is the density of states. We validate using:

- **Mean energy**: $\langle E \rangle = \frac{1}{N}\sum_{t=1}^N E(\theta_t)$
- **Energy variance**: $\text{Var}(E) = \langle E^2 \rangle - \langle E \rangle^2 = C_V k_B T^2$
- **Heat capacity**: $C_V = \frac{\partial \langle E \rangle}{\partial T} = \frac{\text{Var}(E)}{k_B T^2}$

#### Autocorrelation Analysis

The integrated autocorrelation time is computed as:

$\tau_{\text{int}} = \frac{1}{2} + \sum_{k=1}^{\infty} \rho(k)$

where $\rho(k) = \frac{\text{Cov}(\theta_t, \theta_{t+k})}{\text{Var}(\theta_t)}$ is the autocorrelation at lag $k$. The effective sample size is then:

$ESS = \frac{N}{\tau_{\text{int}}}$

### Algorithm Implementation Details

```python
def mcmc_optimization(theta_initial, energy_function, n_steps=10000, sigma=10.0, temperature=300):
    """
    MCMC optimization using symmetric random walk proposal.
    
    Parameters:
    -----------
    theta_initial : float
        Initial theta angle in degrees
    energy_function : callable
        Function that returns energy for given theta
    n_steps : int
        Number of MCMC steps
    sigma : float
        Proposal step size in degrees
    temperature : float
        Temperature in Kelvin
    
    Returns:
    --------
    theta_best : float
        Best theta angle found
    acceptance_rate : float
        Fraction of accepted proposals
    energy_history : list
        History of energy values
    """
    # Initialize
    theta_current = np.deg2rad(theta_initial)
    energy_current = energy_function(theta_current)
    
    best_theta = theta_current
    best_energy = energy_current
    
    acceptance_count = 0
    energy_history = [energy_current]
    theta_history = [theta_current]
    
    # Precompute Boltzmann factor scaling
    beta = 1.0 / (KB * temperature)
    
    for step in range(n_steps):
        # Propose new theta
        delta_theta = np.random.normal(0, sigma_rad)
        theta_proposed = theta_current + delta_theta
        
        # Apply periodic boundary conditions
        theta_proposed = np.mod(theta_proposed, 2 * np.pi)
        
        # Calculate energy of proposed state
        energy_proposed = energy_function(theta_proposed)
        
        # Metropolis acceptance criterion
        delta_energy = energy_proposed - energy_current
        
        if delta_energy < 0:
            accept = True
        else:
            acceptance_probability = np.exp(-beta * delta_energy)
            accept = np.random.random() < acceptance_probability
        
        # Update state
        if accept:
            theta_current = theta_proposed
            energy_current = energy_proposed
            acceptance_count += 1
            
            # Track best state
            if energy_current < best_energy:
                best_theta = theta_current
                best_energy = energy_current
        
        # Record history
        energy_history.append(energy_current)
        theta_history.append(theta_current)
    
    acceptance_rate = acceptance_count / n_steps
    return np.rad2deg(best_theta), best_energy, energy_history, acceptance_rate
```

## Installation

### Prerequisites

- **Python 3.7 or higher** - The pipeline requires Python 3.7+ for full compatibility with type hints and dataclasses
- **pip package manager** - For installing Python dependencies
- **GROMACS (optional)** - Only required if you plan to run molecular dynamics simulations with the generated topology files

### Python Dependencies

The pipeline requires the following Python packages:

| Package | Minimum Version | Purpose |
|---------|----------------|---------|
| `numpy` | 1.19.0+ | Numerical operations, rotation matrices, distance calculations |
| `scipy` | 1.5.0+ | Scientific computing (optional, used for advanced statistics) |
| `pandas` | 1.1.0+ | Reading TSV files containing glycosylation site tables |
| `matplotlib` | 3.3.0+ | SNFG (Symbol Nomenclature for Glycans) visualization |
| `mdtraj` | 1.9.0+ | Molecular trajectory manipulation (required by glycosylator) |
| `glycosylator` | 0.1.0+ | Core library for protein glycosylation operations |
| `tqdm` | 4.50.0+ | Progress bars for MCMC simulations |

### System Requirements

- **RAM**: 
  - Minimum: 2 GB for single chain proteins
  - Recommended: 8 GB for 4 parallel workers
  - Large systems (>1000 residues): 16+ GB

- **CPU**: 
  - Multi-core processor recommended for parallel MCMC processing
  - The pipeline automatically detects and uses available CPU cores

- **Disk Space**: 
  - Minimum: 500 MB for CHARMM36 force field and dependencies
  - Additional space required for output files (varies with system size)

### Install Dependencies

```bash
# Install all required packages via pip
pip install numpy>=1.19.0 scipy>=1.5.0 pandas>=1.1.0 matplotlib>=3.3.0 mdtraj>=1.9.0 tqdm>=4.50.0 glycosylator

# Or use the requirements file
pip install -r requirements.txt

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
    │    ├── Asparagine side chain optimization (HD22 orientation)
    │    ├── Glycan attachment to ASN (N-linked) or SER/THR (O-linked)
    │    ├── Chain and residue renumbering
    │    └── Glycan coordinate extraction
    │
    ├──> Step 2: Parametrization
    │    ├── CHARMM36 force field setup
    │    ├── JSON generation for each glycan
    │    ├── RTP file generation for GROMACS
    │    ├── Topology unification
    │    └── Glycoprotein construction
    │
    └──> Step 3: Carbohydrate Orientation
         ├── PDB to JSON conversion with CHARMM36 parameters
         ├── Single-axis grid search (0-360°)
         ├── MCMC refinement with Metropolis-Hastings
         ├── Optional parallel tempering
         └── Convergence diagnostics
```

## Command Line Interface

### Step 1: Glycosylation Preparation

Adds glycans to protein structure and optimizes asparagine side chain orientation to maximize HD22 distance from neighboring atoms.

```bash
glyco-prep -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input`: Input PDB file (protein structure without glycans)
- `-o, --output-dir`: Output directory for results

**Optional arguments:**

- `--asn-tsv`: TSV file with glycosylation sites containing columns: `site`, `iupac_glycosylator`, `residue_number`, `protein_chain`. If not provided, auto-detects from structure.

- `--rotate-atoms`: Comma-separated list of atoms to rotate during asparagine optimization. Rotation is performed around the CA-CB axis to maximize HD22 distance. Default: "OD1,CG,ND2,HD22,HD21,HB2,HB3"

- `--fixed-atom`: Atom that serves as the fixed point for rotation (default: "CB"). The rotation axis is defined by the vector from CA to this atom.

- `--center-atom`: Atom that defines the center of the neighbor detection sphere (default: "CA"). All atoms within the sphere radius are considered for distance calculation.

- `--radius`: Radius in Angstroms for neighbor detection during asparagine optimization. Atoms within this distance from the center atom are considered for clash detection (default: 30.0).

- `--rotation-step`: Rotation step in degrees for asparagine side chain optimization. The algorithm tests rotations from 0 to 360 degrees at this increment (default: 1).

- `--protein-residue-start`: Starting residue number for protein renumbering. Used to calculate chain offsets (default: 1).

- `--keep-temp`: Keep temporary files generated during processing.

### Step 2: Parametrization

Generates force field parameters and topology files for GROMACS simulations using CHARMM36 force field.

```bash
glyco-param -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input-pdb`: Input PDB file from Step 1 (glycosylated_protein_renumbered.pdb)
- `-o, --output-dir`: Output directory for topology files

**Optional arguments:**

- `--download-charmm`: Automatically download CHARMM36 force field from the official repository.

- `--charmm-url`: Custom URL for CHARMM36 force field download. Use this if the default mirror is unavailable.

- `--force-download`: Force re-download of CHARMM36 even if a local backup exists.

- `--n-cpus`: Number of CPU cores for parallel processing of glycan parameter generation (default: 1).

- `--keep-intermediate`: Keep intermediate files including individual glycan JSON, RTP, and HDB files.

### Step 3: Carbohydrate Orientation

Optimizes glycan orientations using a two-stage approach: single-axis grid search followed by MCMC refinement around the optimal axis.

```bash
glyco-orient -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input-pdb`: Input PDB file from Step 2
- `-o, --output-dir`: Output directory for optimized structures

**MCMC Parameters:**

- `--theta-step`: Step size in degrees for grid search over the rotation axis (0 to 360 degrees). Larger values give faster but coarser search (default: 10).

- `--n-steps`: Number of MCMC refinement steps per glycan per cycle. More steps improve convergence but increase computation time (default: 10000).

- `--max-cycles`: Maximum number of optimization cycles. The pipeline iterates until all glycans converge or this limit is reached (default: 5).

- `--radius`: Interaction radius in Angstroms for local energy calculation. Only atoms within this sphere centered on the glycan's center of mass are considered (default: 300.0).

- `--use-coulomb`: Include Coulomb electrostatic energy in the total energy calculation (default: false).

**Parallel Processing:**

- `--n-workers`: Number of CPU workers for parallel MCMC processing. Each worker processes different theta angles simultaneously during grid search (default: 4).

**Output Options:**

- `--save-individual-glycans`: Save each optimized glycan as a separate PDB file in the output directory.

- `--save-before-after`: Save PDB files before and after optimization for each glycan to compare orientations.

- `--verbose`: Print detailed progress information including energy values, acceptance rates, and coordinate changes.

- `--report-file`: Custom path for the optimization report file containing energy history and convergence data.

**Force Field Options:**

- `--charmm-dir`: Custom directory containing CHARMM36 force field files (atomtypes.atp, ffnonbonded.itp, and .rtp files).

### Complete Pipeline

Runs all three steps sequentially with coordinated parameters.

```bash
glyco-all -i INPUT_PDB -o OUTPUT_DIR [OPTIONS]
```

**Required arguments:**
- `-i, --input`: Input PDB file
- `-o, --output-dir`: Output directory for all results

**General Options:**

- `--asn-tsv`: TSV file with glycosylation sites

- `--download-charmm`: Download CHARMM36 force field automatically

- `--n-cpus`: CPUs for parametrization parallel processing (default: 1)

- `--n-workers`: Workers for MCMC parallel processing (default: 1)

- `--keep-temp`: Keep temporary files from all steps

- `--verbose`: Verbose output for all steps

**Asparagine Orientation Options (passed to Step 1):**

- `--rotate-atoms`: Atoms to rotate (default: "OD1,CG,ND2,HD22,HD21,HB2,HB3")

- `--fixed-atom`: Fixed atom for rotation (default: "CB")

- `--center-atom`: Center atom for neighbor sphere (default: "CA")

- `--radius`: Neighbor detection radius in Angstroms (default: 30.0)

- `--rotation-step`: Rotation step in degrees (default: 1)

**MCMC Options (passed to Step 3):**

- `--theta-step`: Grid search step size in degrees (default: 10)

- `--n-steps`: MCMC refinement steps per cycle (default: 10000)

- `--max-cycles`: Maximum optimization cycles (default: 5)

- `--mcmc-radius`: Interaction radius for energy calculation in Angstroms (default: 300.0)

- `--use-coulomb`: Include Coulomb electrostatics (default: false)

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
    --n-steps 20000 \
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
    --n-steps 50000 \
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
│   ├── glycosylated_protein_final_optimized_initial.pdb
│   ├── report.txt           # MCMC convergence report
│   └── PDB_CARB_ONLY/       # Individual optimized glycans
│       ├── [glycan].pdb
│       ├── [glycan]_before_[cycle].pdb
│       ├── [glycan]_after_[cycle].pdb
│       └── [glycan]_trajectory/
│           ├── step_*.pdb
│           └── energy_*.dat
```

## Mathematical Details

### Rotation Axis Definition

For each glycan, the rotation axis is defined as the vector from the protein attachment atom (ND2 for N-linked ASN, OG/OG1 for O-linked SER/THR) to the C1 carbon of the first glycan residue. Rotation is performed around this axis with the C1 atom as the pivot point.

### Grid Search

A full 360° grid search is performed at $\theta_{\text{step}}$ increments to find the orientation that minimizes:

$E(\theta) = E_{\text{vdW}}(\theta) + E_{\text{coulomb}}(\theta)$

### MCMC Refinement

After identifying the optimal angle $\theta_{\text{best}}$ from grid search, MCMC refinement explores the local energy landscape:

$\theta_{\text{proposed}} = \theta_{\text{current}} + \Delta\theta$, where $\Delta\theta \sim \mathcal{U}(-\theta_{\text{step}}, \theta_{\text{step}})$

Acceptance probability:

$\alpha = \min\left(1, \exp\left(-\frac{E(\theta_{\text{proposed}}) - E(\theta_{\text{current}})}{k_B T}\right)\right)$

### Convergence Criteria

A glycan is considered converged when the energy improvement between cycles is less than 1.0 kJ/mol.

### Detailed Balance Proof

The Metropolis-Hastings algorithm satisfies detailed balance:

$\pi(\theta) P(\theta \rightarrow \theta') = \pi(\theta') P(\theta' \rightarrow \theta)$

For our symmetric proposal:

$P(\theta \rightarrow \theta') = q(\theta'|\theta) \alpha(\theta,\theta') = \frac{1}{\sigma\sqrt{2\pi}} e^{-(\theta'-\theta)^2/2\sigma^2} \cdot \min(1, e^{-(E(\theta')-E(\theta))/k_B T})$

This ensures the stationary distribution is exactly $\pi(\theta)$.

## Troubleshooting

### Common Issues and Solutions

**Issue**: CHARMM force field download fails

**Solution**: Manually download using a mirror URL
```bash
glyco-param -i input.pdb -o output --charmm-url "your_mirror_url"
```

**Issue**: MCMC not converging (low acceptance rate)

**Solution**: Adjust step size or increase temperature
```bash
glyco-orient -i input.pdb -o output --n-steps 50000 --theta-step 5
```

**Issue**: High autocorrelation in MCMC samples

**Solution**: Increase thinning interval or run longer chains
```bash
glyco-orient -i input.pdb -o output --n-steps 100000 --max-cycles 10
```

**Issue**: Memory error during parametrization

**Solution**: Reduce parallel workers
```bash
glyco-param -i input.pdb -o output --n-cpus 2
```

**Issue**: Glycosylation site not found

**Solution**: Ensure TSV file has correct columns and residue numbers match the protein numbering

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{automated_glycosylation_2026,
  author = {Silva de Souza, Anacleto},
  title = {Automated Glycosylation Pipeline for Glycoproteins},
  year = {2026},
  url = {https://github.com/anacletosouza/automated_glycosylation}
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
