# Conda Environment Setup for MGC Predictor

This document provides instructions for setting up a conda environment for the MGC Predictor tool.

## Quick Start

```bash
# Create environment from environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate mgc_predictor

# Install BLAST+ (required for KEGG annotation)
conda install -c bioconda blast
```

## Manual Setup

If you prefer to set up the environment manually:

```bash
# Create a new conda environment with Python
conda create -n mgc_predictor python=3.10

# Activate the environment
conda activate mgc_predictor

# Install core dependencies from conda-forge
conda install -c conda-forge pandas numpy scipy biopython scikit-learn joblib xgboost

# Install BLAST+ from bioconda (required for KEGG annotation)
conda install -c bioconda blast

# Or install from requirements.txt using pip (alternative)
pip install -r requirements.txt
# Then install BLAST+ separately: conda install -c bioconda blast
```

## Alternative: pip-only Installation

If you prefer to use pip only:

```bash
# Create a new conda environment
conda create -n mgc_predictor python=3.10

# Activate the environment
conda activate mgc_predictor

# Install all Python dependencies
pip install -r requirements.txt

# Install BLAST+ via conda (recommended) or system package manager
conda install -c bioconda blast
```

## Verify Installation

After installation, verify that all dependencies are available:

```bash
# Activate environment
conda activate mgc_predictor

# Verify Python packages
python -c "import pandas, numpy, scipy, Bio, xgboost, sklearn, joblib; print('All Python packages installed successfully')"

# Verify BLAST+
blastp -version
```

## External Dependencies

The following external tools are optional but may be needed for full feature extraction:

- **LightDock**: For docking feature extraction (if using pre-computed docking results)
- **Foldseek**: For structural matching feature extraction (if using pre-computed Foldseek results)
- **CladePP**: For phylogenetic profiling feature extraction (if using pre-computed CladePP results)
- **E2P2**: For enzyme classification feature extraction (if using pre-computed E2P2 results)

These tools are not required if you are using pre-computed feature results files, as the tool can parse existing result files directly.

## Troubleshooting

### BLAST+ Installation Issues

If you have trouble installing BLAST+:

1. **Via conda**: `conda install -c bioconda blast`
2. **Via system package manager**:
   - Ubuntu/Debian: `sudo apt-get install ncbi-blast+`
   - macOS: `brew install blast`
   - CentOS/RHEL: `sudo yum install ncbi-blast+`
3. **Manual download**: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

### XGBoost Installation Issues

If XGBoost fails to install:

```bash
# Try installing from conda-forge first
conda install -c conda-forge xgboost

# If that fails, try pip
pip install xgboost
```

### Biopython Installation Issues

```bash
# Install from conda-forge (recommended)
conda install -c conda-forge biopython

# Or via pip
pip install biopython
```

