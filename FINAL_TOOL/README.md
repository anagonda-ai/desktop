# MGC Predictor - Final Tool

A comprehensive tool for predicting Metabolic Gene Clusters (MGCs) from genome annotations.

## Overview

This tool processes GFF3 annotation files and corresponding genome FASTA files to:
1. Parse and order genes by position
2. Annotate genes using KEGG metabolic database (via BLAST)
3. Identify potential MGC candidates using sliding windows
4. Extract features (promoter, docking, CladePP, Foldseek, E2P2)
5. Predict MGC probability using trained machine learning models

## Installation

### Dependencies

- Python 3.7+
- Required Python packages:
  - pandas
  - numpy
  - scikit-learn
  - scipy
  - BioPython
- External tools (optional, for full feature extraction):
  - BLAST+ (required for KEGG annotation)
  - LightDock (for docking features)
  - CladePP (for phylogenetic profiling features)
  - Foldseek (for structural matching features)
  - E2P2 (for enzyme classification features)

### Setup

1. **Run the setup script** to verify and prepare all required data:
   ```bash
   python FINAL_TOOL/setup.py
   ```
   This will check:
   - BLAST+ installation
   - KEGG FASTA database (in `FINAL_TOOL/data/`)
   - Model files (in `FINAL_TOOL/models/`)
   
   All required data files are already included in the package.

2. **Alternative manual setup**:
   - Install dependencies: `conda env create -f environment.yml` or `pip install -r requirements.txt`
   - Install BLAST+: `conda install -c bioconda blast`
   - All data files should be in `FINAL_TOOL/data/` and `FINAL_TOOL/models/`

## Usage

### Basic Usage

```bash
python mgc_predictor.py <gff3_file> <fasta_file> -o <output_csv>
```

### Example

```bash
python mgc_predictor.py genome.gff3 genome.fasta -o results.csv
```

### Advanced Options

```bash
python mgc_predictor.py genome.gff3 genome.fasta \
    -o results.csv \
    --window-size 15 \
    --min-genes 4 \
    --threshold 0.6 \
    --temp-dir /tmp/mgc_work \
    --kegg-db /path/to/kegg.fasta
```

### Arguments

- `gff3_file`: Path to GFF3 annotation file
- `fasta_file`: Path to genome FASTA file
- `-o, --output`: Output CSV file path (required)
- `--window-size`: Sliding window size in genes (default: 10)
- `--min-genes`: Minimum metabolic genes per window (default: 3)
- `--threshold`: Prediction threshold for MGC classification (default: 0.5)
- `--temp-dir`: Temporary directory for intermediate files
- `--kegg-db`: Path to KEGG metabolic pathways FASTA database

## Output

The output CSV file contains:
- Candidate information: candidate_id, chromosome, pathway, genes, coordinates
- Feature values: All extracted features (promoter, docking, CladePP, Foldseek, E2P2)
- Predictions: 
  - Category probabilities (promoter_probability, docking_probability, etc.)
  - Combined probability (combined_probability)
  - MGC prediction (is_mgc: True/False)

## Notes

- Some feature extractors (docking, CladePP, Foldseek, E2P2) require external tools. If these are not available, the corresponding features will be NaN, but predictions will still be made using available features.
- The promoter feature extractor uses a simplified k-mer based approach. Full implementation would use TFBS pattern databases.
- Model weights are loaded from CSV files. Intercepts may need to be loaded from model summary files for full accuracy.

## Directory Structure

```
FINAL_TOOL/
├── mgc_predictor.py          # Main entry point
├── config.py                 # Configuration
├── parsers/                  # GFF3 parsing
├── annotation/               # KEGG annotation
├── candidate_finder/         # Sliding window candidate finding
├── feature_extraction/       # Feature extractors
├── prediction/               # Model prediction
└── utils/                    # Utility functions
```


