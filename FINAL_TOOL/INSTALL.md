# Installation and Setup Guide

## Prerequisites

1. **Python 3.7+** with the following packages:
   ```bash
   pip install pandas numpy scikit-learn scipy biopython
   ```

2. **BLAST+** (required for KEGG annotation)
   - Install BLAST+ from NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
   - Ensure `blastp` and `makeblastdb` are in your PATH

3. **KEGG Database** (required)
   - The KEGG metabolic pathways FASTA file should be located at:
     `/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/fasta/merged_metabolic_pathways.fasta`
   - Or update the path in `config.py`

4. **Model Weights** (required)
   - Model weight files should be in:
     `/groups/itay_mayrose/alongonda/desktop/mibig_validate/`
   - Should contain subdirectories: `promoter/`, `docking/`, `cladepp/`, `foldseek/`, `e2p2/`
   - Each should contain `*_multi_feature_weights.csv` files

## Optional Dependencies

For full feature extraction (currently placeholders return NaN if unavailable):

- **LightDock**: For docking features
- **CladePP**: For phylogenetic profiling features  
- **Foldseek**: For structural matching features
- **E2P2**: For enzyme classification features

## Setup Steps

1. **Clone or navigate to the FINAL_TOOL directory**
   ```bash
   cd /groups/itay_mayrose/alongonda/desktop/FINAL_TOOL
   ```

2. **Verify configuration**
   - Check `config.py` and update paths if needed:
     - `KEGG_DB`: Path to KEGG FASTA database
     - `MODEL_WEIGHTS_DIR`: Path to model weights directory

3. **Test installation**
   ```bash
   python3 mgc_predictor.py --help
   ```

## Usage Example

```bash
python3 mgc_predictor.py \
    genome_annotation.gff3 \
    genome_sequence.fasta \
    -o predictions.csv \
    --window-size 10 \
    --min-genes 3 \
    --threshold 0.5
```

## Troubleshooting

### Import Errors
- Ensure you're running from the correct directory
- Check that all Python dependencies are installed
- Verify Python version is 3.7+

### BLAST Errors
- Verify BLAST+ is installed: `which blastp`
- Check that KEGG database file exists and is readable
- Ensure sufficient disk space for temporary files

### Missing Model Weights
- Check that weight files exist in the configured directory
- Verify file names match expected patterns (see config.py)
- Some categories may work even if others are missing (will use NaN values)

### Memory Issues
- Large genomes may require significant memory
- Consider processing chromosomes separately if needed
- Adjust temporary directory location if disk space is limited


