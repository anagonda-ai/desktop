#!/bin/bash
#SBATCH --job-name=kegg_organism_download
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00


python /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/sliding_windows/kegg_scanner_min_genes_based.py