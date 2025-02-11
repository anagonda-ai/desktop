#!/bin/bash
#SBATCH --job-name=statistical_scanner
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym

conda init
conda activate base

# Run your script
python /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/girvan_newman_graph_clustering.py
