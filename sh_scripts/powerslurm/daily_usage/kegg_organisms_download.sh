#!/bin/bash
#SBATCH --job-name=kegg_organism_download
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

# Get the input FASTA file from the command line argument
CODE=$1

# Define the output directory
ROOT_FOLDER=$2

METABOLIC_MODULES_PATH=$3

python /groups/itay_mayrose/alongonda/desktop/python_scripts/process_organism.py "$CODE" "$ROOT_FOLDER" "$METABOLIC_MODULES_PATH"