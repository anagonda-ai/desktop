#!/bin/bash
#SBATCH --job-name=kegg_organism_download_${SLURM_ARRAY_TASK_ID}
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

conda deactivate
conda activate colabfold_fixed

# Get the input FASTA file from the command line argument
fasta_path=$1

# Define the output directory
pdb_out_dir=$2

colabfold_batch $fasta_path $pdb_out_dir --msa-mode single_sequence --num-relax 0 --num-recycle 1