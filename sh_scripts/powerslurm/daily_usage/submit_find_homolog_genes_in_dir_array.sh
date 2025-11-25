#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=itaym-pool
#SBATCH --time=6-24:00:00

# Array job script for homolog gene finding
# Usage: sbatch --array=1-N submit_find_homolog_genes_in_dir_array.sh mgc_list_file.txt

MGC_LIST_FILE=$1

# Get the MGC path for this array task
CLUSTRER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MGC_LIST_FILE")

conda run -n base python /groups/itay_mayrose/alongonda/desktop/python_scripts/asaph_aharoni/find_homolog_genes_for_asaph_cluster.py --example_mgc "$CLUSTRER"