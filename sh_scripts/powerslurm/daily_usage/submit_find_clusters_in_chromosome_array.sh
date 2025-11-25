#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym-pool
#SBATCH --time=1-00:00:00

# Array job script for cluster finding
# Usage: sbatch --array=1-N submit_find_clusters_in_chromosome_array.sh mgc_list_file.txt

MGC_LIST_FILE=$1

# Get the MGC path for this array task
CLUSTRER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MGC_LIST_FILE")

conda run -n base python /groups/itay_mayrose/alongonda/desktop/python_scripts/asaph_aharoni/find_cluster_in_chromosomes.py --example_mgc "$CLUSTRER"