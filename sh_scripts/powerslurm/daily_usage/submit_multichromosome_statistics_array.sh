#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=6-24:00:00

# Array job script for multichromosome statistics
# Usage: sbatch --array=1-N submit_multichromosome_statistics_array.sh mgc_list_file.txt

MGC_LIST_FILE=$1

# Get the MGC path for this array task
CLUSTRER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MGC_LIST_FILE")

python /groups/itay_mayrose/alongonda/desktop/python_scripts/asaph_aharoni/single_multi_chromosome_cluster_statistics.py --example_mgc "$CLUSTRER"