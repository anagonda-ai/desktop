#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=3-24:00:00

# Get the Cluster directory code from main script
CLUSTRER=$1

python /groups/itay_mayrose/alongonda/desktop/python_scripts/asaph_aharoni/find_cluster_in_chromosomes.py --example_mgc "$CLUSTRER" 