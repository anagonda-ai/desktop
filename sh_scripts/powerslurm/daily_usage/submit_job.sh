#!/bin/bash
#SBATCH --job-name=kegg_random_cladepp_mgc_analysis
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/kegg_random_cladepp_mgc_analysis.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/kegg_random_cladepp_mgc_analysis.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym-pool
#SBATCH --time=6-24:00:00

conda run -n base python /groups/itay_mayrose/alongonda/desktop/python_scripts/features/cladepp_phylo_profiling/run_analysis.py