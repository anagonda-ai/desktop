#!/bin/bash
#SBATCH --job-name=e2p2_tagging
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out_e2p2_tagging.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error_e2p2_tagging.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=7-00:00:00

python /groups/itay_mayrose/alongonda/desktop/python_scripts/features/e2p2_tagging/e2p2_tagging.py