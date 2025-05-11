#!/bin/bash
#SBATCH --job-name=statistical_scanner
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

conda init
conda deactivate
conda activate base


python /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_processing/duplication_filter.py