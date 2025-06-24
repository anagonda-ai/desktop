#!/bin/bash
#SBATCH --job-name=duplication_filter
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=3-24:00:00

python /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_processing/duplication_filter.py