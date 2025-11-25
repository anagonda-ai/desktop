#!/bin/bash
#SBATCH --job-name=comparison_csv_pipeline
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/comparison_csv_pipeline.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/comparison_csv_pipeline.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym-pool
#SBATCH --time=6-24:00:00

conda run -n base python /groups/itay_mayrose/alongonda/desktop/python_scripts/asaph_aharoni/comparison_csv_pipeline.py