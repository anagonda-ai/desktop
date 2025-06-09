#!/bin/bash
#SBATCH --job-name=colabfold_fixed
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

conda deactivate
conda activate colabfold_fixed

python /groups/itay_mayrose/alongonda/desktop/python_scripts/features/foldseek_matching/alphafold_predictions_candidates.py