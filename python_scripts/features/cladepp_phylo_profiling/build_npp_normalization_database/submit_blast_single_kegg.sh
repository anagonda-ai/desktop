#!/bin/bash
#SBATCH --job-name=blast_organisms_on_aglu_selected_1000
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/blast_organisms_on_aglu_selected_1000.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/blast_organisms_on_aglu_selected_1000.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym-pool
#SBATCH --time=6-24:00:00

python /groups/itay_mayrose/alongonda/desktop/python_scripts/features/cladepp_phylo_profiling/build_npp_normalization_database/blast_organisms_on_single_kegg.py --example_mgc /groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/aglu_selected_1000