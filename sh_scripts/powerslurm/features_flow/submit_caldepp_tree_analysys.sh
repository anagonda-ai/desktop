#!/bin/bash
#SBATCH --job-name=cladepp_tree_analysis
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out_cladepp_tree_analysis.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error_cladepp_tree_analysis.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym-pool
#SBATCH --time=7-00:00:00

tree_path="/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/species.nwk"
# mgc_candidates_dir="/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_dir_fixed"
mgc_candidates_dir="/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/test_dir"
mapping_file="/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping.csv"

python /groups/itay_mayrose/alongonda/desktop/python_scripts/features/cladepp_phylo_profiling/run_analysis.py --tree_path "$tree_path" --mgc_candidates_dir "$mgc_candidates_dir" --mapping_file "$mapping_file" --compute_gain_loss_coevolution