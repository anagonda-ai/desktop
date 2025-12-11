import argparse
from cladepp_phylo_profiling_helpfuncs.tree_analysis import analyze_tree_clades

if __name__ == "__main__":

    analyze_tree_clades(
        tree_path="/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/species.nwk",
        mgc_candidates_dir="/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_dir_fixed",
        mapping_file="/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping.csv",
        compute_gain_loss_coevolution=False
    )
