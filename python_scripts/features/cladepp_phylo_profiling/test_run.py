import os
from cladepp_phylo_profiling_helpfuncs.analyze_tree_clades_dynamic import analyze_tree_clades_dynamic

TREE_PATH = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/species.nwk"
COMPARISON_CSV = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_final/MGC_CANDIDATE_1546/comparison_results.csv"
MAPPING_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping.csv"
OUTPUT_PREFIX = "/groups/itay_mayrose/alongonda/desktop/tree_test/clade_analysis_test"

os.makedirs(os.path.dirname(OUTPUT_PREFIX), exist_ok=True)

analyze_tree_clades_dynamic(
    tree_path=TREE_PATH,
    comparison_csv=COMPARISON_CSV,
    mgc_output_dir=OUTPUT_PREFIX,
    mapping_file=MAPPING_FILE,
    compute_gain_loss_coevolution=False,
    max_workers=32,
    use_processes=False
)

output_file = os.path.join(OUTPUT_PREFIX, "global_score.txt")
assert os.path.exists(output_file), "Output file not created"
with open(output_file) as f:
    lines = f.readlines()
    assert len(lines) > 1, "Output file is empty"
print("✅ Test passed")
