import argparse
from cladepp_phylo_profiling_helpfuncs.tree_analysis import analyze_tree_clades

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CladePP analysis with optional gain/loss coevolution")
    parser.add_argument("--tree_path", required=True, help="Path to species tree in Newick format")
    parser.add_argument("--mgc_candidates_dir", required=True, help="Path to mgc_candidates_dir")
    parser.add_argument("--mapping_file", default=None, help="Optional mapping CSV")
    parser.add_argument("--compute_gain_loss_coevolution", action="store_true", help="Compute gain/loss coevolution scores")

    args = parser.parse_args()

    analyze_tree_clades(
        tree_path=args.tree_path,
        mgc_candidates_dir=args.mgc_candidates_dir,
        mapping_file=args.mapping_file,
        compute_gain_loss_coevolution=args.compute_gain_loss_coevolution
    )
