import argparse
from cladepp_phylo_profiling_helpfuncs.tree_analysis import analyze_tree_clades_dynamic

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CladePP analysis with optional gain/loss coevolution")
    parser.add_argument("--tree-path", required=True, help="Path to species tree in Newick format")
    parser.add_argument("--comparison_csv", required=True, help="Path to comparison_results.csv")
    parser.add_argument("--output_prefix", default="clade_analysis", help="Prefix for output files")
    parser.add_argument("--mapping_file", default=None, help="Optional mapping CSV")
    parser.add_argument("--compute_gain_loss_coevolution", action="store_true", help="Compute gain/loss coevolution scores")

    args = parser.parse_args()

    analyze_tree_clades_dynamic(
        tree_path=args.tree_path,
        comparison_csv=args.comparison_csv,
        output_prefix=args.output_prefix,
        mapping_file=args.mapping_file,
        compute_gain_loss_coevolution=args.compute_gain_loss_coevolution
    )
