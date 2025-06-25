import argparse
from cladepp_phylo_profiling_helpfuncs.tree_analysis import analyze_tree_clades_dynamic

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree-path", required=True)
    parser.add_argument("--comparison_csv", required=True)
    parser.add_argument("--anchor_genes", nargs="+", required=True)
    parser.add_argument("--output_prefix", default="clade_analysis")
    parser.add_argument("--mapping_file", default=None)
    args = parser.parse_args()

    analyze_tree_clades_dynamic(
        tree_path=args.tree_path,
        comparison_csv=args.comparison_csv,
        anchor_genes=args.anchor_genes,
        output_prefix=args.output_prefix,
        mapping_file=args.mapping_file
    )
