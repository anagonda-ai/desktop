import pandas as pd
from Bio import Phylo
import numpy as np
from itertools import combinations
import os
import matplotlib.pyplot as plt
import seaborn as sns
from cladepp_phylo_profiling_helpfuncs.cladepp_core import normalize_npp, build_profile_matrix
from cladepp_phylo_profiling_helpfuncs.io_utils import load_selected_blast_results, load_mapping_if_exists

def normalize_name(name):
    return name.strip().lower().replace(" ", "_").replace("-", "_").replace(".", "").replace("(", "").replace(")", "")

def match_tree_to_comparison(tree_tips, comparison_df, mapping_df):
    tree_names = {normalize_name(name): name for name in tree_tips}
    matched_rows = []
    for _, row in comparison_df.iterrows():
        organism = row["Directory"].split("/")[-2]
        mapping_match = mapping_df[mapping_df["Original Filename"].str.contains(organism)]
        if not mapping_match.empty:
            organism_name = mapping_match.iloc[0]["Organism"]
            norm_org_name = normalize_name(organism_name)
            if norm_org_name in tree_names:
                matched_rows.append((row.to_dict(), tree_names[norm_org_name]))
    return matched_rows

def compute_anchor_corr_stats(submatrix):
    genes = submatrix.index
    anchor_corrs = []
    for g1, g2 in combinations(genes, 2):
        anchor_corrs.append(submatrix.loc[g1].corr(submatrix.loc[g2]))
    if anchor_corrs:
        return {
            "mean_anchor_corr": np.mean(anchor_corrs),
            "std_anchor_corr": np.std(anchor_corrs),
            "max_anchor_corr": np.max(anchor_corrs),
            "min_anchor_corr": np.min(anchor_corrs)
        }
    else:
        return {
            "mean_anchor_corr": np.nan,
            "std_anchor_corr": np.nan,
            "max_anchor_corr": np.nan,
            "min_anchor_corr": np.nan
        }

def plot_and_save_heatmap(matrix, title, filename):
    plt.figure(figsize=(8, 6))
    sns.heatmap(matrix, annot=True, cmap="coolwarm", center=0)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def analyze_tree_clades_dynamic(tree_path, comparison_csv, anchor_genes, output_prefix="clade_analysis", mapping_file=None):
    print(f"Loading tree from {tree_path}")
    tree = Phylo.read(tree_path, "newick")
    terminals = [term.name for term in tree.get_terminals()]

    print("Loading comparison table")
    comp_df = pd.read_csv(comparison_csv)
    comp_df = comp_df.dropna(subset=["Directory", "Largest Chromosome File"])

    mapping_df = load_mapping_if_exists(mapping_file)
    match_list = match_tree_to_comparison(terminals, comp_df, mapping_df)

    print(f"Matched {len(match_list)} entries from comparison table to tree")

    tree_name_map = {t[1]: t[0] for t in match_list}

    result_rows = []
    clade_id = 1
    os.makedirs(f"{output_prefix}_clades", exist_ok=True)

    for clade in tree.get_nonterminals():
        tips = clade.get_terminals()
        tip_names = [t.name for t in tips if t.name in tree_name_map]
        if len(tip_names) < 3:
            continue

        sub_comparisons = [tree_name_map[name] for name in tip_names]
        sub_df = pd.DataFrame(sub_comparisons).reset_index(drop=True)

        try:
            blast_df = load_selected_blast_results(sub_df)
            raw_matrix = build_profile_matrix(blast_df, anchor_genes)
            npp_matrix = normalize_npp(raw_matrix)
            corr_stats = compute_anchor_corr_stats(npp_matrix)

            result_rows.append({
                "clade_id": clade_id,
                "clade_size": len(tip_names),
                "tip_names": ",".join(tip_names[:5]) + "..." if len(tip_names) > 5 else ",".join(tip_names),
                **corr_stats
            })

            # Save matrix files
            raw_path = f"{output_prefix}_clades/matrix_raw_clade_{clade_id}.csv"
            npp_path = f"{output_prefix}_clades/matrix_npp_clade_{clade_id}.csv"
            raw_matrix.to_csv(raw_path)
            npp_matrix.to_csv(npp_path)

            # Save heatmap
            plot_and_save_heatmap(npp_matrix, f"Normalized Profile Clade {clade_id}", f"{output_prefix}_clades/heatmap_clade_{clade_id}.png")

        except Exception as e:
            print(f"Skipping clade {clade_id} due to error: {e}")
        clade_id += 1

    df = pd.DataFrame(result_rows)
    df.to_csv(f"{output_prefix}_summary.csv", index=False)
    print(f"Saved clade correlation summary to {output_prefix}_summary.csv")
