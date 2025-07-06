import pandas as pd
from Bio import Phylo
import numpy as np
from itertools import combinations
import os
import matplotlib.pyplot as plt
import seaborn as sns
from cladepp_phylo_profiling_helpfuncs.cladepp_core import normalize_npp, build_profile_matrix, compute_gain_loss_coevolution_copap_style
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


def compute_anchor_corr_stats(submatrix, anchor_genes=None):
    genes = [g for g in anchor_genes if g in submatrix.index]
    anchor_corrs = []
    corr_gene_pairs = []
    for g1, g2 in combinations(genes, 2):
        corr = submatrix.loc[g1].corr(submatrix.loc[g2])
        anchor_corrs.append(corr)
        if corr > 0:
            corr_gene_pairs.append((g1, g2, corr))
    stats = {
        "mean_anchor_corr": np.mean(anchor_corrs) if anchor_corrs else np.nan,
        "std_anchor_corr": np.std(anchor_corrs) if anchor_corrs else np.nan,
        "max_anchor_corr": np.max(anchor_corrs) if anchor_corrs else np.nan,
        "min_anchor_corr": np.min(anchor_corrs) if anchor_corrs else np.nan,
        "positive_corr_pairs": "; ".join(f"{g1}-{g2}:{corr:.2f}" for g1, g2, corr in corr_gene_pairs)
    }
    return stats, anchor_corrs


def save_clade_heatmap(npp_matrix, clade_id, tip_names, output_dir):
    plt.figure(figsize=(5, 4))
    corr_matrix = npp_matrix.T.corr()
    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1)
    plt.title(f"Clade {clade_id} (n={len(tip_names)}): {', '.join(tip_names[:3])}...")
    os.makedirs(output_dir, exist_ok=True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"clade_{clade_id}_corr_heatmap.png"))
    plt.close()


def compute_cladepp_score(npp_matrix, anchor_genes):
    """
    Compute CladePP score based on pairwise correlation between anchor genes.
    """
    from itertools import combinations
    existing_genes = [g for g in anchor_genes if g in npp_matrix.index]

    corr_matrix = npp_matrix.T.corr()

    corr_values = []
    for g1, g2 in combinations(existing_genes, 2):
        corr = corr_matrix.loc[g1, g2]
        corr_values.append(corr)

    return np.mean(corr_values) if corr_values else np.nan


def compute_cladepp_global_score(npp_matrix: pd.DataFrame, anchor_genes: list[str]) -> float:
    existing_genes = [g for g in anchor_genes if g in npp_matrix.index]

    scores = []
    for g1, g2 in combinations(existing_genes, 2):
        corr = npp_matrix.loc[g1].corr(npp_matrix.loc[g2])
        scores.append(corr)
    
    return float(np.mean(scores)) if scores else np.nan

def get_anchor_genes_from_comparison_dir(comparison_csv):
    """
    Extract anchor genes from the directory of the comparison CSV file.
    Assumes that anchor genes are fasta files in the same directory.
    """
    comp_dir = os.path.dirname(comparison_csv)
    anchor_genes = [
        os.path.splitext(f)[0]
        for f in os.listdir(comp_dir)
        if f.endswith(".fasta") or f.endswith(".fa") or f.endswith(".faa")
    ]
    return anchor_genes

def analyze_tree_clades_dynamic(
    tree_path, 
    comparison_csv, 
    output_prefix="clade_analysis", 
    mapping_file=None, 
    compute_gain_loss_coevolution=False
):
    # Find fasta files in the comparison_csv dir and use them as anchor genes
    anchor_genes = get_anchor_genes_from_comparison_dir(comparison_csv)
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
    heatmap_dir = f"{output_prefix}_clade_figures"
    os.makedirs(heatmap_dir, exist_ok=True)

    # 🟦 Compute Global CladePP before clade loop
    try:
        all_matched_comparisons = [tree_name_map[name] for name in terminals if name in tree_name_map]
        all_df = pd.DataFrame(all_matched_comparisons).reset_index(drop=True)

        blast_df_global = load_selected_blast_results(all_df)
        raw_matrix_global = build_profile_matrix(blast_df_global, anchor_genes)
        npp_matrix_global = normalize_npp(raw_matrix_global)

        global_score = compute_cladepp_global_score(npp_matrix_global, anchor_genes)
        print(f"✅ Global CladePP Score: {global_score}")

        # שמור את המטריצה והציון
        npp_matrix_global.to_csv(f"{output_prefix}_matrix_npp_global.csv")
        raw_matrix_global.to_csv(f"{output_prefix}_matrix_raw_global.csv")
        with open(f"{output_prefix}_global_score.txt", "w") as f:
            f.write(f"CladePP Global Score: {global_score}\n")
    except Exception as e:
        print(f"⚠️ Skipping global CladePP due to error: {e}")
        global_score = None

    # 🔵 Clade analysis
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
            corr_stats, anchor_corrs = compute_anchor_corr_stats(npp_matrix, anchor_genes=anchor_genes)
            cladepp_score = compute_cladepp_score(npp_matrix, anchor_genes)

            # Gain/Loss Coevolution if requested
            gain_loss_score = None
            if compute_gain_loss_coevolution:
                gain_loss_scores = compute_gain_loss_coevolution_copap_style(raw_matrix, tree)
                gain_loss_score = np.mean(list(gain_loss_scores.values()))

            result_rows.append({
                "clade_id": clade_id,
                "clade_size": len(tip_names),
                "tip_names": ",".join(tip_names[:5]) + "..." if len(tip_names) > 5 else ",".join(tip_names),
                **corr_stats,
                "cladepp_score": cladepp_score,
                "gain_loss_score": gain_loss_score
            })

            save_clade_heatmap(npp_matrix, clade_id, tip_names, heatmap_dir)

        except Exception as e:
            print(f"Skipping clade {clade_id} due to error: {e}")
        clade_id += 1

    df = pd.DataFrame(result_rows)
    df.to_csv(f"{output_prefix}_summary.csv", index=False)
    print(f"Saved clade correlation summary to {output_prefix}_summary.csv")

