import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple
from scipy.stats import mannwhitneyu

# === Configuration === #
BASE = "/groups/itay_mayrose/alongonda"
RESULTS_DIR = os.path.join(BASE, "Plant_MGC/results/analysis_metabolic_enrichment")
os.makedirs(RESULTS_DIR, exist_ok=True)

PREDICTED_PATH = os.path.join(BASE, "Plant_MGC/kegg_metabolic_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10.csv")
RANDOM_PATH = os.path.join(BASE, "Plant_MGC/kegg_metabolic_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/random_control_windows.csv")
MIBIG_SPLIT_DIR = os.path.join(BASE, "Plant_MGC/mibig_metabolic_output/annotated_genomes_metabolic/merged_annotation_with_chloroplast/split_by_mgc")
MIBIG_STATS_PATH = os.path.join(BASE, "Plant_MGC/mibig_metabolic_output/annotated_genomes_metabolic/merged_annotation_with_chloroplast/mgc_best_sources_stats.csv")

# === Functions === #

def compute_ratio_kegg(df: pd.DataFrame, label: str) -> pd.DataFrame:
    df = df.copy()
    df["total_genes"] = df["genes"].apply(lambda x: len(str(x).split(",")))
    df["metabolic_genes_count"] = df["metabolic_genes"].apply(lambda x: len(str(x).split(",")) if pd.notna(x) and x else 0)
    df["metabolic_ratio"] = df["metabolic_genes_count"] / df["total_genes"]
    df["cluster_id"] = df.index.astype(str)
    df["group"] = label
    return df[["cluster_id", "total_genes", "metabolic_genes_count", "metabolic_ratio", "group"]]

def compute_ratio_mibig(split_dir: str, stats_path: str) -> pd.DataFrame:
    stats_df = pd.read_csv(stats_path)
    all_records = []
    for mgc_id in stats_df["mgc_id"]:
        mgc_file = os.path.join(split_dir, f"{mgc_id}.csv")
        if not os.path.exists(mgc_file):
            continue
        df = pd.read_csv(mgc_file)
        metabolic_count = df["matched_gene"].notna().sum()
        total_genes = df["matched_gene"].count()
        ratio = metabolic_count / total_genes if total_genes > 0 else 0
        all_records.append({
            "cluster_id": mgc_id,
            "total_genes": total_genes,
            "metabolic_genes_count": metabolic_count,
            "metabolic_ratio": ratio,
            "group": "mibig"
        })
    return pd.DataFrame(all_records)

def compare_group_ratios(data: pd.DataFrame, pairs: List[Tuple[str, str]]) -> pd.DataFrame:
    results = []
    for g1, g2 in pairs:
        r1 = data[data["group"] == g1]["metabolic_ratio"]
        r2 = data[data["group"] == g2]["metabolic_ratio"]
        u, p = mannwhitneyu(r1, r2, alternative="greater")
        results.append({
            "comparison": f"{g1} > {g2}",
            "U_statistic": round(u, 2),
            "p_value": p
        })
    return pd.DataFrame(results)

def plot_ratio_distributions(df: pd.DataFrame, output_path: str):
    sns.set(style="whitegrid")
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df, x="group", y="metabolic_ratio", palette="Set2")
    plt.title("Metabolic Gene Ratio per Cluster Type", fontsize=14)
    plt.ylabel("Ratio of Metabolic Genes")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# === Main === #

def main():
    print("▶ Processing KEGG predicted...")
    predicted_raw = pd.read_csv(PREDICTED_PATH)
    predicted_df = compute_ratio_kegg(predicted_raw, "predicted")

    print("▶ Processing random windows...")
    random_raw = pd.read_csv(RANDOM_PATH)
    random_df = compute_ratio_kegg(random_raw, "random")

    print("▶ Processing MIBiG data...")
    mibig_df = compute_ratio_mibig(MIBIG_SPLIT_DIR, MIBIG_STATS_PATH)

    print("▶ Combining datasets...")
    combined = pd.concat([mibig_df, predicted_df, random_df], ignore_index=True)
    combined.to_csv(os.path.join(RESULTS_DIR, "all_metabolic_ratios.csv"), index=False)

    print("▶ Running statistical tests...")
    stats = compare_group_ratios(combined, [
        ("mibig", "random"),
        ("mibig", "predicted"),
        ("predicted", "random")
    ])
    stats.to_csv(os.path.join(RESULTS_DIR, "statistical_tests.csv"), index=False)

    print("▶ Plotting results...")
    plot_ratio_distributions(combined, os.path.join(RESULTS_DIR, "metabolic_ratio_boxplot.png"))

    print("✅ Done.")
    print(stats)

if __name__ == "__main__":
    main()
