import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import partial
import pandas as pd
import subprocess
import threading
import argparse
from Bio import Phylo
import numpy as np

# Add the parent directory to Python path to find the modules
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

# Now try to import the modules with different approaches
try:
    # First try: assume we're running from the parent directory
    from cladepp_phylo_profiling_helpfuncs.cladepp_core import normalize_npp, build_profile_matrix, compute_gain_loss_coevolution_copap_style
    from cladepp_phylo_profiling_helpfuncs.io_utils import load_selected_blast_results, load_mapping_if_exists
except ModuleNotFoundError:
    try:
        # Second try: import directly from current directory
        from cladepp_core import normalize_npp, build_profile_matrix, compute_gain_loss_coevolution_copap_style
        from io_utils import load_selected_blast_results, load_mapping_if_exists
    except ModuleNotFoundError:
        # Third try: add current directory to path and import
        sys.path.insert(0, script_dir)
        from cladepp_core import normalize_npp, build_profile_matrix, compute_gain_loss_coevolution_copap_style
        from io_utils import load_selected_blast_results, load_mapping_if_exists


def save_clade_heatmap(npp_matrix, clade_id, tip_names, output_dir):
    # Thread-safe matplotlib usage
    plt.figure(figsize=(5, 4))
    corr_matrix = npp_matrix.T.corr()
    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1)
    plt.title(f"Clade {clade_id} (n={len(tip_names)}): {', '.join(tip_names[:3])}...")
    os.makedirs(output_dir, exist_ok=True)
    plt.tight_layout()
    
    # Thread-safe file saving
    with threading.Lock():
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

def normalize_name(name):
    return name.strip().lower().replace(" ", "_").replace("-", "_").replace(".", "").replace("(", "").replace(")", "")

def process_single_clade(clade_data, output_dir):
    """
    Process a single clade - designed for concurrent execution.
    """
    clade_id, tip_names, tree_name_map, anchor_genes, heatmap_dir, compute_gain_loss_coevolution, tree = clade_data
    
    try:
        sub_comparisons = [tree_name_map[name] for name in tip_names]
        sub_df = pd.DataFrame(sub_comparisons)

        blast_df = load_selected_blast_results(sub_df)
        raw_matrix = build_profile_matrix(blast_df, anchor_genes, output_dir)
        
        if raw_matrix.shape[1] < 2:  # Less than 2 organisms
            return {
                "status": "skipped", 
                "clade_id": clade_id,
                "reason": f"only {raw_matrix.shape[1]} organism(s)"
            }
        
        npp_matrix = normalize_npp(raw_matrix, output_dir)
        corr_stats, anchor_corrs = compute_anchor_corr_stats(npp_matrix, anchor_genes=anchor_genes)
        cladepp_score = compute_cladepp_score(npp_matrix, anchor_genes)

        # Gain/Loss Coevolution if requested
        gain_loss_score = None
        if compute_gain_loss_coevolution and tree is not None:
            gain_loss_scores = compute_gain_loss_coevolution_copap_style(raw_matrix, tree)
            gain_loss_score = np.mean(list(gain_loss_scores.values()))

        # Save heatmap
        save_clade_heatmap(npp_matrix, clade_id, tip_names, heatmap_dir)

        return {
            "status": "success",
            "clade_id": clade_id,
            "clade_size": len(tip_names),
            "tip_names": ",".join(tip_names[:5]) + "..." if len(tip_names) > 5 else ",".join(tip_names),
            **corr_stats,
            "cladepp_score": cladepp_score,
            "gain_loss_score": gain_loss_score
        }

    except Exception as e:
        return {
            "status": "error", 
            "clade_id": clade_id,
            "error": str(e)
        }

def match_tree_to_comparison(tree_tips, comparison_df, mapping_df):
    tree_names = {normalize_name(name): name for name in tree_tips}
    matched_rows = []
    
    for i, row in comparison_df.iterrows():
        organism = row["Directory"].split("/")[-2]
        print(organism)
        mapping_match = mapping_df[mapping_df["Original Filename"].str.contains(organism)]
        if not mapping_match.empty:
            organism_name = mapping_match.iloc[0]["Organism"]
            norm_org_name = normalize_name(organism_name)
            if norm_org_name in tree_names:
                matched_rows.append((row.to_dict(), tree_names[norm_org_name]))
    
    return matched_rows

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

def analyze_tree_clades_dynamic(
    tree_path, 
    comparison_csv, 
    mgc_output_dir, 
    mapping_file, 
    compute_gain_loss_coevolution,
    max_workers=8,
    use_processes=False
):
    """
    Analyze tree clades with concurrent processing.
    
    Parameters:
    - max_workers: Number of concurrent workers (default: CPU count)
    - use_processes: If True, use ProcessPoolExecutor; if False, use ThreadPoolExecutor
    """
    os.makedirs(mgc_output_dir, exist_ok=True)
    
    # Find fasta files in the comparison_csv dir and use them as anchor genes
    anchor_genes = get_anchor_genes_from_comparison_dir(comparison_csv)
    print(f"Loading tree from {tree_path}")
    tree = Phylo.read(tree_path, "newick")
    terminals = [term.name for term in tree.get_terminals()]

    print("Loading comparison table")
    comp_df = pd.read_csv(comparison_csv)
    print(comp_df)
    comp_df = comp_df.dropna(subset=["Directory", "Largest Chromosome File"])

    mapping_df = load_mapping_if_exists(mapping_file)
    match_list = match_tree_to_comparison(terminals, comp_df, mapping_df)

    print(f"Matched {len(match_list)} entries from comparison table to tree")

    tree_name_map = {t[1]: t[0] for t in match_list}
    heatmap_dir = os.path.join(mgc_output_dir, "clade_figures")
    os.makedirs(heatmap_dir, exist_ok=True)

    # 🟦 Compute Global CladePP before clade loop
    try:
        print("Computing global CladePP score...")
        all_matched_comparisons = [tree_name_map[name] for name in terminals if name in tree_name_map]
        all_df = pd.DataFrame(all_matched_comparisons)

        blast_df_global = load_selected_blast_results(all_df)
        raw_matrix_global = build_profile_matrix(blast_df_global, anchor_genes, mgc_output_dir)
        npp_matrix_global = normalize_npp(raw_matrix_global)

        global_score = compute_cladepp_global_score(npp_matrix_global, anchor_genes)
        print(f"✅ Global CladePP Score: {global_score}")

        # Save matrix and score
        npp_matrix_global.to_csv(os.path.join(mgc_output_dir, "matrix_npp_global.csv"))
        raw_matrix_global.to_csv(os.path.join(mgc_output_dir, "matrix_raw_global.csv"))
        with open(os.path.join(mgc_output_dir, "global_score.txt"), "w") as f:
            f.write(f"CladePP Global Score: {global_score}\n")
    except Exception as e:
        print(f"⚠️ Skipping global CladePP due to error: {e}")
        global_score = None

    # 🔵 Prepare clade data for concurrent processing
    print("Preparing clade data for concurrent analysis...")
    clade_tasks = []
    clade_id = 1
    
    for clade in tree.get_nonterminals():
        tips = clade.get_terminals()
        tip_names = [t.name for t in tips if t.name in tree_name_map]
        if len(tip_names) < 3:
            continue
            
        clade_data = (
            clade_id, tip_names, tree_name_map, anchor_genes, 
            heatmap_dir, compute_gain_loss_coevolution, 
            tree if compute_gain_loss_coevolution else None
        )
        clade_tasks.append(clade_data)
        clade_id += 1

    print(f"Processing {len(clade_tasks)} clades concurrently...")
    
    # Choose executor type
    executor_class = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    
    result_rows = []
    completed_count = 0
    
    # Process clades concurrently
    with executor_class(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_clade = {
            executor.submit(process_single_clade, clade_data, mgc_output_dir): clade_data[0] 
            for clade_data in clade_tasks
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_clade):
            clade_id = future_to_clade[future]
            completed_count += 1
            
            try:
                result = future.result()
                
                if result["status"] == "success":
                    result_rows.append(result)
                    print(f"✅ Completed clade {result['clade_id']} ({completed_count}/{len(clade_tasks)})")
                elif result["status"] == "skipped":
                    print(f"⏭️  Skipped clade {result['clade_id']}: {result['reason']} ({completed_count}/{len(clade_tasks)})")
                else:  # error
                    print(f"❌ Failed clade {result['clade_id']}: {result['error']} ({completed_count}/{len(clade_tasks)})")
                    
            except Exception as exc:
                print(f"❌ Clade {clade_id} generated an exception: {exc} ({completed_count}/{len(clade_tasks)})")

    # Save results
    if result_rows:
        df = pd.DataFrame(result_rows)
        df.to_csv(os.path.join(mgc_output_dir, f"summary.csv"), index=False)
        print(f"✅ Saved {len(result_rows)} successful clade analyses to summary.csv")
    else:
        print("⚠️ No successful clade analyses to save")
    
    print(f"🎉 Analysis complete! Results saved in {mgc_output_dir}")
    
def main():
    parser = argparse.ArgumentParser(
        description="Analyze tree clades for phylogenetic profiling",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--tree_path", 
        required=True, 
        help="Path to the phylogenetic tree file"
    )
    
    parser.add_argument(
        "--comparison_csv", 
        required=True, 
        help="Path to the comparison CSV file"
    )
    
    parser.add_argument(
        "--mgc_output_dir", 
        required=True, 
        help="Output directory for this MGC analysis"
    )
    
    parser.add_argument(
        "--mapping_file", 
        required=True, 
        help="Path to organism mapping file"
    )
    
    parser.add_argument(
        "--compute_gain_loss_coevolution", 
        type=str, 
        choices=['true', 'false'], 
        required=True, 
        help="Whether to compute gain/loss coevolution (true/false)"
    )
    
    parser.add_argument(
        "--max_workers", 
        type=int, 
        default=8, 
        help="Number of concurrent workers"
    )
    
    parser.add_argument(
        "--use_processes", 
        type=str, 
        choices=['true', 'false'], 
        default='false', 
        help="Whether to use processes instead of threads (true/false)"
    )
    
    args = parser.parse_args()
    
    # Convert string boolean arguments to actual booleans
    compute_gain_loss = args.compute_gain_loss_coevolution.lower() == 'true'
    use_processes = args.use_processes.lower() == 'true'
    
    # Print configuration
    print(f"Starting tree clade analysis with the following parameters:")
    print(f"  Tree path: {args.tree_path}")
    print(f"  Comparison CSV: {args.comparison_csv}")
    print(f"  Output directory: {args.mgc_output_dir}")
    print(f"  Mapping file: {args.mapping_file}")
    print(f"  Compute gain/loss coevolution: {compute_gain_loss}")
    print(f"  Max workers: {args.max_workers}")
    print(f"  Use processes: {use_processes}")
    print(f"  MGC name: {os.path.basename(args.mgc_output_dir)}")
    print()
    
    try:
        # Call the dynamic analysis function
        analyze_tree_clades_dynamic(
            tree_path=args.tree_path,
            comparison_csv=args.comparison_csv,
            mgc_output_dir=args.mgc_output_dir,
            mapping_file=args.mapping_file,
            compute_gain_loss_coevolution=compute_gain_loss,
            max_workers=args.max_workers,
            use_processes=use_processes
        )
        
        print(f"✅ Successfully completed analysis for {os.path.basename(args.mgc_output_dir)}")
        
    except Exception as e:
        print(f"❌ Failed analysis for {os.path.basename(args.mgc_output_dir)}: {str(e)}")
        
        # Create error log
        mgc_name = os.path.basename(args.mgc_output_dir)
        error_log_path = os.path.join(args.mgc_output_dir, f"{mgc_name}_error.log")
        os.makedirs(os.path.dirname(error_log_path), exist_ok=True)
        
        with open(error_log_path, 'w') as f:
            f.write(f"Error analyzing {mgc_name}:\n")
            f.write(f"Error: {str(e)}\n")
            f.write(f"Comparison CSV: {args.comparison_csv}\n")
            f.write(f"Tree path: {args.tree_path}\n")
            f.write(f"Mapping file: {args.mapping_file}\n")
            f.write(f"Compute gain/loss coevolution: {compute_gain_loss}\n")
            f.write(f"Max workers: {args.max_workers}\n")
            f.write(f"Use processes: {use_processes}\n")
        
        sys.exit(1)

if __name__ == "__main__":
    main()