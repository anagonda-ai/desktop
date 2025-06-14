import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, fcluster
from cladepp_phylo_profiling_helpfuncs.io_utils import load_selected_blast_results
from cladepp_phylo_profiling_helpfuncs.plotting import plot_dendrogram, plot_heatmap

def build_profile_matrix(blast_df, anchors, min_coverage=0.2):
    blast_df = blast_df[blast_df["origin_file"].isin(anchors)]
    pivot = blast_df.pivot_table(
        index="origin_file",
        columns="organism",
        values="bit_score",
        aggfunc="max"
    ).fillna(0)
    pivot[pivot < 24.6] = 0

    min_orgs = int(min_coverage * pivot.shape[1])
    filtered = pivot[(pivot > 0).sum(axis=1) >= min_orgs]
    filtered.to_csv("matrix_raw.csv")
    return filtered

def normalize_npp(matrix):
    npp = matrix.apply(zscore, axis=0).fillna(0)
    npp.to_csv("matrix_npp.csv")
    return npp

def compute_clusters(npp_matrix, threshold):
    linkage_matrix = linkage(npp_matrix, method="average", metric="correlation")
    cluster_ids = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return pd.Series(cluster_ids, index=npp_matrix.index), linkage_matrix

def coevolution_score(cluster_series, anchor_names):
    scores = {}
    clusters = cluster_series.groupby(cluster_series)
    for clust_id, genes in clusters:
        in_cluster = list(genes.index)
        n_anchors = len([g for g in in_cluster if g in anchor_names])
        if n_anchors >= 2:
            for g in in_cluster:
                scores[g] = n_anchors / len(in_cluster)
    return scores

def cladepp_from_comparison_table(comparison_csv, anchor_genes, output_file,
                                   threshold=1.0, min_coverage=0.2,
                                   n_threads=4, plot=False):
    print("Loading selected BLAST results...")
    blast_df = load_selected_blast_results(comparison_csv, n_threads)

    print("Building profile matrix...")
    matrix = build_profile_matrix(blast_df, anchor_genes, min_coverage)

    print("Normalizing (NPP)...")
    npp = normalize_npp(matrix)

    print("Clustering...")
    print("Anchor correlation matrix:")
    corr = npp.T.corr()
    print(corr.round(2))
    corr.to_csv("anchor_correlation.csv")
    clusters, linkage_matrix = compute_clusters(npp, threshold)

    print("Scoring coevolution...")
    scores = coevolution_score(clusters, anchor_genes)

    print("Saving output...")
    df_scores = pd.DataFrame.from_dict(scores, orient='index', columns=['coevolution_score'])
    df_scores.index.name = 'gene'
    df_scores.to_csv(output_file)

    if plot:
        print("Generating plots...")
        plot_dendrogram(linkage_matrix, labels=npp.index.tolist())
        plot_heatmap(npp)

    print("Done. Results saved to", output_file)
    print("Top scoring genes:")
    print(df_scores.sort_values(by="coevolution_score", ascending=False).head(10))
