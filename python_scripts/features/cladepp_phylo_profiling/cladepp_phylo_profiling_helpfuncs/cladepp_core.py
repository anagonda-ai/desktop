import pandas as pd
from itertools import combinations
import numpy as np
import os
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, fcluster

def build_profile_matrix(blast_df, anchors, output_dir, min_coverage=0.2):
    print(f"🔍 Input blast_df shape: {blast_df.shape}")
    print(f"🔍 Unique organisms in blast_df: {blast_df['organism'].nunique()}")
    print(f"🔍 Organism list: {blast_df['organism'].unique()[:10]}")  # First 10
    
    blast_df = blast_df[blast_df["origin_file"].isin(anchors)]
    print(f"🔍 After anchor filtering: {blast_df.shape}")
    print(f"🔍 Remaining organisms after anchor filter: {blast_df['organism'].nunique()}")
    
    pivot = blast_df.pivot_table(
        index="origin_file",
        columns="organism",
        values="bit_score",
        aggfunc="max"
    ).fillna(0)
    print(f"🔍 Pivot shape (genes x organisms): {pivot.shape}")
    
    pivot[pivot < 24.6] = 0
    print(f"🔍 After score threshold: {pivot.shape}")
    
    min_orgs = int(min_coverage * pivot.shape[1])
    print(f"🔍 Min organisms required per gene: {min_orgs} (of {pivot.shape[1]})")
    
    filtered = pivot[(pivot > 0).sum(axis=1) >= min_orgs]
    print(f"🔍 Final matrix shape after gene filtering: {filtered.shape}")
    print(f"🔍 Final organism count: {filtered.shape[1]}")
    
    filtered.to_csv(os.path.join(output_dir, "matrix_raw.csv"))
    return filtered

def normalize_npp(matrix, output_dir):
    from scipy.stats import zscore
    
    print(f"Input type: {type(matrix)}")
    print(f"Input shape: {matrix.shape}")
    
    # Check if we have enough organisms for correlation
    if matrix.shape[1] < 2:
        raise ValueError(f"Need at least 2 organisms for correlation, got {matrix.shape[1]}")
    
    # Apply zscore directly on the numpy array to avoid pandas issues
    # axis=1 means normalize each gene (row) across organisms (columns)
    normalized_values = zscore(matrix.values, axis=1, nan_policy='omit')
    
    # Create new DataFrame with same structure
    npp = pd.DataFrame(normalized_values, index=matrix.index, columns=matrix.columns)
    npp = npp.fillna(0)
    
    print(f"Final type: {type(npp)}")
    print(f"Final shape: {npp.shape}")
    
    npp.to_csv(os.path.join(output_dir, "matrix_npp.csv"))
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

from Bio import Phylo

def compute_gain_loss_coevolution_copap_style(presence_absence_matrix, tree):
    """
    Compute gain/loss coevolution score using a parsimony approximation on the tree.
    For each gene pair, checks if presence/absence patterns map similarly to the tree.
    """
    def reconstruct_parsimony_states(gene_vector):
        states = {}
        terminals = {t.name for t in tree.get_terminals()}

        for tip in terminals:
            states[tip] = set([1]) if gene_vector.get(tip, 0) == 1 else set([0])

        tree_clades = list(tree.find_clades(order="postorder"))

        for clade in tree_clades:
            if clade.is_terminal():
                continue
            child_states = [states[child.name] if child.is_terminal() else states[child] for child in clade.clades]
            intersect = set.intersection(*child_states)
            if intersect:
                states[clade] = intersect
            else:
                states[clade] = set.union(*child_states)

        return states

    def pattern_to_dict(gene_row):
        return (gene_row > 0).astype(int).to_dict()

    gene_pairs = list(combinations(presence_absence_matrix.index, 2))
    scores = {}

    for g1, g2 in gene_pairs:
        vec1 = pattern_to_dict(presence_absence_matrix.loc[g1])
        vec2 = pattern_to_dict(presence_absence_matrix.loc[g2])

        states1 = reconstruct_parsimony_states(vec1)
        states2 = reconstruct_parsimony_states(vec2)

        shared_changes = 0
        total_changes = 0

        for clade in tree.get_nonterminals():
            child1, child2 = clade.clades
            for state_dict in [states1, states2]:
                parent_state = list(state_dict[clade])[0]
                for child in [child1, child2]:
                    child_state = list(state_dict[child.name] if child.is_terminal() else state_dict[child])[0]
                    if parent_state != child_state:
                        total_changes += 1

        for clade in tree.get_nonterminals():
            child1, child2 = clade.clades

            parent1 = list(states1[clade])[0]
            parent2 = list(states2[clade])[0]

            for child in [child1, child2]:
                child1_state = list(states1[child.name] if child.is_terminal() else states1[child])[0]
                child2_state = list(states2[child.name] if child.is_terminal() else states2[child])[0]

                change1 = parent1 != child1_state
                change2 = parent2 != child2_state

                if change1 and change2:
                    shared_changes += 1

        score = shared_changes / total_changes if total_changes > 0 else np.nan
        scores[(g1, g2)] = score

    return scores
