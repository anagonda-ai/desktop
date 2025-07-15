
class CladeppCoreProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: cladepp_core.py."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass


from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import pandas as pd
from itertools import combinations
import numpy as np
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, fcluster

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
