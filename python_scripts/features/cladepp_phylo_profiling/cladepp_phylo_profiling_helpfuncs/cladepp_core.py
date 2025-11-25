import pandas as pd
from itertools import combinations
import numpy as np
import os

def build_profile_matrix(blast_df, anchors, output_dir, min_coverage=0.2, clade_id=None):
    output_path = os.path.join(output_dir, f"matrix_raw_{clade_id}.csv")
    
    if os.path.exists(output_path):
        print(f"⚠️ Matrix for clade {clade_id} already exists at {output_path}, skipping build.")
        output_file = pd.read_csv(output_path, index_col=0)
        return output_file
    
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
    
    filtered.to_csv(output_path)
    return filtered

def get_source_file(comparison_csv):
    comp_dir = os.path.dirname(comparison_csv)
    mgc_name = os.path.basename(comp_dir)
    random_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_csv_files/"
    verified_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_csv_files/"
    mibig_organisms = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/gbk_files/organisms.csv"
    mibig_organisms_df = pd.read_csv(mibig_organisms)
    csv_file = os.path.join(random_dir, mgc_name + ".csv") if "RANDOM_MGC" in mgc_name else os.path.join(verified_dir, mgc_name + ".csv")
    if os.path.exists(csv_file):
        if "BGC" in mgc_name:
            source_file = mibig_organisms_df[mibig_organisms_df['MGC'] == mgc_name]['filtered_path'].iloc[0]
        df = pd.read_csv(csv_file)
        if 'source_file' in df.columns:
            source_file = df['source_file'].iloc[0]
        else:
            source_file = None
    else:
        source_file = None
    return source_file

def get_kegg_mean_std(comparison_csv):
    source_file = get_source_file(comparison_csv)
    source_file_basename = os.path.basename(source_file).replace("_annotated", "")
    source_kegg_mapping = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/dataset_organism_mapping_with_fasta.csv"
    kegg_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes"
    source_kegg_mapping_df = pd.read_csv(source_kegg_mapping)
    if source_file is None:
        return None
    # Ensure source_file_basename is just the basename of 'filtered_path' column (should always be)
    source_kegg_mapping_df['filtered_path_basename'] = source_kegg_mapping_df['filtered_path'].apply(os.path.basename)
    kegg_fasta_file = source_kegg_mapping_df[source_kegg_mapping_df['filtered_path_basename'] == source_file_basename]['kegg_fasta'].iloc[0]
    keg_organism = os.path.basename(kegg_fasta_file).split("_")[0]
    kegg_organism_dir = os.path.join(kegg_dir, keg_organism)
    if os.path.exists(kegg_organism_dir):
        return os.path.join(kegg_organism_dir, "best_hits_fixed", "LPP", "LPP_statistics_summary.csv") 
    return None

def normalize_npp(
    matrix,
    output_dir,
    clade_id,
    self_alignment_scores,
    comparison_csv,
):
    """
    Normalize raw BLAST bit-score matrix in two steps:

    1. **Length bias correction (LPP)**:
       LPP_{i,j} = BS_{i,j} / BS_{i,human}
       Implemented by dividing each gene's row by its self-alignment bit score.

    2. **Phylogenetic distance correction (NPP)**:
       For each genome/organism (column j), transform to Z-scores:
       NPP_{i,j} = (LPP_{i,j} - μ_j) / σ_j
       where μ_j and σ_j are the mean and std of column j over all genes.

    Parameters
    ----------
    matrix : pd.DataFrame
        Raw bit-score matrix (genes x organisms).
    output_dir : str
        Directory to write intermediate and final matrices.
    clade_id : Any
        Identifier used in output filenames.
    self_alignment_scores : dict
        Mapping from gene (row index) to its self-alignment bit score BS_{i,human}.
    comparison_csv : str
        Path to comparison CSV file.
    """

    print(f"Input type: {type(matrix)}")
    print(f"Input shape (genes x organisms): {matrix.shape}")

    # Check if we have enough organisms for correlation / Z-scoring
    if matrix.shape[1] < 2:
        raise ValueError(f"Need at least 2 organisms for correlation, got {matrix.shape[1]}")

    # ----- Step 1: correct for protein length (compute LPP) -----
    # denominator[i] = BS_{i,human}; fall back to NaN if missing
    denominator = matrix.index.map(lambda g: self_alignment_scores.get(g, np.nan))
    denominator = pd.Series(denominator, index=matrix.index)

    with np.errstate(divide='ignore', invalid='ignore'):
        lpp_values = matrix.div(denominator, axis=0)

    lpp = pd.DataFrame(lpp_values, index=matrix.index, columns=matrix.columns)
    lpp = lpp.replace([np.inf, -np.inf], np.nan)

    # Optionally save the intermediate LPP matrix for inspection
    lpp_filled = lpp.fillna(0)
    lpp_path = os.path.join(output_dir, f"matrix_lpp_{clade_id}.csv")
    lpp_filled.to_csv(lpp_path)
    print(f"Saved length-normalized LPP matrix to {lpp_path}")

    # ----- Step 2: correct for phylogenetic distance (column-wise Z-score) -----
    # Compute mean and std for each organism (column) across genes,
    LPP_statistics_summary_csv = get_kegg_mean_std(comparison_csv)
    if LPP_statistics_summary_csv is not None:
        LPP_statistics_summary_df = pd.read_csv(LPP_statistics_summary_csv)
    # For each organism column in the input matrix, use the organism's mean and std from summary df for normalization
    # The summary df must have columns "organism", "mean", "std" (or similar)
    # We'll infer the mapping between matrix columns and summary rows by common organism name

    # Try several possible column name variations for robustness
    organism_col_candidates = ["organism"]  # Add more as needed
    mean_col_candidates = ["mean_lpp"]
    std_col_candidates = ["std_lpp"]

    # Find summary df columns for 'organism', 'mean', 'std'
    def find_column(df, candidates):
        for col in candidates:
            if col in df.columns:
                return col
        raise ValueError(f"Could not find column in DataFrame from options: {candidates}")
    
    organism_col = find_column(LPP_statistics_summary_df, organism_col_candidates)
    mean_col = find_column(LPP_statistics_summary_df, mean_col_candidates)
    std_col = find_column(LPP_statistics_summary_df, std_col_candidates)

    # Set up means and stds for each column in the matrix
    col_means = pd.Series(index=matrix.columns, dtype=float)
    col_stds = pd.Series(index=matrix.columns, dtype=float)
    for col in matrix.columns:
        # Try to find summary row matching this organism/column
        match = None
        for idx, row in LPP_statistics_summary_df.iterrows():
            if str(col) in str(row[organism_col]):
                match = row
                break
        if match is not None:
            col_means[col] = float(match[mean_col])
            col_stds[col] = float(match[std_col])
        else:
            # Fall back to column mean/std from this matrix if not found
            col_means[col] = lpp[col].mean()
            col_stds[col] = lpp[col].std(ddof=0)

    # Avoid division by zero: columns with std == 0 will become 0 after fillna
    col_stds_replaced = col_stds.replace(0, np.nan)

    with np.errstate(divide='ignore', invalid='ignore'):
        npp_values = (lpp.sub(col_means, axis=1)).div(col_stds_replaced, axis=1)

    npp = pd.DataFrame(npp_values, index=matrix.index, columns=matrix.columns)
    npp = npp.replace([np.inf, -np.inf], np.nan).fillna(0)

    print(f"Final NPP type: {type(npp)}")
    print(f"Final NPP shape (genes x organisms): {npp.shape}")

    npp_path = os.path.join(output_dir, f"matrix_npp_{clade_id}.csv")
    npp.to_csv(npp_path)
    print(f"Saved NPP matrix (length + phylogenetic distance normalized) to {npp_path}")

    return npp

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
