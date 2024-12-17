import json
import pandas as pd
import random

def load_csv(filepath):
    """
    Load CSV file with robust error handling
    """
    try:
        return pd.read_csv(filepath, low_memory=False)
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None

def extract_gene_ids(row):
    """
    Extract gene IDs from a row, handling different formats
    """
    if pd.isna(row):
        return set()
    
    # Handle different possible formats
    genes = str(row).split(';')
    genes = [gene.strip().split(':')[0].strip().lower() for gene in genes]
    return set(genes)

def find_filtered_out_candidates(before_filtering_path, after_filtering_path):
    """
    Find candidates that were present before filtering but removed after filtering
    """
    # Load both dataframes
    df_before = load_csv(before_filtering_path)
    df_after = load_csv(after_filtering_path)
    
    if df_before is None or df_after is None:
        print("Could not load one or both files.")
        return []
    
    # Print initial dataset sizes
    print(f"Before filtering dataset size: {len(df_before)}")
    print(f"After filtering dataset size: {len(df_after)}")
    
    # Extract pathway and gene information from both dataframes
    before_pathways = df_before['pathway'].tolist()
    before_genes = df_before['genes'].apply(lambda x: set(str(x).lower().split(','))).tolist()
    before_window_cluster_genes = df_before['window_cluster_genes'].apply(lambda x: set(str(x).split(','))).tolist()
    
    after_pathways = df_after['Pathway (Occurrence)'].tolist()
    after_genes_lists = df_after['Gene IDs'].apply(extract_gene_ids).tolist()
    
    # Find filtered out candidates
    filtered_out = []
    remaining_candidates = []
    
    for i, (pathway, gene_set, window_cluster_genes) in enumerate(zip(before_pathways, before_genes, before_window_cluster_genes)):
        # Flag to check if this candidate is filtered out
        is_filtered_out = True
        
        # Check for substring matches
        for j, (after_pathway, after_genes) in enumerate(zip(after_pathways, after_genes_lists)):
            # Check if ALL after_genes are substrings of genes in gene_set
            substring_match = all(
                any(after_gene in before_gene for before_gene in gene_set) 
                for after_gene in after_genes
            )
            
            if substring_match:
                is_filtered_out = False
                remaining_candidates.append({
                    'bafore_pathway': pathway,
                    'before_genes': list(gene_set),
                    'after_pathway': after_pathway,
                    'after_genes': list(after_genes)
                })
                break
        
        # If no match found, it's truly filtered out
        if is_filtered_out:
            filtered_out.append({
                'pathway': pathway,
                'genes': list(gene_set),
                'window_cluster_genes': list(window_cluster_genes),
                'source_file': df_before.loc[i, 'source_file'] if 'source_file' in df_before.columns else 'Unknown'
            })
    
    # Print detailed filtering information
    print(f"Candidates completely filtered out: {len(filtered_out)}")
    print(f"Candidates with partial/matching representation: {len(remaining_candidates)}")
    
    return filtered_out

def main():
    # Paths to your files - PLEASE VERIFY THESE PATHS
    before_filtering_path = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/unique_clusters_start_end.csv'
    after_filtering_path = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/candidates.csv'
    output_path = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/filtered_out_candidates_random.json'
    
    # Find filtered out candidates
    filtered_out_candidates = find_filtered_out_candidates(before_filtering_path, after_filtering_path)
    
    # If more than 5 filtered out, randomly select 5
    if len(filtered_out_candidates) > 5:
        filtered_out_candidates = random.sample(filtered_out_candidates, 5)
        

    # Write filtered out candidates to output file
    with open(output_path, 'w') as f:
        json.dump(filtered_out_candidates, f, indent=2)
    
    # Print results of completely filtered out candidates
    print("\nCompletely Filtered Out Candidates:")
    for candidate in filtered_out_candidates:
        print("\nPathway:", candidate['pathway'])
        print("Genes:", candidate['genes'])
        print("Window Cluster Genes:", candidate['window_cluster_genes'])
        print("Source File:", candidate['source_file'])

if __name__ == "__main__":
    main()