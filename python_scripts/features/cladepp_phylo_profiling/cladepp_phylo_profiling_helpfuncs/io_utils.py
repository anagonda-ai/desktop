import os
import glob
import pandas as pd
from joblib import Parallel, delayed

def load_mapping_if_exists(mapping_file):
    if mapping_file and os.path.exists(mapping_file):
        mapping_df = pd.read_csv(mapping_file)
        mapping_df["Organism"] = mapping_df["Organism"].str.replace(" ", "_").str.strip()
        return mapping_df
    return None

def extract_organism_from_subject(subject_gene_str):
    try:
        return subject_gene_str.split(";")[0].split("=")[1].split(".")[0]
    except:
        return "unknown"

def load_one_blast(path):
    if not os.path.exists(path):
        print(f"Warning: file not found: {path}")
        return None
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df["organism"] = df["subject_gene"].map(extract_organism_from_subject)
    return df

def load_selected_blast_results(comp_df, n_threads=4):
    """
    Load BLAST results from directories specified in comp_df.
    Each row contains Directory path and organism metadata.
    """
    all_blast_results = []
    
    for idx, row in comp_df.iterrows():
        directory = row['Directory']
        
        # Extract organism name from directory path
        # The organism name is typically the second-to-last directory component
        organism_dir_name = directory.split('/')[-2]
        
        # Remove common suffixes to get clean organism name
        organism_name = organism_dir_name.replace('_updated', '').replace('.filtered_updated', '')
        
        # Find all CSV files in the directory (BLAST results)
        csv_files = glob.glob(os.path.join(directory, '*.csv'))
        
        for csv_file in csv_files:
            try:
                # Load the BLAST CSV file
                blast_data = pd.read_csv(csv_file)
                
                # Add organism column to identify which organism this data belongs to
                blast_data['organism'] = organism_name
                
                # Add source file info
                blast_data['source_file'] = os.path.basename(csv_file)
                blast_data['source_directory'] = directory
                
                all_blast_results.append(blast_data)
                
            except Exception as e:
                print(f"Warning: Could not load {csv_file}: {e}")
                continue
    
    if not all_blast_results:
        raise ValueError("No BLAST results could be loaded from the provided directories")
    
    # Combine all BLAST results into one DataFrame
    combined_blast_df = pd.concat(all_blast_results, ignore_index=True)
    
    print(f"Loaded BLAST data: {combined_blast_df.shape[0]} rows from {combined_blast_df['organism'].nunique()} organisms")
    print(f"Organisms found: {sorted(combined_blast_df['organism'].unique())}")
    
    return combined_blast_df
