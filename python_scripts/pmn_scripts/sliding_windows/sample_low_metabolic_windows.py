import os
import pandas as pd
import ahocorasick
from tqdm import tqdm
from collections import defaultdict
from random import sample
from concurrent.futures import ThreadPoolExecutor, as_completed

# Aho-Corasick Implementation to match multiple patterns in text
def build_aho_corasick(pathway_dict):
    automaton = ahocorasick.Automaton()
    for pathway, gene_lists in pathway_dict.items():
        for occurrence_idx, gene_list in enumerate(gene_lists):
            unique_occurrence_id = f"{pathway}__{occurrence_idx}"
            for gene in gene_list:
                automaton.add_word(gene.lower(), (unique_occurrence_id, gene.lower()))
    automaton.make_automaton()
    return automaton

# Create pathways dictionary from pathways file
def create_pathways_dict(pathways_file):
    pathways_dict = defaultdict(list)
    chunk_size = 10000
    for chunk in pd.read_csv(pathways_file, chunksize=chunk_size, low_memory=False):
        gene_id_cols = [col for col in chunk.columns if col.startswith('GENE-ID.')]
        for _, row in chunk.iterrows():
            unique_id = row['UNIQUE-ID']
            gene_ids = [str(val).lower() for col in gene_id_cols if pd.notna(val := row[col]) and val]
            if gene_ids:
                pathways_dict[unique_id].append(gene_ids)
    return dict(pathways_dict)

# Process a sliding window of size `window_size` and sample windows with fewer than 2 metabolic genes
def process_file(file_path, automaton, sampled_windows, sample_size, window_size):
    df = pd.read_csv(file_path, usecols=['id', 'start', 'end'])
    df = df.sort_values(by=['start', 'end'])
    num_windows = len(df) - window_size + 1
    local_sample = []

    with tqdm(total=num_windows, desc=f"Processing {os.path.basename(file_path)}", unit="window") as pbar:
        for i in range(num_windows):
            window = df.iloc[i:i + window_size]
            gene_ids = [gene.lower() for gene in window['id']]
            
            # Use Aho-Corasick to find metabolic genes in the window
            metabolic_genes = set()
            for _, (_, gene) in automaton.iter(" ".join(gene_ids)):
                metabolic_genes.add(gene)
            
            if len(metabolic_genes) < 2:
                local_sample.append({
                    'window_genes': ','.join(window['id']),
                    'start': window['start'].min(),
                    'end': window['end'].max(),
                    'metabolic_gene_count': len(metabolic_genes),
                    'source_file': file_path
                })

            # Limit local samples to prevent excessive memory usage
            if len(local_sample) > sample_size * 2:
                sampled_windows.extend(sample(local_sample, min(sample_size, len(local_sample))))
                local_sample = []

            pbar.update(1)

    # Add any remaining samples from this file
    sampled_windows.extend(sample(local_sample, min(sample_size - len(sampled_windows), len(local_sample))))

# Process all files in a genome directory and save sampled windows
def process_genome_dir(genome_dir, automaton, output_file, sample_size, window_size, max_workers):
    file_paths = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir) if f.endswith('.csv')]
    sampled_windows = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_file, file_path, automaton, sampled_windows, sample_size, window_size)
            for file_path in file_paths
        ]
        for future in as_completed(futures):
            future.result()
    
    # Limit to the final desired sample size
    final_sample = sample(sampled_windows, min(sample_size, len(sampled_windows)))

    # Save results
    pd.DataFrame(final_sample).to_csv(output_file, index=False)
    print(f"Sampled {len(final_sample)} windows saved to: {output_file}")

def main():
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations"
    ]
    pathways_file = "/groups/itay_mayrose_nosnap/alongonda/plantcyc/all_organisms/merged_pathways.csv"
    output_dir = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/sliding_window_outputs"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    max_workers = min(32, os.cpu_count())
    window_size = 10
    sample_size = 5000
    sampled_output_file = os.path.join(output_dir, "sampled_low_metabolic_windows.csv")
    
    print("Loading pathways...")
    pathway_dict = create_pathways_dict(pathways_file)
    automaton = build_aho_corasick(pathway_dict)
    print(f"Loaded {len(pathway_dict)} pathways.")
    
    print("Processing genome directories...")
    for genome_dir in genome_dirs:
        process_genome_dir(genome_dir, automaton, sampled_output_file, sample_size, window_size, max_workers)

if __name__ == "__main__":
    main()
