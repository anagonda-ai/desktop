from threading import Lock
import os
from typing import List
import pandas as pd
import hashlib
import random
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import ahocorasick
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy import stats

def generate_group_hash(group):
    """Generate a unique hash for a group of genes."""
    sorted_genes = sorted(group['pathway_cluster_genes'].split(','))
    hash_string = f"{'|'.join(sorted_genes)}"
    return hashlib.md5(hash_string.encode()).hexdigest()

class WindowStatistics:
    """Track statistics for window configurations."""
    def __init__(self):
        self._lock = Lock()
        self.total_windows = 0
        self.matching_windows = 0
        self.random_matches = 0
        self.random_total = 0
    
    def update_windows(self, total, matching):
        with self._lock:
            self.total_windows += total
            self.matching_windows += matching
    
    def update_random(self, total, matching):
        with self._lock:
            self.random_total += total
            self.random_matches += matching
    
    def get_statistics(self):
        actual_ratio = self.matching_windows / self.total_windows if self.total_windows > 0 else 0
        random_ratio = self.random_matches / self.random_total if self.random_total > 0 else 0
        
        return {
            'actual_ratio': actual_ratio,
            'random_ratio': random_ratio,
            'enrichment_ratio': actual_ratio / random_ratio if random_ratio > 0 else 0,
            'total_windows': self.total_windows,
            'matching_windows': self.matching_windows,
            'random_total': self.random_total,
            'random_matches': self.random_matches
        }

class UniqueMatchTracker:
    """Track unique matches across windows."""
    def __init__(self):
        self._seen_matches = []
        self._seen_hashes = set()
        self._lock = Lock()

    def is_unique_match(self, pathway_cluster):
        new_cluster_genes = set(pathway_cluster['pathway_cluster_genes'].split(','))
        group_hash = generate_group_hash(pathway_cluster)

        with self._lock:
            if group_hash in self._seen_hashes:
                return False

            for match in self._seen_matches:
                existing_genes = set(match.split(','))
                if new_cluster_genes.issubset(existing_genes):
                    return False

            self._seen_matches.append(pathway_cluster['pathway_cluster_genes'])
            self._seen_hashes.add(group_hash)
            return True

def build_aho_corasick(pathway_dict):
    """Build Aho-Corasick automaton for efficient string matching."""
    automaton = ahocorasick.Automaton()
    for pathway, gene_lists in pathway_dict.items():
        for occurrence_idx, gene_list in enumerate(gene_lists):
            unique_occurrence_id = f"{pathway}__{occurrence_idx}"
            for gene in gene_list:
                automaton.add_word(gene.lower(), (unique_occurrence_id, gene.lower()))
    automaton.make_automaton()
    return automaton

def create_pathways_dict(pathways_file):
    """Create dictionary of pathways and their associated genes."""
    pathways_dict = defaultdict(list)
    chunk_size = 10000
    for chunk in pd.read_csv(pathways_file, chunksize=chunk_size, low_memory=False):
        gene_id_cols = [col for col in chunk.columns if col.startswith('GENE-ID.')]
        
        for _, row in chunk.iterrows():
            unique_id = row['UNIQUE-ID']
            gene_ids = [str(val).lower() for col in gene_id_cols 
                       if pd.notna(val := row[col]) and val]
            if gene_ids:
                pathways_dict[unique_id].append(gene_ids)
    return dict(pathways_dict)

def process_window_with_aho_corasick(window_data, automaton, unique_tracker, output_file, file_lock, file_path, min_genes):
    """Process a single window using Aho-Corasick algorithm."""
    gene_ids = [gene.lower() for gene in window_data['id']]
    matches_found = 0
    matched_pathways = defaultdict(lambda: defaultdict(set))
    gene_positions = defaultdict(dict)

    for end_pos, (unique_occurrence_id, gene) in automaton.iter(" ".join(gene_ids)):
        pathway, occurrence_idx = unique_occurrence_id.split('__')
        matched_pathways[pathway][occurrence_idx].add(gene)

    for pathway, occurrences in matched_pathways.items():
        for occurrence_idx, genes in occurrences.items():
            if len(genes) >= min_genes:
                matched_window_genes = []
                for window_gene, start, end in zip(window_data['id'], window_data['start'], window_data['end']):
                    if any(gene in window_gene.lower() for gene in genes):
                        matched_window_genes.append(window_gene)
                        gene_positions[window_gene] = {'start': start, 'end': end}

                unique_genes = set(genes)
                matched_positions = [
                    f"{gene}({gene_positions[gene]['start']}-{gene_positions[gene]['end']})"
                    for gene in matched_window_genes if gene in gene_positions
                ]
                
                group = {
                    'pathway': f"{pathway} (Occurrence {occurrence_idx})",
                    'genes': ','.join(window_data['id']),
                    'pathway_cluster_genes': ','.join(sorted(unique_genes)),
                    'window_cluster_genes': ','.join(sorted(matched_window_genes)),
                    'gene_positions': ','.join(matched_positions),
                    'start': window_data['start'].min(),
                    'end': window_data['end'].max(),
                    'distance': f"{(window_data['end'].max() - window_data['start'].min())/1000}kbp",
                    'source_file': file_path
                }
                
                pathway_cluster = {'pathway_cluster_genes': group['pathway_cluster_genes']}
                if unique_tracker.is_unique_match(pathway_cluster):
                    matches_found += 1
                    if output_file:
                        with file_lock:
                            mode = 'a' if os.path.exists(output_file) else 'w'
                            header = not os.path.exists(output_file)
                            pd.DataFrame([group]).to_csv(output_file, mode=mode, header=header, index=False)
    
    return matches_found

def process_window_with_statistics(window_data, automaton, unique_tracker, output_file, file_lock, file_path, min_genes, stats_tracker, is_random=False):
    """Process a window and update statistics."""
    matches = process_window_with_aho_corasick(
        window_data, automaton, unique_tracker, 
        output_file if not is_random else None,
        file_lock, file_path, min_genes
    )
    
    if is_random:
        stats_tracker.update_random(1, 1 if matches > 0 else 0)
    else:
        stats_tracker.update_windows(1, 1 if matches > 0 else 0)
    
    return matches

def process_file_with_statistics(file_path, pathway_tries, output_file, unique_tracker, file_lock, window_size, min_genes, stats_tracker):
    """Process a single file and collect statistics."""
    print(f"Processing file: {file_path} with window size: {window_size}")
    total_matches = 0
    df = pd.read_csv(file_path, usecols=['id', 'start', 'end'])
    df = df.sort_values(by=['start', 'end'])
    
    # Process actual windows
    num_windows = len(df) - window_size + 1
    with tqdm(total=num_windows, desc=f"File: {os.path.basename(file_path)}", unit="window") as pbar:
        for i in range(num_windows):
            window = df.iloc[i:i + window_size]
            window_matches = process_window_with_statistics(
                window, pathway_tries, unique_tracker, output_file,
                file_lock, file_path, min_genes, stats_tracker, is_random=False
            )
            total_matches += window_matches
            pbar.update(1)
    
    # Process random windows
    all_genes = df['id'].tolist()
    with tqdm(total=num_windows, desc=f"Random trials: {os.path.basename(file_path)}", unit="trial") as pbar:
        for _ in range(num_windows):
            random_genes = random.sample(all_genes, window_size)
            random_window = pd.DataFrame({
                'id': random_genes,
                'start': 1,
                'end': 1
            })
            process_window_with_statistics(
                random_window, pathway_tries, unique_tracker, output_file,
                file_lock, file_path, min_genes, stats_tracker, is_random=True
            )
            pbar.update(1)
    
    print(f"Completed file: {file_path}, Matches Found: {total_matches}")
    return total_matches

def process_genome_dirs(genome_dirs: List[str],
                       pathway_tries: ahocorasick.Automaton,
                       output_file: str,
                       max_workers: int,
                       window_size: int,
                       min_genes: int,
                       stats_tracker: WindowStatistics) -> int:
    """Process all genome directories in parallel."""
    # Collect all files from all directories
    all_files = []
    for genome_dir in genome_dirs:
        files = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir)
                if f.endswith('.csv') and os.path.isfile(os.path.join(genome_dir, f))]
        all_files.extend(files)
    
    unique_tracker = UniqueMatchTracker()
    file_lock = Lock()
    total_matches = 0
    
    print(f"Processing {len(all_files)} total files across all directories")
    
    with tqdm(total=len(all_files), desc="Processing files", unit="file") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    process_file_with_statistics,
                    file_path, pathway_tries, output_file,
                    unique_tracker, file_lock, window_size,
                    min_genes, stats_tracker
                )
                for file_path in all_files
            ]
            
            for future in as_completed(futures):
                total_matches += future.result()
                pbar.update(1)
    
    return total_matches

def create_output_subdir(output_dir, min_genes):
    """Create a subdirectory for each min_genes value."""
    subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    return subdir

def enrichment_ratio_of_clusters_for_windowsize_clustersize(root_directory):
    file_pattern = "statistics"

    # Recursively find all CSV files matching the pattern
    csv_files = [
        os.path.join(dirpath, file)
        for dirpath, _, files in os.walk(root_directory)
        for file in files
        if file.startswith(file_pattern) and file.endswith(".csv")
    ]

    if not csv_files:
        print(f"No files matching pattern '{file_pattern}*.csv' found in {root_directory}")
        return

    # Consolidate data from all files
    rankings = []
    data = []
    for file in csv_files:
        df = pd.read_csv(file)
        if "enrichment_ratio" in df.columns:
            file_name = os.path.basename(file)
            median_enrichment = df["enrichment_ratio"].median()
            rankings.append({"file_name": file_name, "median_enrichment": median_enrichment})
            for enrichment_ratio in df["enrichment_ratio"]:
                data.append({"file_name": file_name, "enrichment_ratio": enrichment_ratio})

    # Create a DataFrame for ranking
    ranking_df = pd.DataFrame(rankings)
    ranking_df = ranking_df.sort_values(by="median_enrichment", ascending=False)

    # Save the ranking to a CSV file
    output_ranking_path = os.path.join(root_directory, "file_rankings.csv")
    ranking_df.to_csv(output_ranking_path, index=False)
    print(f"Ranking saved to {output_ranking_path}")

    # Create a consolidated DataFrame for plotting
    consolidated_df = pd.DataFrame(data)

    # Plot the enrichment_ratio vs file_name (ordered by ranking)
    ordered_files = ranking_df["file_name"]
    plt.figure(figsize=(12, 6))
    plt.boxplot(
        [
            consolidated_df[consolidated_df["file_name"] == name]["enrichment_ratio"].values
            for name in ordered_files
        ],
        labels=ordered_files,
        showfliers=False,
    )
    plt.xticks(rotation=90)
    plt.title("Enrichment Ratio by File Name (Ordered by Median)")
    plt.xlabel("File Name")
    plt.ylabel("Enrichment Ratio")
    plt.tight_layout()

    # Save the plot
    output_plot_path = output_ranking_path = os.path.join(root_directory, "enrichment_ratio_plot.png")
    plt.savefig(output_plot_path)
    print(f"Plot saved to {output_plot_path}")
    plt.show()


def main():
    """Main execution function with improved parallel processing."""
    # Input directories and files
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations_with_chromosomes"
    ]
    pathways_file = "/groups/itay_mayrose_nosnap/alongonda/plantcyc/all_organisms/merged_pathways.csv"
    output_dir = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/sliding_window_outputs_with_statistics"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Use maximum available workers
    max_workers = os.cpu_count()  # Use all CPU cores
    print(f"Using {max_workers} workers")
    
    print("Loading pathways...")
    pathway_dict = create_pathways_dict(pathways_file)
    pathway_tries = build_aho_corasick(pathway_dict)
    print(f"Loaded {len(pathway_dict)} pathways with Aho-Corasick automaton")
    
    results = []
    # window_sizes = range(5, 21)
    window_sizes = [10]
    total_configs = sum(((size // 2) + 1) - 2 for size in window_sizes)
    
    with tqdm(total=total_configs, desc="Processing configurations", unit="config") as config_pbar:
        for window_size in window_sizes:
            max_min_genes = (window_size // 2) + 1
            
            for min_genes in range(3, max_min_genes + 1):
                print(f"\nProcessing window_size={window_size}, min_genes={min_genes}")
                
                min_genes_subdir = create_output_subdir(output_dir, min_genes)
                output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}.csv")
                
                if os.path.exists(output_file):
                    os.remove(output_file)
                
                stats_tracker = WindowStatistics()
                
                # Process all directories at once
                total_matches = process_genome_dirs(
                    genome_dirs, pathway_tries, output_file,
                    max_workers, window_size, min_genes, stats_tracker
                )
                
                # Collect and save statistics
                stats = stats_tracker.get_statistics()
                stats.update({
                    'window_size': window_size,
                    'min_genes': min_genes,
                    'total_matches': total_matches
                })
                
                results.append(stats)
                
                # Save current configuration results
                stats_output_file = os.path.join(
                    min_genes_subdir,
                    f"statistics_w{window_size}_m{min_genes}.csv"
                )
                pd.DataFrame([stats]).to_csv(stats_output_file, index=False)
                
                print(f"\nConfiguration results:")
                print(f"Window size: {window_size}, Min genes: {min_genes}")
                print(f"Actual ratio: {stats['actual_ratio']:.4f}")
                print(f"Random ratio: {stats['random_ratio']:.4f}")
                print(f"Enrichment ratio: {stats['enrichment_ratio']:.4f}")
                print(f"Total matches: {total_matches}")
                print(f"Results saved to: {output_file}")
                
                config_pbar.update(1)
    
    # Save final results
    results_df = pd.DataFrame(results)
    stats_output = os.path.join(output_dir, 'window_configuration_statistics.csv')
    results_df.to_csv(stats_output, index=False)
    
    # Print best configurations
    best_configs = results_df.nlargest(5, 'enrichment_ratio')
    print("\nTop 5 configurations by enrichment ratio:")
    print(best_configs[['window_size', 'min_genes', 'enrichment_ratio', 'actual_ratio', 'random_ratio']])
    
    # Generate final visualization
    enrichment_ratio_of_clusters_for_windowsize_clustersize(output_dir)

if __name__ == "__main__":
    main()