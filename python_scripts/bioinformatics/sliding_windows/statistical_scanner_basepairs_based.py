from threading import Lock
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def create_pathways_dict(pathways_file):
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

def clean_and_process_genes_and_pathways(genes_and_pathways):
    processed = {}
    for gene, values in genes_and_pathways.items():
        stripped_values = [''.join(c for c in v if c not in "[]'\"").split(",") for v in values]
        processed[gene] = stripped_values
    return processed

def find_first_common_element(genes_and_pathways, min_genes):
    from itertools import product, combinations
    processed = clean_and_process_genes_and_pathways(genes_and_pathways)
    
    gene_sets = {}
    for gene, paths in processed.items():
        flattened = set(item.strip() for sublist in paths for item in sublist)
        gene_sets[gene] = flattened
    
    gene_combinations = combinations(gene_sets.keys(), min_genes)
    
    for selected_genes in gene_combinations:
        selected_sets = [gene_sets[gene] for gene in selected_genes]
        all_combinations = product(*selected_sets)
        
        for combination in all_combinations:
            if len(set(combination)) == 1:
                return combination[0], list(selected_genes)
    
    return None, []

def process_df(filtered_df, max_kbp, min_genes, file_path, file_lock, output_file):
    total_matches = 0
    chromosomes = filtered_df['chromosome'].unique()
    
    for chromosome in chromosomes:
        chromosome_data = filtered_df[filtered_df['chromosome'] == chromosome]
        num_genes = len(chromosome_data)
        i = 0
        while i < num_genes:
            genes_and_pathways = dict()
            window = [chromosome_data.iloc[i]]
            genes_and_pathways[chromosome_data.iloc[i]['metabolic_gene']] = [chromosome_data.iloc[i]['pathway']]
            start_pos = chromosome_data.iloc[i]['start']
            end_pos = chromosome_data.iloc[i]['end']
            
            for j in range(i+1, num_genes):
                current_end = chromosome_data.iloc[j]['end']
                if (current_end - start_pos <= max_kbp * 1000) and (current_end > end_pos):
                    window.append(chromosome_data.iloc[j])
                    if (chromosome_data.iloc[j]['metabolic_gene'] not in genes_and_pathways.keys()):
                        genes_and_pathways[chromosome_data.iloc[j]['metabolic_gene']] = [chromosome_data.iloc[j]['pathway']]
                    else:
                        genes_and_pathways[chromosome_data.iloc[j]['metabolic_gene']].append(chromosome_data.iloc[j]['pathway'])
                    end_pos = current_end
                else:
                    break
                    
            if len(genes_and_pathways.keys()) >= min_genes:
                window_df = pd.DataFrame(window)
                pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_genes)
                if pathway:
                    group = {
                        'pathway': pathway,
                        'genes': ','.join(window_df['id']),
                        'metabolic_genes': ','.join(metabolic_genes),
                        'start': window_df['start'].min(),
                        'end': window_df['end'].max(),
                        'source_file': file_path
                    }
                    with file_lock:
                        mode = 'a' if os.path.exists(output_file) else 'w'
                        header = not os.path.exists(output_file)
                        pd.DataFrame([group]).to_csv(output_file, mode=mode, header=header, index=False)
                        total_matches += 1
            i += 1
    return total_matches

def process_file(file_path, output_file, file_lock, max_kbp, min_genes):
    print(f"Processing file: {file_path} with max kbp: {max_kbp}")
    df = pd.read_csv(file_path)
    df["index"] = df.index
    filtered_df = df[df["pathway"].notna()]
    
    total_matches = process_df(filtered_df, max_kbp, min_genes, file_path, file_lock, output_file)
    
    random_advantage = 10
    shuffled_df = df.copy()
    total_random_matches = 0
    
    for i in range(random_advantage):
        permutation = np.random.permutation(len(shuffled_df))
        shuffled_df['pathway'] = shuffled_df['pathway'].values[permutation]
        shuffled_df['metabolic_gene'] = shuffled_df['metabolic_gene'].values[permutation]
        filtered_shuffled_df = shuffled_df[shuffled_df["pathway"].notna()]
        
        total_random_matches += process_df(filtered_shuffled_df, max_kbp, min_genes, file_path, file_lock, output_file) / random_advantage
    
    print(f"Completed file: {file_path}, Matches Found: {total_matches}")
    return total_matches, total_random_matches

def process_genome_dir(genome_dir, output_file, max_workers, max_kbp, min_genes):
    file_paths = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir) 
                  if f.endswith('.csv') and os.path.isfile(os.path.join(genome_dir, f))]
    file_lock = Lock()
    total_matches = 0
    total_random_matches = 0
    
    with tqdm(total=len(file_paths), desc=f"Directory: {os.path.basename(genome_dir)}", unit="file") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    process_file, 
                    file_path, 
                    output_file, 
                    file_lock,
                    max_kbp,
                    min_genes
                )
                for file_path in file_paths
            ]
            for future in as_completed(futures):
                future_matches, future_random_matches = future.result()
                total_matches += future_matches
                total_random_matches += future_random_matches
                pbar.update(1)
                
    print(f"Directory {genome_dir} - Total Matches: {total_matches}")
    print(f"Directory {genome_dir} - Total Random Matches: {total_random_matches}")
    return total_matches, total_random_matches

def create_output_subdir(output_dir, min_genes):
    subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    return subdir

def enrichment_analysis(total_random_matches, total_matches, max_kbp, min_genes, enrichment_output_file):
    if total_random_matches > 0:
        enrichment_ratio = total_matches / total_random_matches
    else:
        enrichment_ratio = float('inf')

    enrichment_stats = {
        'max_kbp': max_kbp,
        'min_genes': min_genes,
        'total_matches': total_matches,
        'total_random_matches': total_random_matches,
        'enrichment_ratio': enrichment_ratio
    }

    enrichment_df = pd.DataFrame([enrichment_stats])
    mode = 'a' if os.path.exists(enrichment_output_file) else 'w'
    header = not os.path.exists(enrichment_output_file)
    enrichment_df.to_csv(enrichment_output_file, mode=mode, header=header, index=False)

def main():
    genome_dirs = ["/groups/itay_mayrose/alongonda/full_genomes/mgc_enriched_files"]
    pathways_file = "/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms/merged_pathways.csv"
    output_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/basepairs_based_10_most_enriched_genomes_statistical_tests"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")

    min_genes = 3
    max_kbp_values = [20, 40, 50, 60, 80, 100]
    
    for max_kbp in max_kbp_values:
        min_genes_subdir = create_output_subdir(output_dir, min_genes)
        output_file = os.path.join(min_genes_subdir, f"potential_groups_w{max_kbp}kbp.csv")
        enrichment_output_file = os.path.join(min_genes_subdir, f"enrichment_w{max_kbp}kbp.csv")
        
        if os.path.exists(output_file):
            os.remove(output_file)

        total_matches = 0
        total_random_matches = 0
        for genome_dir in genome_dirs:
            print(f"Processing genome directory: {genome_dir} with max kbp: {max_kbp} and min_genes: {min_genes}")
            matches, random_matches = process_genome_dir(genome_dir, output_file, max_workers, max_kbp, min_genes)
            total_matches += matches
            total_random_matches += random_matches
        
        print(f"TOTAL MATCHES FOUND for max kbp {max_kbp} and min_genes {min_genes}: {total_matches}")
        print(f"TOTAL RANDOM MATCHES FOUND for max kbp {max_kbp} and min_genes {min_genes}: {total_random_matches}")
        
        enrichment_analysis(total_random_matches, total_matches, max_kbp, min_genes, enrichment_output_file)
        
        print(f"Results saved to: {output_file}")
        print(f"Enrichment statistics saved to: {enrichment_output_file}")

if __name__ == "__main__":
    main()