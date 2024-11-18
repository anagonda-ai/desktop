import os
import csv
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

def create_pathways_dict(pathways_file):
    pathways_dict = {}
    with open(pathways_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            unique_id = row['UNIQUE-ID']
            gene_ids = [row[col].lower() for col in row if col.startswith('GENE-ID') and row[col]]
            pathways_dict[unique_id] = gene_ids
    return pathways_dict

def read_csv_files(input_csv_dir):
    file_paths = []
    for file_name in os.listdir(input_csv_dir):
        file_path = os.path.join(input_csv_dir, file_name)
        if os.path.isfile(file_path) and file_name.endswith('.csv'):
            file_paths.append(file_path)
            print(f"Found CSV file: {file_path}")
    return file_paths

def process_file(file_path, pathway_dict, window_size=10, min_genes=3):
    print(f"Processing file: {file_path}")
    potential_groups = set()
    df = pd.read_csv(file_path)
    
    # Sort the DataFrame by 'start' and 'end' columns
    df_sorted = df.sort_values(by=['start', 'end'])
    
    # Use a sliding window approach to find potential groups
    for i in range(len(df_sorted) - window_size + 1):
        window = df_sorted.iloc[i:i + window_size]
        gene_ids = window['id'].tolist()
        
        # Check if at least 3 genes in the window belong to the same pathway
        for pathway, genes in pathway_dict.items():
            cluster_genes = [gene for gene in genes if any(gene in g for g in gene_ids)]
            if len(cluster_genes) >= min_genes:
                group = {
                    'pathway': pathway,
                    'genes': tuple(sorted(gene_ids)),
                    'cluster_genes': tuple(sorted(cluster_genes)),
                    'start': window['start'].min(),
                    'end': window['end'].max()
                }
                potential_groups.add(tuple(sorted(group.items())))

    return potential_groups

def save_potential_groups_to_file(potential_groups, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['pathway', 'genes', 'cluster_genes', 'start', 'end']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for group in potential_groups:
            group_dict = dict(group)
            writer.writerow({
                'pathway': group_dict['pathway'],
                'genes': ','.join(group_dict['genes']),
                'cluster_genes': ','.join(group_dict['cluster_genes']),
                'start': group_dict['start'],
                'end': group_dict['end']
            })

def process_csv_files_concurrently(genome_dirs, pathway_dict, max_workers=4):
    all_potential_groups = set()
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for genome_dir in genome_dirs:
            file_paths = read_csv_files(genome_dir)
            for file_path in file_paths:
                futures.append(executor.submit(process_file, file_path, pathway_dict))
        
        for future in as_completed(futures):
            try:
                potential_groups = future.result()
                all_potential_groups.update(potential_groups)
            except Exception as exc:
                print(f"Error processing file: {exc}")
    
    return all_potential_groups

def main():
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    output_file = "/groups/itay_mayrose/alongonda/desktop/potential_groups.csv"

    pathway_dict = create_pathways_dict(pathways_file)
    
    all_potential_groups = process_csv_files_concurrently(genome_dirs, pathway_dict, max_workers=32)
    
    # Save the potential groups to a file
    save_potential_groups_to_file(all_potential_groups, output_file)
    print(f"Saved potential groups to: {output_file}")

if __name__ == "__main__":
    main()