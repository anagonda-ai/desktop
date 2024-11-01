import os
import csv
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    print(f"Reading CSV files from: {input_csv_dir}")
    file_paths = []
    for file_name in os.listdir(input_csv_dir):
        file_path = os.path.join(input_csv_dir, file_name)
        if os.path.isfile(file_path) and file_name.endswith('.csv'):
            file_paths.append(file_path)
    return file_paths

def process_file(file_path, pathway_dict, max_distance=10, min_genes=3):
    potential_groups = set()
    df = pd.read_csv(file_path)
    
    # Sort the DataFrame by 'start' and 'end' columns
    df_sorted = df.sort_values(by=['start', 'end'])
    
    # Use a more efficient approach to find potential groups
    for pathway, genes in pathway_dict.items():
        gene_set = set(genes)
        indices = df_sorted.index[df_sorted['id'].isin(gene_set)].tolist()
        
        for i in range(len(indices)):
            for j in range(i + min_genes - 1, min(i + max_distance + min_genes, len(indices))):
                window_indices = indices[i:j + 1]
                window = df_sorted.loc[window_indices]
                gene_ids = window['id'].tolist()
                
                if len(gene_ids) >= min_genes:
                    group = {
                        'pathway': pathway,
                        'genes': tuple(sorted(gene_ids)),
                        'start': window['start'].min(),
                        'end': window['end'].max()
                    }
                    potential_groups.add(tuple(sorted(group.items())))
                    print(f"Found potential group: {group}")

    return potential_groups

def save_potential_groups_to_file(potential_groups, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['pathway', 'genes', 'start', 'end']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for group in potential_groups:
            group_dict = dict(group)
            writer.writerow({
                'pathway': group_dict['pathway'],
                'genes': ','.join(group_dict['genes']),
                'start': group_dict['start'],
                'end': group_dict['end']
            })

def main():
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations_sorted"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    output_file = "/groups/itay_mayrose/alongonda/desktop/potential_groups.csv"

    pathway_dict = create_pathways_dict(pathways_file)
    
    all_potential_groups = set()
    
    with ProcessPoolExecutor() as executor:
        futures = []
        for genome_dir in genome_dirs:
            file_paths = read_csv_files(genome_dir)
            for file_path in file_paths:
                futures.append(executor.submit(process_file, file_path, pathway_dict))
        
        for future in as_completed(futures):
            all_potential_groups.update(future.result())
    
    # Save the potential groups to a file
    save_potential_groups_to_file(all_potential_groups, output_file)
    print(f"Saved potential groups to: {output_file}")

if __name__ == "__main__":
    main()