import os
import csv
import pandas as pd

def load_csv_to_dict(csv_files):
    gene_dict = {}
    for csv_file in csv_files:
        with open(csv_file, mode='r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                transcript = row['id'].replace("cds-", "").lower()
                gene_id = transcript[:transcript.rfind('.')] if '.' in transcript else transcript
                start = int(row['start'])
                end = int(row['end'])
                gene_index = int(row['gene_index'])
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = {transcript: [start, end, gene_index]}
                else:
                    gene_dict[gene_id][transcript] = [start, end, gene_index]
    return gene_dict

def create_pathways_dict(pathways_file):
    pathways_dict = {}
    with open(pathways_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            unique_id = row['UNIQUE-ID']
            gene_ids = [row[col].lower() for col in row if col.startswith('GENE-ID') and row[col]]
            pathways_dict[unique_id] = gene_ids
    return pathways_dict

def find_gene_positions(pathways_dict, gene_dict):
    not_found = 0
    total_genes = 0
    
    found_genes = []
    not_found_genes = []
    
    # Initialize dictionary to store min and max distances
    pathway_distances = {}
    pathway_neighboring_distances = {}

    for unique_id, gene_ids in pathways_dict.items():
        gene_positions = []
        
        for gene_id in gene_ids:
            total_genes += 1
            if gene_id in gene_dict:
                transcripts = gene_dict[gene_id]
                print(f"Pathway: {unique_id}, Gene: {gene_id}, Transcripts: {transcripts}")
                found_genes.append([unique_id, gene_id, transcripts])
                
                # Get the longest transcript (the difference between start and end position)
                longest_transcript = max(transcripts.values(), key=lambda t: int(t[1]) - int(t[0]))
                gene_positions.append((gene_id, longest_transcript[0], longest_transcript[1], longest_transcript[2]))
            else:
                print(f"Pathway: {unique_id}, Gene: {gene_id} not found in dictionary")
                not_found_genes.append([unique_id, gene_id])
                not_found += 1
        
        # Sort genes by the start position of their longest transcript
        gene_positions.sort(key=lambda x: x[1])
        
        # Calculate distances between consecutive genes
        distances = [gene_positions[i+1][1] - gene_positions[i][2] for i in range(len(gene_positions) - 1)]
        neighbor_relationships = [gene_positions[i+1][3] - gene_positions[i][3] for i in range(len(gene_positions) - 1)]
        
        if distances:
            min_distance = min(distances)
            max_distance = max(distances)
            pathway_distances[unique_id] = {'min_distance': min_distance, 'max_distance': max_distance}
        else:
            pathway_distances[unique_id] = {'min_distance': None, 'max_distance': None}
            
        if neighbor_relationships:
            min_distance = min(neighbor_relationships)
            max_distance = max(neighbor_relationships)
            pathway_neighboring_distances[unique_id] = {'min_distance': min_distance, 'max_distance': max_distance}
        else:
            pathway_neighboring_distances[unique_id] = {'min_distance': None, 'max_distance': None}


    found_genes_df = pd.DataFrame(found_genes, columns=['Pathway', 'Gene', 'Transcripts'])
    not_found_genes_df = pd.DataFrame(not_found_genes, columns=['Pathway', 'Gene'])
    
    new_dir = "/groups/itay_mayrose/alongonda/full_genomes/annotations/test_dir/pmn_genomes"
    os.makedirs(new_dir, exist_ok=True)
    
    # Write the distances results to a new file
    with open(os.path.join(new_dir, 'pathway_distances.csv'), 'w') as f:
        f.write('Pathway,Min_Distance,Max_Distance\n')
        for pathway, distances in pathway_distances.items():
            if distances['min_distance'] is not None and distances['max_distance'] is not None:
                f.write(f"{pathway},{distances['min_distance']},{distances['max_distance']}\n")
                
    # Write the neighboring results to a new file
    with open(os.path.join(new_dir, 'pathway_neighboring_distances.csv'), 'w') as f:
        f.write('Pathway,Min_Neighboring_Distance,Max_Neighboring_Distance\n')
        for pathway, distances in pathway_neighboring_distances.items():
            if distances['min_distance'] is not None and distances['max_distance'] is not None:
                f.write(f"{pathway},{distances['min_distance']},{distances['max_distance']}\n")
    
    found_genes_df.to_csv(os.path.join(new_dir, "found_genes.csv"), index=False)
    not_found_genes_df.to_csv(os.path.join(new_dir, "not_found_genes.csv"), index=False)
    
    print(f"Total number of genes: {total_genes}")
    print(f"Number of genes not found: {not_found}")

def main():
    gene_files = ["/groups/itay_mayrose/alongonda/full_genomes/annotations/test_dir/merged_annotations.csv"]
    pathways_file = "/groups/itay_mayrose/alongonda/plantcyc/all_organisms/merged_pathways.csv"
    
    gene_dict = load_csv_to_dict(gene_files)
    pathways_dict = create_pathways_dict(pathways_file)
    find_gene_positions(pathways_dict, gene_dict)

if __name__ == "__main__":
    main()