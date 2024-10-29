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
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = {transcript:[start, end]}
                else:
                    gene_dict[gene_id][transcript] = [start, end]
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
    
    for unique_id, gene_ids in pathways_dict.items():
        for gene_id in gene_ids:
            total_genes += 1
            if gene_id in gene_dict:
                transcripts = gene_dict[gene_id]
                print(f"Pathway: {unique_id}, Gene: {gene_id}, Transcripts: {transcripts}")
                found_genes.append([unique_id, gene_id, transcripts])
            else:
                print(f"Pathway: {unique_id}, Gene: {gene_id} not found in dictionary")
                not_found_genes.append([unique_id, gene_id])
                not_found += 1
                
    found_genes_df = pd.DataFrame(found_genes, columns=['Pathway', 'Gene', 'Transcripts'])
    not_found_genes_df = pd.DataFrame(not_found_genes, columns=['Pathway', 'Gene'])
    
    new_dir = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/pmn_genomes"
    os.makedirs(new_dir, exist_ok=True)
    
    found_genes_df.to_csv(os.path.join(new_dir, "found_genes.csv"), index=False)
    not_found_genes_df.to_csv(os.path.join(new_dir, "not_found_genes.csv"), index=False)
    
    print(f"Total number of genes: {total_genes}")
    print(f"Number of genes not found: {not_found}")

def main():
    gene_files = ["/groups/itay_mayrose_nosnap/alongonda/full_genomes/annotations/ensembl_annotations.csv",
                  "/groups/itay_mayrose_nosnap/alongonda/full_genomes/annotations/plaza_annotations.csv",
                  "/groups/itay_mayrose_nosnap/alongonda/full_genomes/annotations/phytozome_annotations.csv"]
    pathways_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    
    gene_dict = load_csv_to_dict(gene_files)
    pathways_dict = create_pathways_dict(pathways_file)
    find_gene_positions(pathways_dict, gene_dict)

if __name__ == "__main__":
    main()