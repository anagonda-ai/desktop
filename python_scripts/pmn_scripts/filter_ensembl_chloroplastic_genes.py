import os
from Bio import SeqIO
import csv
from pathlib import Path
import pandas as pd

def extract_gene_id(description):
    # Find the part that starts with 'gene:'
    if 'locus=' in description:
        # Split by space and find the element containing 'gene:'
        for part in description.split():
            if part.startswith('locus='):
                return part.split('locus=')[1]
    return None

def get_genes_from_fasta(fasta_dir):
    """Extract all gene IDs from FASTA files in the directory"""
    gene_ids = set()
    for root, dirs, files in os.walk(fasta_dir):
        for file in files:
            if file.endswith('.fasta'):
                fasta_file = os.path.join(root, file)
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    # Assuming the ID in the FASTA header matches the annotation ID
                    gene_id = record.id.split()[0]
                    if gene_id:
                        gene_ids.add(gene_id.lower())
    return gene_ids

def process_annotation_file(annotation_file, allowed_genes):
    """Filter annotation file to keep only rows with genes present in allowed_genes"""
    # filtered_rows = []
    
    # with open(annotation_file, 'r') as f:
    #     content = f.read()
    #     rows = content.strip().split('\n')
        
    #     for row in rows:
    #         if rows.index(row) == 0:
    #             continue
    #         fields = row.split(';')
    #         # Extract gene ID from the first column
    #         id = fields[0]
    #         gene_id = id.split("ID=")[1]
            
    #         if gene_id.lower() in allowed_genes:
    #             filtered_rows.append(row)
    
    # return filtered_rows
    """
    Filter the annotation file to keep only rows with genes present in allowed_genes.
    :param annotation_file: Path to the annotation file (CSV format with ',' delimiter).
    :param allowed_genes: A set or list of gene IDs to filter.
    :return: A filtered pandas DataFrame.
    """
    # Read the CSV file into a pandas DataFrame, handling quoted fields
    df = pd.read_csv(annotation_file, delimiter=',', quotechar='"', dtype=str)
    
    # Extract the gene ID from the 'ID=' part of the 'id' column
    df['gene_id'] = df['id'].str.extract(r'ID=([^;]+)', expand=False)
    
    # Filter the DataFrame based on allowed_genes (case-insensitive)
    filtered_df = df[df['gene_id'].str.lower().isin(allowed_genes)]
    
    # Drop the temporary 'gene_id' column before returning
    filtered_df = filtered_df.drop(columns=['gene_id'])
    
    return filtered_df

def main():
    # Directories
    fasta_dir = '/groups/itay_mayrose/alongonda/full_genomes/plaza_without_chloroplast'
    annotations_dir = '/groups/itay_mayrose/alongonda/full_genomes/plaza/files_to_filter'
    
    # Get all genes from FASTA files
    print("Reading FASTA files...")
    allowed_genes = get_genes_from_fasta(fasta_dir)
    print(f"Found {len(allowed_genes)} genes in FASTA files")
    
    import concurrent.futures

    def process_and_write_file(annotation_file, allowed_genes):
        print(f"Processing {annotation_file.name}")
        
        # Filter the file
        filtered_df = process_annotation_file(annotation_file, allowed_genes)
        
        # Create output filename
        output_file = annotation_file.with_stem(annotation_file.stem + '_filtered')

        # Write filtered content as CSV using pandas
        # df = pd.DataFrame([row.split(',') for row in filtered_rows])
        # df.columns = ['id', 'start', 'end', 'chromosome']
        filtered_df.to_csv(output_file, index=False, header=True, sep=',')
        # # Write filtered content
        # with open(output_file, 'w') as f:
        #     for row in filtered_rows:
        #         f.write(row + '\n')
        
        print(f"Wrote filtered annotations to {output_file}")

    # Process each annotation file using threads
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_and_write_file, annotation_file, allowed_genes) for annotation_file in Path(annotations_dir).glob('*.csv')]
        for future in concurrent.futures.as_completed(futures):
            future.result()

if __name__ == "__main__":
    main()