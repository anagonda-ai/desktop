from Bio import SeqIO
import os
import pandas as pd

def get_adjacent_genes(fasta_file, target_gene_id):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, fasta_file)
    
    # Parse the FASTA file and extract gene information
    gene_list = []

    for record in SeqIO.parse(file_path, "fasta"):
        gene_id = record.id.split(".")[0]  # Gene ID from FASTA header
        gene_list.append(gene_id)

    # Find the index of the target gene
    try:
        target_index = gene_list.index(target_gene_id)
    except ValueError:
        print(f"Gene ID {target_gene_id} not found.")
        return None

    # Get adjacent genes
    adjacent_genes = {}

    if target_index > 0:
        adjacent_genes['upstream'] = gene_list[target_index - 1]
    else:
        adjacent_genes['upstream'] = None  # No upstream gene (first gene in list)

    if target_index < len(gene_list) - 1:
        adjacent_genes['downstream'] = gene_list[target_index + 1]
    else:
        adjacent_genes['downstream'] = None  # No downstream gene (last gene in list)

    return adjacent_genes

# Example usage
fasta_file = "/groups/itay_mayrose/alongonda/desktop/arabidopsis/adjacent_genes_coexpression.py"
target_gene_id = input("Enter target gene: \n")  # Replace with your target gene ID
adjacent_genes = get_adjacent_genes(fasta_file, target_gene_id)

if adjacent_genes:
    upstream = adjacent_genes['upstream']
    downstream = adjacent_genes['downstream']
    # Get the directory where the Co-expression files are located
    co_expression_dir = input("Enter co-expression dir: \n")

    # List all files in the script directory
    files_in_script_dir = os.listdir(co_expression_dir)

    # Filter to get only .txt files that start with "Co-expression"
    txt_files = [file for file in files_in_script_dir if file.startswith("Co-expression") and file.endswith(".txt")]

    # List to store DataFrames
    dfs = []

    # Empty DataFrame to hold the merged result
    merged_df = pd.DataFrame()

    # Read and append each .txt file to the list (assuming tab-delimited files)
    for file in txt_files:
        df = pd.read_csv(os.path.join(co_expression_dir, file), delimiter='\t')
        
        # Find upstream and downstream genes
        target_rows = df[df['Gene_A'] == target_gene_id]
        
        # Collect upstream and downstream weights
        upstream_weights = target_rows[target_rows['Gene_B'] == upstream]
        downstream_weights = target_rows[target_rows['Gene_B'] == downstream]
        
        if not upstream_weights.empty and not downstream_weights.empty:
            break
    print(upstream_weights)
    print(downstream_weights)