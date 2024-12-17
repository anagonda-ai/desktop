import os
import json
import subprocess
from Bio import SeqIO
import pandas as pd
import seaborn as sns
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
import math

# Input files
json_file = "/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/random_candidates.json"  # JSON file containing pathway genes
fasta_file = "/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/pmn.fasta"       # Original protein FASTA file
output_dir = "/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/filtered_out_blastp_self_comparison_results"  # Base directory for results

# STEP 1: Load JSON file and extract pathways with 'window_cluster_genes'
def load_pathways(json_file):
    with open(json_file, "r") as file:
        pathways = json.load(file)
        return [{
            "pathway": pathway["Pathway (Occurrence)"],
            "Gene IDs": [gene.lower() for gene in pathway["Gene IDs"]]
        } for pathway in pathways if pathway["Gene IDs"]]

# STEP 2: Extract sequences for the window_cluster_genes
def extract_window_genes(input_fasta, output_fasta, target_genes):
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id.lower() in target_genes:
                SeqIO.write(record, out_handle, "fasta")

# STEP 3: Run BLASTp for self-comparison
def run_self_blastp(input_fasta, output_file):
    db_name = input_fasta + "_db"
    subprocess.run(["makeblastdb", "-in", input_fasta, "-dbtype", "prot", "-out", db_name], check=True)
    subprocess.run([
        "blastp", 
        "-query", input_fasta, 
        "-db", db_name, 
        "-out", output_file, 
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ], check=True)
    
# Function to parse BLAST results
def parse_blast_results(file_path):
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    # Read the file with space delimiter in chunks
    return pd.read_csv(file_path, header=None, names=columns, skiprows=1, engine='python', sep="\t")

def calculate_composite_score(args):
    row, max_bit_score = args
    e_value = row['evalue']
    bit_score = row['bitscore']
    alignment_length = row['length']
    query_length = row['qend'] - row['qstart'] + 1
    percent_identity = row['pident']

    if e_value == 0:
        e_value = 1e-180  # Replace zero to avoid math errors

    # Normalized bit score
    normalized_bit_score = bit_score / max_bit_score
    log_e_value = -math.log10(e_value) if e_value > 0 else 0
    score_per_length = bit_score / alignment_length
    score_per_query_length = bit_score / query_length
    percent_identity_weight = percent_identity / 100
    coverage = alignment_length / query_length

    composite_score = (
        normalized_bit_score * log_e_value * score_per_length * 
        score_per_query_length * percent_identity_weight * coverage
    )

    return composite_score


def process_chunk(chunk, max_bit_score):
    chunk['composite_score'] = chunk.apply(lambda row: calculate_composite_score((row, max_bit_score)), axis=1)
    print(chunk)
    return chunk

def normalize_scores(file_path):
    df = parse_blast_results(file_path)
    
    max_bit_score = df['bitscore'].max()
    df = process_chunk(df, max_bit_score)
    
    min_score = df['composite_score'].min()
    max_score = df['composite_score'].max()

    # Normalize the composite score to range between 0 and 1
    df['normalized_composite_score'] = (
        (df['composite_score'] - min_score) / (max_score - min_score)
    ).fillna(0)  # Handle division by zero if all scores are the same
    
    df = df[['qseqid', 'sseqid', 'normalized_composite_score']]
    
    return df

def parse_blast_to_matrix(blast_output):
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df = pd.read_csv(blast_output, sep="\t", names=columns)
    
    # Normalize the results
    df = normalize_scores(blast_output)

    # Create a symmetric matrix using the normalized composite score
    matrix = df.pivot_table(index="qseqid", columns="sseqid", values="normalized_composite_score", fill_value=0)
    matrix = matrix.add(matrix.T, fill_value=0) / 2  # Ensure symmetry
    return matrix

# STEP 5: Plot heatmap for pairwise similarities
def plot_score_heatmap(matrix, pathway_name, output_path):
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, cmap="viridis", annot=True, fmt=".1f", linewidths=0.5, cbar_kws={'shrink': 0.5})
    plt.title(f"Pairwise BLASTp Self-Comparison Heatmap\n{pathway_name}")
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# Main Workflow
def process_pathways(json_file, fasta_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    pathways = load_pathways(json_file)

    if not pathways:
        print("No valid pathways with 'window_cluster_genes' found.")
        return

    for pathway in pathways:
        pathway_name_clean = pathway["pathway"].replace(" ", "_").replace("(", "").replace(")", "")
        window_genes = pathway["Gene IDs"]
        
        # Create directory for this pathway
        pathway_dir = os.path.join(output_dir, pathway_name_clean)
        os.makedirs(pathway_dir, exist_ok=True)
        
        print(f"Processing pathway: {pathway['pathway']} ({len(window_genes)} genes)...")
        
        # Define file paths
        subset_fasta = os.path.join(pathway_dir, f"{pathway_name_clean}_window_genes.fasta")
        blast_output = os.path.join(pathway_dir, f"{pathway_name_clean}_blast_results.txt")
        matrix_file = os.path.join(pathway_dir, f"{pathway_name_clean}_pairwise_matrix.csv")
        heatmap_file = os.path.join(pathway_dir, f"{pathway_name_clean}_heatmap.png")

        # Steps
        extract_window_genes(fasta_file, subset_fasta, window_genes)
        run_self_blastp(subset_fasta, blast_output)
        matrix = parse_blast_to_matrix(blast_output)
        
        # Save results
        matrix.to_csv(matrix_file)
        plot_score_heatmap(matrix, pathway["pathway"], heatmap_file)
        
        print(f"Results saved in '{pathway_dir}'.")

    print("Processing complete. All results are stored in:", output_dir)

# Run the process
process_pathways(json_file, fasta_file, output_dir)
