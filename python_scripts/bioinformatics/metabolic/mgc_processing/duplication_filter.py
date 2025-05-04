import os
import csv
from Bio import SeqIO
from collections import defaultdict
from networkx import Graph
from networkx.algorithms.matching import max_weight_matching

def read_sequences_from_csv(file_path):
    sequences = set()
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            sequence = row['Translation'].replace('*', '')  # Remove '*' characters
            sequences.add(sequence)
    return sequences

def read_sequences_from_fasta(file_path):
    sequences = set()
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq).replace('*', '')  # Remove '*' characters
        sequences.add(sequence)
    return sequences

def build_bipartite_graph(candidate_sequences, mgc_sequences):
    graph = Graph()
    for candidate_seq in candidate_sequences:
        for mgc_seq in mgc_sequences:
            if candidate_seq in mgc_seq or mgc_seq in candidate_seq:
                graph.add_edge(candidate_seq, mgc_seq, weight=1)
    return graph

def check_identity(mgc_directory, candidate_directory, output_file):
    unmatched_candidates = []
    with open(output_file, 'w') as out_file:
        for fasta_filename in os.listdir(candidate_directory):
            if fasta_filename.endswith(".fasta") or fasta_filename.endswith(".fa"):
                fasta_path = os.path.join(candidate_directory, fasta_filename)
                candidate_sequences = read_sequences_from_fasta(fasta_path)
                all_match_found = False
                
                for csv_filename in os.listdir(mgc_directory):
                    if csv_filename.endswith(".csv"):
                        csv_path = os.path.join(mgc_directory, csv_filename)
                        mgc_sequences = read_sequences_from_csv(csv_path)
                        
                        graph = build_bipartite_graph(candidate_sequences, mgc_sequences)
                        matching = max_weight_matching(graph, maxcardinality=True)
                        
                        if len(matching) == len(candidate_sequences):
                            out_file.write(f"The candidate {fasta_path} is matched to {csv_path} cluster.\n")
                            all_match_found = True
                            break
                
                if not all_match_found:
                    out_file.write(f"Not all sequences in {fasta_path} have matches in any MGC CSV file.\n")
                    unmatched_candidates.append(fasta_path)
    return unmatched_candidates

def main():
    mgc_directory = "/groups/itay_mayrose/alongonda/datasets/MIBIG/csv_files"
    candidate_directory = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    output_file = os.path.join(candidate_directory, "comparison_results.txt")
    
    # Check the identity of the candidates
    unmatched_candidates = check_identity(mgc_directory, candidate_directory, output_file)
    
    # Get the MGC filepaths
    mgc_filepaths = os.listdir(mgc_directory)
    # Prepend the full path to MGC filepaths
    mgc_filepaths = [os.path.join(mgc_directory, filename) for filename in mgc_filepaths]
    
    # Merge the unmatched candidates with the MGC filepaths
    merged_list = unmatched_candidates + mgc_filepaths
    
    # Write the merged list to a file
    with open(os.path.join(candidate_directory, "merged_list.txt"), 'w') as merged_file:
        for item in merged_list:
            merged_file.write(f"{item}\n")

if __name__ == "__main__":
    main()