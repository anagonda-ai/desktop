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

def main():
    mgc_directory = "/groups/itay_mayrose/alongonda/desktop/MGCs/Plant_MGC/csv_files"
    candidate_directory = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_without_e2p2_filtered"
    output_file = os.path.join(candidate_directory, "comparison_results.txt")
    check_identity(mgc_directory, candidate_directory, output_file)

if __name__ == "__main__":
    main()