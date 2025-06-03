import csv
from Bio import SeqIO
import os

import pandas as pd

def extract_genes_and_transcription(gbk_file, output_csv):
    # Create a CSV writer object
    with open(output_csv, mode='w', newline='') as output:
        cluster_start = 0
        cluster_end = 0
        organism = ""
        writer = csv.writer(output)
        
        # Write the header row
        writer.writerow(["Feature_Type", "Gene", "Locus_Tag", "Protein_ID", "Product", "Note", "Gene_Functions", "Translation", "Start", "End"])
        
        # Parse the GenBank file
        for record in SeqIO.parse(gbk_file, "genbank"):
            print(f"Processing record {record.id}")
            cluster_start = int(record.annotations.get("structured_comment").get("antiSMASH-Data").get("Orig. start"))
            cluster_end = int(record.annotations.get("structured_comment").get("antiSMASH-Data").get("Orig. end"))
            organism = record.annotations.get("organism")
            # Extract the cluster information
            for feature in record.features:
                # Process gene and CDS features
                if feature.type == "CDS":
                    feature_type = feature.type
                    gene = feature.qualifiers.get("gene", [""])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    protein_id = feature.qualifiers.get("protein_id", [""])[0]
                    product = feature.qualifiers.get("product", [""])[0]
                    note = feature.qualifiers.get("note", [""])[0]
                    gene_functions = feature.qualifiers.get("gene_functions", [""])[0]
                    translation = feature.qualifiers.get("translation", [""])[0] if feature.type == "CDS" else ""
                    start = int(feature.location.parts[0].start)
                    end = int(feature.location.parts[-1].end)
                    
                    # Write the extracted information to the CSV file
                    writer.writerow([feature_type, gene, locus_tag, protein_id, product, note, gene_functions, translation, start, end])
        return organism, cluster_start, cluster_end

# Process all GenBank files in a directory
def process_directory(gbk_directory, distance_file, organisms_file):
    # Create a directory for the output CSV files
    distances = dict()
    organisms = dict()
    
    for filename in os.listdir(gbk_directory):
        if filename.endswith(".gbk"):
            gbk_file = os.path.join(gbk_directory, filename)
            output_csv_dir = os.path.join(gbk_directory, "csv_files_test")
            os.makedirs(output_csv_dir, exist_ok=True)
            output_csv = os.path.join(output_csv_dir, os.path.basename(gbk_file).replace(".gbk", ".csv"))
            organism, cluster_start, cluster_end = extract_genes_and_transcription(gbk_file, output_csv)
            cluster = os.path.basename(filename).replace(".gbk", "")
            distances[cluster] = (cluster_start, cluster_end)
            organisms[cluster] = organism
            
    
    distances_df = pd.DataFrame.from_dict(distances, orient='index', columns=['Start', 'End'])
    distances_df["Cluster_Size"] = distances_df["End"] - distances_df["Start"]
    distances_df.sort_index(inplace=True)
    print("num of clusters with size < 50000:", len(distances_df[distances_df["Cluster_Size"] < 50000]))
    print("num of clusters with size >= 50000:", len(distances_df[distances_df["Cluster_Size"] >= 50000]))
    distances_df.to_csv(distance_file, index_label='MGC')
    
    organisms_df = pd.DataFrame.from_dict(organisms, orient='index', columns=['Organism'])
    # Drop rows where Organism value is "."
    organisms_df = organisms_df[organisms_df["Organism"] != "."]
    organisms_df.to_csv(organisms_file, index_label='MGC')

def main():
    gbk_directory = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/gbk_files"  # Replace with your directory path
    distance_file = os.path.join(gbk_directory, "distances.csv")
    organisms_file = os.path.join(gbk_directory, "organisms.csv")
    process_directory(gbk_directory, distance_file, organisms_file)

if __name__ == "__main__":
    main()
