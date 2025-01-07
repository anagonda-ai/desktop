import csv
from Bio import SeqIO
import os

def extract_genes_and_transcription(gbk_file, output_csv):
    # Open the output file for writing using the csv module
    with open(output_csv, mode='w', newline='') as output:
        writer = csv.writer(output)
        
        # Write the header row
        writer.writerow(["Feature_Type", "Gene", "Locus_Tag", "Product", "Note", "Gene_Functions", "Translation"])
        
        # Parse the GenBank file
        for record in SeqIO.parse(gbk_file, "genbank"):
            print(f"Processing record {record.id}")
            for feature in record.features:
                # Process gene and CDS features
                if feature.type == "CDS":
                    feature_type = feature.type
                    gene = feature.qualifiers.get("gene", [""])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    product = feature.qualifiers.get("product", [""])[0]
                    note = feature.qualifiers.get("note", [""])[0]
                    gene_functions = feature.qualifiers.get("gene_functions", [""])[0]
                    translation = feature.qualifiers.get("translation", [""])[0] if feature.type == "CDS" else ""
                    
                    # Write the extracted information to the CSV file
                    writer.writerow([feature_type, gene, locus_tag, product, note, gene_functions, translation])

# Process all GenBank files in a directory
def process_directory(gbk_directory):
    for filename in os.listdir(gbk_directory):
        if filename.endswith(".gbk"):
            gbk_file = os.path.join(gbk_directory, filename)
            output_csv_dir = os.path.join(gbk_directory, "csv_files")
            os.makedirs(output_csv_dir, exist_ok=True)
            output_csv = os.path.join(output_csv_dir, os.path.basename(gbk_file).replace(".gbk", ".csv"))
            extract_genes_and_transcription(gbk_file, output_csv)

# Example usage
gbk_directory = "/groups/itay_mayrose_nosnap/alongonda/desktop/MGCs/Plant_MGC"  # Replace with your directory path
process_directory(gbk_directory)
