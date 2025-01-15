import csv
import os

# Define the input directory and output directory
INPUT_DIR = '/groups/itay_mayrose/alongonda/Plant_MGC/csv_files'
OUTPUT_DIR = '/groups/itay_mayrose/alongonda/Plant_MGC/fasta_files'

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_cluster_name(filename):
    """Extract the cluster name from the file name."""
    return os.path.basename(filename).split('.')[0]

def convert_csv_to_fasta(input_csv, output_fasta, cluster_name):
    """Convert a CSV file to a FASTA file."""
    with open(input_csv, 'r') as csv_file, open(output_fasta, 'w') as fasta_file:
        csv_reader = csv.DictReader(csv_file)
        
        for row in csv_reader:
            locus_tag = row['Locus_Tag']
            translation = row['Translation']
            
            if not locus_tag or not translation:
                continue
            
            fasta_file.write(f'>{locus_tag} | {cluster_name}\n')
            fasta_file.write(f'{translation}\n')

def process_files(input_dir, output_dir):
    """Process all CSV files in the input directory and convert them to FASTA files."""
    for filename in os.listdir(input_dir):
        if filename.endswith('.csv'):
            input_csv = os.path.join(input_dir, filename)
            output_fasta = os.path.join(output_dir, filename.replace('.csv', '.fasta'))
            cluster_name = get_cluster_name(input_csv)
            convert_csv_to_fasta(input_csv, output_fasta, cluster_name)

def main():
    process_files(INPUT_DIR, OUTPUT_DIR)
    print(f"FASTA files created in folder: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()