import csv
import os

# Define the input directory and output directory
input_dir = '/groups/itay_mayrose/alongonda/Plant_MGC/csv_files'
output_dir = '/groups/itay_mayrose/alongonda/Plant_MGC/fasta_files'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate through each file in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.csv'):
        input_csv = os.path.join(input_dir, filename)
        output_fasta = os.path.join(output_dir, filename.replace('.csv', '.fasta'))
        
        # Extract the cluster name from the file name
        cluster_name = os.path.basename(input_csv).split('.')[0]
        
        # Open the input CSV file and the output FASTA file
        with open(input_csv, 'r') as csv_file, open(output_fasta, 'w') as fasta_file:
            csv_reader = csv.DictReader(csv_file)
            
            # Iterate through each row in the CSV file
            for row in csv_reader:
                locus_tag = row['Locus_Tag']
                translation = row['Translation']
                
                # Skip rows without a locus tag or translation
                if not locus_tag or not translation:
                    continue
                
                # Write the FASTA formatted entry with the cluster name in the header
                fasta_file.write(f'>{locus_tag} | {cluster_name}\n')
                fasta_file.write(f'{translation}\n')