import sys
import pandas as pd

def convert_blast_output_to_csv(blast_output_file, csv_output_file):
    # Define the column headers
    headers = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    # Read the BLAST output file into a DataFrame
    df = pd.read_csv(blast_output_file, sep='\t', header=None, names=headers)
    
    # Convert relevant columns to float, coercing errors to NaN
    float_columns = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                     'sstart', 'send', 'evalue', 'bitscore']
    df[float_columns] = df[float_columns].apply(pd.to_numeric, errors='coerce')
    
    
    # Save the DataFrame to a CSV file
    df.to_csv(csv_output_file,index=False)

def main():
    blast_output_file = sys.argv[1]
    csv_output_file = sys.argv[2]
    
    convert_blast_output_to_csv(blast_output_file, csv_output_file)
    print(f"Converted BLAST output to CSV: {csv_output_file}")

if __name__ == '__main__':
    main()