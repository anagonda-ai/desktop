import pandas as pd

def load_file_into_df(file_path):
    genes = []
    sequences = []
    with open(file_path, 'r') as file:
        gene = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if gene is not None:
                    # Append the previous gene and its sequence
                    genes.append(gene)
                    sequences.append(''.join(sequence))
                gene = line[1:]  # Remove '>'
                sequence = []
            else:
                sequence.append(line)
        # Append the last gene and its sequence
        if gene is not None:
            genes.append(gene)
            sequences.append(''.join(sequence))
    
    # Create DataFrame
    df = pd.DataFrame({'gene': genes, 'sequence': sequences})
    return df

def main():
    duplicates = '/groups/itay_mayrose/alongonda/plantcyc/gene_by_db_seq/PHYTOZOME/unique_ids_seq_test.txt'
    unique_ids_seq = '/groups/itay_mayrose/alongonda/plantcyc/gene_by_db_seq/PHYTOZOME/no_duplicates.fasta'
    
    # Load file into DataFrame
    df = load_file_into_df(duplicates)
    
    # Process DataFrame (example: remove duplicates and keep longest sequence)
    df['sequence_length'] = df['sequence'].apply(len)
    df = df.sort_values('sequence_length', ascending=False).drop_duplicates('gene').drop(columns='sequence_length')
    
    # Write to file
    with open(unique_ids_seq, 'w') as unique_ids_seq_file:
        for _, row in df.iterrows():
            unique_ids_seq_file.write(f">{row['gene']}\n{row['sequence']}\n")

if __name__ == '__main__':
    main()