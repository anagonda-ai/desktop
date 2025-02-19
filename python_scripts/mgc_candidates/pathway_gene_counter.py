import pandas as pd

def count_genes_in_unique_ids(csv_file):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_file)

    # Print the column names for debugging
    print("Column names in the CSV file:", df.columns.tolist())

    # Filter columns that match the pattern 'GENE-ID'
    gene_id_columns = [col for col in df.columns if col.startswith('GENE-ID')]

    # Ensure the 'GENE-ID' columns are not empty and count the genes
    df['Gene Count'] = df[gene_id_columns].notna().sum(axis=1)
    
    # Add an occurrence number to each UNIQUE-ID
    df['Occurrence'] = df.groupby('UNIQUE-ID').cumcount() + 1

    # Combine UNIQUE-ID and Occurrence to create a new column
    df['Pathway'] = df['UNIQUE-ID'] + ' (Occurrence ' + df['Occurrence'].astype(str) + ')'

    # Keep UNIQUE-ID along with the Gene Count for each occurrence
    unique_id_gene_count = df[['Pathway', 'Gene Count']]

    return unique_id_gene_count

def main():
    csv_file = '/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms/merged_pathways.csv'
    output_file = '/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms/pathway_gene_count.csv'
    
    unique_id_gene_count = count_genes_in_unique_ids(csv_file)

    # Save the results to a CSV file
    unique_id_gene_count.to_csv(output_file, index=False)

    # Print the results in the desired format
    for _, row in unique_id_gene_count.iterrows():
        print(f"{row['Pathway']} - {row['Gene Count']} genes")

if __name__ == "__main__":
    main()