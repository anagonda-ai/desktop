import pandas as pd

def count_genes_in_pathways(csv_file):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_file)

    # Split the 'Gene IDs and Predicted EC numbers' column by space
    df['Gene IDs List'] = df['Gene IDs and Predicted EC numbers'].str.split()

    # Count the number of genes for each pathway
    df['Gene Count'] = df['Gene IDs List'].apply(len)

    # Group by pathway and sum the gene counts
    pathway_gene_count = df.groupby('Pathway')['Gene Count'].sum().reset_index()

    return pathway_gene_count

def main():
    csv_file = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates.csv'
    output_file = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/candidate_gene_counts.csv'
    
    pathway_gene_count = count_genes_in_pathways(csv_file)

    # Save the results to a CSV file
    pathway_gene_count.to_csv(output_file, index=False)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()