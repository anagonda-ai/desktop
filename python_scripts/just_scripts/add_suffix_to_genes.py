import pandas as pd

def add_unique_suffix_to_ids(csv_file, column_name, output_csv):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    
    # Create a dictionary to keep track of the counters for each gene
    gene_counter = {}
    
    # Iterate through the DataFrame and update the values in the specified column
    for index, row in df.iterrows():
        gene_id = row[column_name]
        if gene_id not in gene_counter:
            gene_counter[gene_id] = 0
        else:
            gene_counter[gene_id] += 1
        df.at[index, column_name] = f"{gene_id}.{gene_counter[gene_id]}"
        print(f"{gene_id}.{gene_counter[gene_id]}")
    
    # Save the modified DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    print(f"Added unique suffix to '{column_name}' and saved to {output_csv}")

def main():
    input_csv = "/groups/itay_mayrose/alongonda/full_genomes/annotations/phytozome_annotations.csv"
    output_csv = "/groups/itay_mayrose/alongonda/full_genomes/annotations/phytozome_annotations_with_suffix.csv"
    column_name_to_count = "id"
    
    # Add unique suffix to each id in the specified column and save to a new CSV file
    add_unique_suffix_to_ids(input_csv, column_name_to_count, output_csv)

if __name__ == "__main__":
    main()