import pandas as pd

def load_and_rename_csv(input_csv, output_csv, new_column_name):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(input_csv)
    
    # Rename the first column
    first_column_name = df.columns[0]
    df.rename(columns={first_column_name: new_column_name}, inplace=True)
    
    # Save the modified DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    print(f"Renamed first column to '{new_column_name}' and saved to {output_csv}")

def main():
    input_csv = "/groups/itay_mayrose/alongonda/full_genomes/annotations/phytozome_annotations.csv"
    output_csv = "/groups/itay_mayrose/alongonda/full_genomes/annotations/phytozome_annotations_fixed.csv"
    new_column_name = "id"
    
    load_and_rename_csv(input_csv, output_csv, new_column_name)

if __name__ == "__main__":
    main()