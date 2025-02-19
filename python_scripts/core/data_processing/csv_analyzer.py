import pandas as pd

def extract_and_print_sorted_distances(file_path):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)
    
    # Filter out negative values in the Min_Distance column
    min_distances = df['Min_Neighboring_Distance'][df['Min_Neighboring_Distance'] > 0]
    
    # Extract, remove duplicates, and sort the Min_Distance column
    sorted_min_distances = sorted(set(min_distances.tolist()))
    
    # Extract, remove duplicates, and sort the Max_Distance column
    sorted_max_distances = sorted(set(df['Max_Neighboring_Distance'].tolist()))
    
    # Print the sorted values
    print("Sorted Min_Distance (ignoring negative values):")
    print(sorted_min_distances[0])
    
    print("Sorted Max_Distance:")
    print(sorted_max_distances[-1])

def main():
    file_path = '/groups/itay_mayrose/alongonda/datasets/full_genomes/pmn_genomes/pathway_neighboring_distances.csv'
    extract_and_print_sorted_distances(file_path)

if __name__ == "__main__":
    main()