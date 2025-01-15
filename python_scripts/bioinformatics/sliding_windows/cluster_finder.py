import json
import pandas as pd

def main():
    # Paths to your files - PLEASE VERIFY THESE PATHS
    unfiltered_candidates_path = '/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/unfiltered_candidates_for_pairwise.csv'
    output_path = '/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/random_candidates.json'
    
    # Load the CSV file
    df = pd.read_csv(unfiltered_candidates_path)
    
    # Choose 5 random rows
    random_candidates = df.sample(n=5).to_dict(orient='records')
    
    # Split "Gene IDs" column to list by "; "
    for candidate in random_candidates:
        candidate['Gene IDs'] = candidate['Gene IDs'].split('; ')
    
    # Write the random candidates to the output file
    with open(output_path, 'w') as f:
        json.dump(random_candidates, f, indent=2)

if __name__ == "__main__":
    main()