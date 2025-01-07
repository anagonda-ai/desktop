import requests
import pandas as pd

def fetch_conservation_score(gene_id):
    """
    Fetch evolutionary conservation scores for a gene from Ensembl API.
    """
    base_url = "https://rest.ensembl.org"
    endpoint = "/lookup/id/{}?expand=1".format(gene_id)
    headers = {"Content-Type": "application/json"}
    
    response = requests.get(base_url + endpoint, headers=headers)
    if not response.ok:
        print(f"Failed to fetch data for gene_id: {gene_id}")
        print(f"Response status code: {response.status_code}")
        print(f"Response content: {response.content}")
        return None

    data = response.json()
    print(f"Response data for gene_id {gene_id}: {data}")  # Log the entire response data

    # Assuming the conservation score is part of the response data
    # This part needs to be adjusted based on the actual API response structure
    conservation_score = data.get("conservation_score", None)
    if conservation_score is not None:
        return conservation_score
    else:
        print(f"Conservation score not found for gene_id: {gene_id}")
        return None

def analyze_conservation(input_file, output_file):
    """
    Analyze conservation scores for genes listed in the input file.
    """
    # Read gene IDs from input file
    genes_df = pd.read_csv(input_file)
    if "Gene_ID" not in genes_df.columns:
        print("Input file must contain a 'Gene_ID' column.")
        return

    # Fetch conservation scores for each gene
    conservation_scores = []
    for gene_id in genes_df["Gene_ID"]:
        print(f"Fetching conservation score for gene: {gene_id}")
        score = fetch_conservation_score(gene_id)
        conservation_scores.append(score)

    # Add scores to the dataframe
    genes_df["Conservation_Score"] = conservation_scores

    # Save results to output file
    genes_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# Example usage
input_file = "/groups/itay_mayrose_nosnap/alongonda/desktop/test.csv"  # Replace with your input file
output_file = "/groups/itay_mayrose_nosnap/alongonda/desktop/conservation_results.csv"  # Replace with your desired output file
analyze_conservation(input_file, output_file)