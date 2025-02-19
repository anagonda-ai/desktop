from collections import defaultdict
import os
import re
import csv
from concurrent.futures import ThreadPoolExecutor

# Define the directory to search in
search_dir = '/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms'

# Define the output CSV file
output_file_name = 'pathway_savi_status.csv'
reliability_output_file_name = 'pathway_reliability.csv'

output_file = os.path.join(search_dir, output_file_name)
reliability_output_file = os.path.join(search_dir, reliability_output_file_name)

# Define the SAVI pipeline categories to search for in the COMMENT field
savi_keywords = {
    'UPP': 'Ubiquitous Plant Pathways',
    'NPP': 'Non-PMN Pathway',
    'CAPP': 'Conditionally Accepted PMN Pathway',
    'AIPP': 'Accept-if-Predicted Pathway',
    'CVP': 'Common Viridiplantae Pathway',
    'MCP': 'Manually Checked Pathway'
}

high_confidence_categories = {'UPP', 'MCP'}
low_confidence_categories = {'AIPP', 'CAPP', 'NPP','CVP'}

# Function to check if the COMMENT contains SAVI pipeline keywords
def get_savi_category(comment):
    for category, keyword in savi_keywords.items():
        if keyword in comment:
            return category
    return None

# Function to process a single pathways.dat file
def process_file(file_path):
    results = []
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            
            # Split the content into pathway sections
            pathways = content.split('//')
            
            for pathway in pathways:
                # Regex to extract UNIQUE-ID and COMMENT fields
                unique_id_match = re.findall(r'UNIQUE-ID\s*-\s*(\S+)', pathway)
                comment_match = re.findall(r'COMMENT\s*-\s*(.*?)(?=\n\s*/|\Z)', pathway, re.DOTALL)
                
                if unique_id_match:
                    unique_id = unique_id_match[0]
                    
                    # Check if COMMENT exists and extract relevant information
                    comment = comment_match[0] if comment_match else None
                    savi_category = get_savi_category(comment) if comment else None
                    results.append([unique_id, 'Yes' if savi_category else 'No', savi_category if savi_category else ''])
    
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
    return results

# Function to traverse the directory and process all pathways.dat files concurrently
def traverse_directory(search_dir):
    all_results = []
    files_to_process = []

    # Collect all paths to 'pathways.dat' files
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if file == 'pathways.dat':
                file_path = os.path.join(root, file)
                files_to_process.append(file_path)

    # Use ThreadPoolExecutor to process files concurrently
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_file, file_path) for file_path in files_to_process]

        # Wait for all futures to complete and collect results
        for future in futures:
            all_results.extend(future.result())

    return all_results

# Function to save the results to a CSV file
def save_to_csv(results, output_file):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['UNIQUE-ID', 'Part of SAVI Pipeline?', 'Category'])
        writer.writerows(results)
        
def classify_pathways(results):
    """Classify UNIQUE-ID as reliable or non-reliable based on SAVI pipeline data."""
    pathway_counts = defaultdict(lambda: {'high': 0, 'low': 0})
    
    for unique_id, is_savi, category in results:
        if is_savi == 'No' or category in high_confidence_categories:
            pathway_counts[unique_id]['high'] += 1
        elif category in low_confidence_categories:
            pathway_counts[unique_id]['low'] += 1
    
    classifications = []
    for unique_id, counts in pathway_counts.items():
        reliability_score = counts['high'] / (counts['low'] + counts['high'])
        classification = 'reliable' if reliability_score > 0.75 else 'non-reliable'
        classifications.append([unique_id, counts['high'], counts['low'], reliability_score, classification])
    
    return classifications

def save_reliability_to_csv(reliability_results, reliability_output_file):
    """Save pathway reliability classification to CSV."""
    with open(reliability_output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['UNIQUE-ID', 'High-Confidence Count', 'Low-Confidence Count', 'Reliability Score (out of 1.0)', 'Classification'])
        writer.writerows(reliability_results)

# Main execution
if __name__ == "__main__":
    # Traverse the directory and gather results concurrently
    results = traverse_directory(search_dir)
    # Save the results to a CSV file
    save_to_csv(results, output_file)
    
    # Classify UNIQUE-IDs
    reliability_results = classify_pathways(results)
    save_reliability_to_csv(reliability_results, reliability_output_file)
    print(f"Results saved to {output_file}")
