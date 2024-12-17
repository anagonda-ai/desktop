import csv
from collections import defaultdict

# File paths
candidates_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/unique_clusters_start_end.csv'
gene_tagging_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/combined_updated_genes.csv'
output_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/candidates.csv'

# Data structure to store gene tagging details from the gene tagging file
gene_tagging = {}

# Step 1: Read the gene tagging file and store EC numbers
with open(gene_tagging_file, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        gene_id = row['Query ID']
        ec_number = row['Predicted EC number']
        gene_tagging[gene_id] = ec_number

# Step 2: Prepare to process candidates
filtered_candidates = []

# Step 3: Process each candidate row
with open(candidates_file, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        pathway = row['pathway']  # Keep the occurrence for clarity
        genes = row['window_cluster_genes'].split(',')
        
        families = set()
        subfamilies = set()
        genes_with_ec = set()

        for gene in genes:
            if gene in gene_tagging:  # Only process genes with EC numbers
                ec_number = gene_tagging[gene]
                family = ec_number.split('.')[0]
                subfamily = '.'.join(ec_number.split('.')[:2])
                
                families.add(family)
                subfamilies.add(subfamily)
                genes_with_ec.add((gene, ec_number))

        # Apply filtering criteria
        if len(families) >= 2:
            if len(subfamilies) >= 3:
                filtered_candidates.append((pathway, genes_with_ec))
            elif len(subfamilies) == 2:
                subfamily_families = defaultdict(set)
                for gene_id, ec_number in genes_with_ec:
                    family = ec_number.split('.')[0]
                    subfamily = '.'.join(ec_number.split('.')[:2])
                    subfamily_families[subfamily].add(family)
                if all(len(families) > 1 for families in subfamily_families.values()):
                    filtered_candidates.append((pathway, genes_with_ec))
                    
# Step 4: Write the filtered candidates to a new CSV file
with open(output_file, mode='w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Pathway (Occurrence)', 'Gene IDs and Predicted EC numbers'])
    for pathway, genes_with_ec in filtered_candidates:
        gene_ec_list = [f"{gene_id}:{ec_number}" for gene_id, ec_number in genes_with_ec]
        writer.writerow([pathway, '; '.join(gene_ec_list)])

print(f"Filtered pathways have been written to {output_file}")
