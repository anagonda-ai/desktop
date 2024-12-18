import csv
from collections import defaultdict

# File paths
candidates_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/unique_clusters_start_end.csv'
gene_tagging_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/combined_updated_genes.csv'
output_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/candidates.csv'
unfiltered_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/negative_training_set.csv'

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
unfiltered_candidates = []

# Step 3: Process each candidate row
with open(candidates_file, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        pathway = row['pathway']  # Keep the occurrence for clarity
        genes = row['window_cluster_genes'].split(',')
        window_genes = row['genes']
        
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
        
        filtered = False
        # Apply filtering criteria
        if len(families) >= 2:
            if len(subfamilies) >= 3:
                filtered_candidates.append((pathway, genes_with_ec))
                filtered = True
            elif len(subfamilies) == 2:
                subfamily_families = defaultdict(set)
                for gene_id, ec_number in genes_with_ec:
                    family = ec_number.split('.')[0]
                    subfamily = '.'.join(ec_number.split('.')[:2])
                    subfamily_families[subfamily].add(family)
                if all(len(families) > 1 for families in subfamily_families.values()):
                    filtered_candidates.append((pathway, genes_with_ec))
                    filtered = True
                    
        if not filtered and len(genes_with_ec) < 2:
            unfiltered_candidates.append((pathway, genes_with_ec, len(subfamilies), window_genes))
                    
# Step 4: Write the filtered candidates to a new CSV file
with open(output_file, mode='w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Pathway (Occurrence)', 'Gene IDs and Predicted EC numbers', 'Gene IDs'])
    for pathway, genes_with_ec in filtered_candidates:
        gene_ec_list = [f"{gene_id}:{ec_number}" for gene_id, ec_number in genes_with_ec]
        gene_list = [f"{gene_id}" for gene_id, _ in genes_with_ec]
        writer.writerow([pathway, '; '.join(gene_ec_list), '; '.join(gene_list)])
        
# Step 5: Write the unfiltered candidates to a new CSV file
with open(unfiltered_file, mode='w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Pathway (Occurrence)', 'Gene IDs and Predicted EC numbers', 'Gene IDs', 'Num of Subfamilies', 'Num of Genes', 'Window Genes'])
    for pathway, genes_with_ec, len_subfamilies, window_genes in unfiltered_candidates:
        gene_ec_list = [f"{gene_id}:{ec_number}" for gene_id, ec_number in genes_with_ec]
        gene_list = [f"{gene_id}" for gene_id, _ in genes_with_ec]
        writer.writerow([pathway, '; '.join(gene_ec_list), '; '.join(gene_list), len_subfamilies, len(genes_with_ec), window_genes])

print(f"Filtered pathways have been written to {output_file}")
print(f"Unfiltered pathways have been written to {unfiltered_file}")
