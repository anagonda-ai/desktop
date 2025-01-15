import csv
from collections import defaultdict

# File paths
candidates_file = '/groups/itay_mayrose/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs_chromosome_sorted/unique_potential_groups_w10.csv'
gene_tagging_file = '/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/results/combined_updated_genes.csv'
output_file = '/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/results/candidates.csv'
small_clusters_file = '/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/results/small_clusters_candidates.csv'

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
        window_genes = row['genes']
        start = int(row['start'])
        end = int(row['end'])
        
        families = set()
        subfamilies = set()
        genes_with_ec = set()

        for gene in genes:
            if gene.lower() in gene_tagging:  # Only process genes with EC numbers
                ec_number = gene_tagging[gene.lower()]
                family = ec_number.split('.')[0]
                subfamily = '.'.join(ec_number.split('.')[:2])
                
                families.add(family)
                subfamilies.add(subfamily)
                genes_with_ec.add((gene, ec_number, start, end))
        
        # Apply filtering criteria
        if len(families) >= 2:
            if len(subfamilies) >= 3:
                filtered_candidates.append((pathway, genes_with_ec))
            elif len(subfamilies) == 2:
                subfamily_families = defaultdict(set)
                for gene_id, ec_number, _, _ in genes_with_ec:
                    family = ec_number.split('.')[0]
                    subfamily = '.'.join(ec_number.split('.')[:2])
                    subfamily_families[subfamily].add(family)
                if all(len(families) > 1 for families in subfamily_families.values()):
                    filtered_candidates.append((pathway, genes_with_ec))

                    
with open(output_file, mode='w', newline='') as outfile, open(small_clusters_file, mode='w', newline='') as small_clusters_outfile:
    writer = csv.writer(outfile)
    writer_small = csv.writer(small_clusters_outfile)
    writer.writerow(['Pathway (Occurrence)', 'Gene IDs and Predicted EC numbers', 'Gene IDs', 'Start', 'End', 'Cluster Size (kbp)'])
    writer_small.writerow(['Pathway (Occurrence)', 'Gene IDs and Predicted EC numbers', 'Gene IDs', 'Start', 'End', 'Cluster Size (kbp)'])
    for pathway, genes_with_ec in filtered_candidates:
        # Calculate the cluster size
        min_start = min(start for _, _, start, _ in genes_with_ec)
        max_end = max(end for _, _, _, end in genes_with_ec)
        cluster_size = (max_end - min_start) / 1000
        # Write the pathway and gene details
        gene_ec_list = [f"{gene_id}:{ec_number}" for gene_id, ec_number, _, _ in genes_with_ec]
        gene_list = [f"{gene_id}" for gene_id, _, _, _ in genes_with_ec]
        writer.writerow([pathway, '; '.join(gene_ec_list), '; '.join(gene_list), min_start, max_end, f"{cluster_size} kbp"])
        if cluster_size < 50:
            writer_small.writerow([pathway, '; '.join(gene_ec_list), '; '.join(gene_list), min_start, max_end, f"{cluster_size} kbp"])

print(f"Filtered pathways have been written to {output_file}")
