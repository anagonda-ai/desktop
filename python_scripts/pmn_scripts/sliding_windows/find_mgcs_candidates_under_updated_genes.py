import csv
from collections import defaultdict

# Read the CSV file
input_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/combined_updated_genes.csv'
output_file = '/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/candidates.csv'

# Data structure to hold the pathways and their EC numbers
pathways = defaultdict(lambda: defaultdict(set))

# Read the input CSV file
with open(input_file, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        pathway = row['Pathway']
        gene_id = row['Query ID']
        ec_number = row['Predicted EC number']
        family = ec_number.split('.')[0]
        subfamily = '.'.join(ec_number.split('.')[:2])
        pathways[pathway]['families'].add(family)
        pathways[pathway]['subfamilies'].add(subfamily)
        pathways[pathway]['genes'].add((gene_id, ec_number))

# Filter pathways based on the criteria
filtered_pathways = {}
for pathway, data in pathways.items():
    if len(data['families']) >= 2:
        if len(data['subfamilies']) >= 3:
            filtered_pathways[pathway] = data
        elif len(data['subfamilies']) == 2:
            subfamily_families = defaultdict(set)
            for gene_id, ec_number in data['genes']:
                family = ec_number.split('.')[0]
                subfamily = '.'.join(ec_number.split('.')[:2])
                subfamily_families[subfamily].add(family)
            if all(len(families) > 1 for families in subfamily_families.values()):
                filtered_pathways[pathway] = data

# Write the filtered pathways to a new CSV file
with open(output_file, mode='w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Pathway', 'Gene IDs and Predicted EC numbers'])
    for pathway, data in filtered_pathways.items():
        gene_ec_list = [f"{gene_id}:{ec_number}" for gene_id, ec_number in data['genes']]
        writer.writerow([pathway, '; '.join(gene_ec_list)])

print(f"Filtered pathways have been written to {output_file}")