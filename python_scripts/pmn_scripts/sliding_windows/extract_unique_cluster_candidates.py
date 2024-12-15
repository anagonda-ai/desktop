import csv
from collections import Counter

def remove_subclusters(clusters):
    unique_clusters = set(clusters)
    for cluster in clusters:
        for other_cluster in clusters:
            if cluster != other_cluster and set(cluster).issubset(set(other_cluster)):
                unique_clusters.discard(cluster)
                break
    return unique_clusters

def normalize_cluster(cluster):
    return tuple(sorted(cluster.split(',')))

def load_csv(file_path):
    clusters = set()
    clusters_with_file_and_pathway = []    
    with open(file_path, mode='r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            if len(row) > 3:
                normalized_cluster = normalize_cluster(row[3])
                pathway = row[0]
                row_path = row[7]
                clusters.add(normalized_cluster)
                for gene in normalized_cluster:
                    clusters_with_file_and_pathway.append([gene, row_path, pathway])
    return clusters, clusters_with_file_and_pathway

def print_clusters(clusters):
        for cluster in clusters:
            print(cluster)

# Example usage
file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/potential_groups_with_window_start_end.csv'
clusters, clusters_with_file_and_pathway = load_csv(file_path)

ensembl_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/ensembl/ensembl_genes.csv'
plaza_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/plaza/plaza_genes.csv'
phytozome_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/phytozome/phytozome_genes.csv'

ensembl_genes = []
plaza_genes = []
phytozome_genes = []

for gene, gene_file_path, pathway in clusters_with_file_and_pathway:
    if 'ensembl' in gene_file_path.lower():
        ensembl_genes.append([gene, gene_file_path, pathway])
    elif 'plaza' in gene_file_path.lower():
        plaza_genes.append([gene, gene_file_path, pathway])
    elif 'phytozome' in gene_file_path.lower():
        phytozome_genes.append([gene, gene_file_path, pathway])

with open(ensembl_file_path, mode='w', newline='') as ensembl_file, open(plaza_file_path, mode='w', newline='') as plaza_file, open(phytozome_file_path, mode='w', newline='') as phytozome_file:
    ensembl_csv_writer = csv.writer(ensembl_file)
    plaza_csv_writer = csv.writer(plaza_file)
    phytozome_csv_writer = csv.writer(phytozome_file)
    for gene, gene_file_path, pathway in clusters_with_file_and_pathway:
        if 'ensembl' in gene_file_path.lower():
            ensembl_csv_writer.writerow([gene, gene_file_path, pathway])
        elif 'plaza' in gene_file_path.lower():
            plaza_csv_writer.writerow([gene, gene_file_path, pathway])
        elif 'phytozome' in gene_file_path.lower():
            phytozome_csv_writer.writerow([gene, gene_file_path, pathway])

unique_clusters = remove_subclusters(clusters)

output_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/unique_clusters_start_end.csv'
with open(file_path, mode='r') as infile, open(output_file_path, mode='w', newline='') as outfile:
    csv_reader = csv.reader(infile)
    csv_writer = csv.writer(outfile)
    for row in csv_reader:
        normalized_cluster = normalize_cluster(row[3])
        if normalized_cluster in unique_clusters:
            csv_writer.writerow(row)
            unique_clusters.remove(normalized_cluster)

