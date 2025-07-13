#!/bin/bash
#SBATCH --job-name=kegg_annot
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=itaym
#SBATCH --time=3-24:00:00

annotated_file=$1
output_file=$2
window_size=$3
min_genes=$4

python3 - <<EOF
import pandas as pd
import os
from itertools import combinations, product
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

annotated_file = "$annotated_file"
output_file = "$output_file"
window_size = int("$window_size")
min_genes = int("$min_genes")

def find_first_common_element(genes_and_pathways, min_genes):
    gene_sets = {gene: set([pathway for pathway in paths if pd.notna(pathway)])
                 for gene, paths in genes_and_pathways.items()}

    for group_size in range(len(gene_sets), min_genes - 1, -1):
        for selected_genes in combinations(gene_sets.keys(), group_size):
            selected_sets = [gene_sets[gene] for gene in selected_genes]
            for combination in product(*selected_sets):
                if len(set(combination)) == 1:
                    return combination[0], list(selected_genes)
    return None, []

def process_chromosome(chr_data, window_size, min_genes, annotated_file):
    matches = []
    num_genes = len(chr_data)
    i = 0
    while i < num_genes:
        window = [chr_data.iloc[i]]
        start_index = chr_data.iloc[i]['index']
        for j in range(i + 1, num_genes):
            end_index = chr_data.iloc[j]['index']
            if (end_index - start_index <= window_size):
                window.append(chr_data.iloc[j])
            else:
                break

        if len(window) >= min_genes:
            window_df = pd.DataFrame(window)
            genes_and_pathways = {row['id']: row['pathway'].split(",") for _, row in window_df.iterrows()}
            genes_and_annotations = {row['id']: row['annotation'].split(",") for _, row in window_df.iterrows()}
            pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_genes)
            if pathway:
                metabolic_annotations = [genes_and_annotations[gene][0] for gene in metabolic_genes]
                group = {
                    'pathway': pathway,
                    'genes': ','.join(window_df['id']),
                    'metabolic_genes': ','.join(metabolic_genes),
                    'metabolic_genes_annotations': ','.join(metabolic_annotations),
                    'start': window_df['start'].min(),
                    'end': window_df['end'].max(),
                    'source_file': annotated_file
                }
                matches.append(group)
        i += 1
    return matches

def process_annotated_file(annotated_file, output_file, window_size, min_genes):
    print(f"Processing: {annotated_file}")
    df = pd.read_csv(annotated_file)
    df["index"] = df.index
    filtered_df = df[df["pathway"].notna()]
    chromosomes = filtered_df['chromosome'].unique()

    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = [
            executor.submit(process_chromosome,
                            filtered_df[filtered_df['chromosome'] == chromosome].copy(),
                            window_size, min_genes, annotated_file)
            for chromosome in chromosomes
        ]
        all_matches = []
        for future in futures:
            result = future.result()
            all_matches.extend(result)

    if all_matches:
        pd.DataFrame(all_matches).drop_duplicates(subset=["metabolic_genes"]).to_csv(
            output_file, index=False
        )
        print(f"✔️ {annotated_file} - Matches: {len(all_matches)}")
    else:
        print(f"⚠️ {annotated_file} - No matches found")

process_annotated_file(annotated_file, output_file, window_size, min_genes)
EOF
