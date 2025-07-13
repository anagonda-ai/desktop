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
from concurrent.futures import ThreadPoolExecutor
import threading

annotated_file = "$annotated_file"
output_file = "$output_file"
window_size = int("$window_size")
min_genes = int("$min_genes")
file_lock = threading.Lock()

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
    prev_matches = set()
    last_match = None
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
            if pathway and not tuple(metabolic_genes) in prev_matches:
                prev_matches.add(tuple(metabolic_genes))
                metabolic_annotations = [genes_and_annotations[gene][0] for gene in metabolic_genes]

                current_match = {
                    'pathway': pathway,
                    'genes': list(window_df['id']),
                    'metabolic_genes': metabolic_genes,
                    'metabolic_genes_annotations': metabolic_annotations,
                    'start': window_df['start'].min(),
                    'end': window_df['end'].max()
                }

                if (
                    last_match and
                    pathway == last_match['pathway'] and
                    metabolic_genes[:2] == last_match['metabolic_genes'][-2:]
                ):
                    new_genes = [g for g in current_match['genes'] if g not in last_match['genes']]
                    new_metabolic_genes = [g for g in current_match['metabolic_genes'] if g not in last_match['metabolic_genes']]
                    new_annotations = [
                        ann for gene, ann in zip(current_match['metabolic_genes'], current_match['metabolic_genes_annotations'])
                        if gene in new_metabolic_genes
                    ]
                    last_match['genes'] += new_genes
                    last_match['metabolic_genes'] += new_metabolic_genes
                    last_match['metabolic_genes_annotations'] += new_annotations
                    last_match['end'] = max(last_match['end'], current_match['end'])
                else:
                    if last_match:
                        matches.append({
                            'pathway': last_match['pathway'],
                            'genes': ','.join(last_match['genes']),
                            'metabolic_genes': ','.join(last_match['metabolic_genes']),
                            'metabolic_genes_annotations': ','.join(last_match['metabolic_genes_annotations']),
                            'start': last_match['start'],
                            'end': last_match['end'],
                            'source_file': annotated_file
                        })
                    last_match = current_match
        i += 1

    if last_match:
        matches.append({
            'pathway': last_match['pathway'],
            'genes': ','.join(last_match['genes']),
            'metabolic_genes': ','.join(last_match['metabolic_genes']),
            'metabolic_genes_annotations': ','.join(last_match['metabolic_genes_annotations']),
            'start': last_match['start'],
            'end': last_match['end'],
            'source_file': annotated_file
        })

    return matches

def process_annotated_file(annotated_file, output_file, window_size, min_genes):
    print(f"Processing: {annotated_file}")
    df = pd.read_csv(annotated_file)
    df["index"] = df.index
    filtered_df = df[df["pathway"].notna()]
    chromosomes = filtered_df['chromosome'].unique()

    def wrapped(chr):
        chr_data = filtered_df[filtered_df['chromosome'] == chr]
        return process_chromosome(chr_data, window_size, min_genes, annotated_file)

    with ThreadPoolExecutor(max_workers=16) as executor:
        results = executor.map(wrapped, chromosomes)

    all_matches = []
    for match_list in results:
        all_matches.extend(match_list)

    if all_matches:
        with file_lock:
            mode = 'a' if os.path.exists(output_file) else 'w'
            header = not os.path.exists(output_file)
            pd.DataFrame(all_matches).to_csv(output_file, mode=mode, header=header, index=False)
            print(f"✔️ {annotated_file} - Matches: {len(all_matches)}")
    else:
        print(f"⚠️ {annotated_file} - No matches found")

process_annotated_file(annotated_file, output_file, window_size, min_genes)
EOF
