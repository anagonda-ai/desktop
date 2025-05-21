import re
from threading import Lock
import os
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from tqdm import tqdm
from Bio.Blast import NCBIXML

# --- KEGG Annotation Functions ---

def blast_and_map_to_kegg(genome_file, kegg_db, temp_dir, identity_threshold=70.0, evalue_threshold=1e-3, coverage_threshold=90.0):
    df = pd.read_csv(genome_file)
    fasta_query = os.path.join(temp_dir, os.path.basename(genome_file).replace('.csv', '.fasta'))
    query_organism_name = os.path.basename(genome_file).split("_filtered")[0].replace("_"," ")

    # Run BLASTP
    blast_output = os.path.join(temp_dir, f"blast_result_{os.path.basename(genome_file)}.xml")
    print(f"Running BLASTP for {genome_file} against KEGG database...")
    if not os.path.exists(blast_output):
        # Save genome sequences to temporary FASTA
        with open(fasta_query, "w") as f:
            for idx, row in df.iterrows():
                f.write(f">{row['id']}\n{row['sequence']}\n")
        subprocess.run([
            "blastp", 
            "-task", "blastp-fast",
            "-query", fasta_query,
            "-db", kegg_db,
            "-out", blast_output,
            "-outfmt", "5",
            "-evalue", str(evalue_threshold),
            "-num_threads", "32",
            "-max_target_seqs", "5",  # <= IMPORTANT: no need to align 1000 hits per gene
        ], check=True)
    print(f"BLASTP completed for {genome_file}.")

    # Parse BLAST results
    id_to_pathway = {}
    id_to_annotation = {}
    if os.path.exists(blast_output):
        print(f"Parsing BLAST results for {blast_output}...")
        with open(blast_output) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                query_id = record.query
                query_length = record.query_length
                best_score = 0
                best_pathway = None
                best_annotation = None
                
                for alignment in record.alignments:
                    header = alignment.hit_def
                    if '$' in header:
                        gene_id_part, annotation_part, pathway_part = header.split('$')
                        pathway_id = pathway_part.strip()
                        annotation = annotation_part.strip()
                        target_organism = re.sub(r'\d', '', pathway_part.strip()).split('_')[0]
                        plants_list = pd.read_csv("/groups/itay_mayrose/alongonda/datasets/KEGG/origin/plants_list.csv")
                        organism_name = plants_list[plants_list["Organism_Code"]==target_organism]["Organism_Name"].iloc[0].lower()
                        
                    else:
                        pathway_id = None
                        annotation = None
                        organism_name = None
                    
                    for hsp in alignment.hsps:
                        coverage = (hsp.align_length / query_length) * 100
                        identity = (hsp.identities / hsp.align_length) * 100
                        bit_score = hsp.bits

                        if query_organism_name in organism_name and coverage >= coverage_threshold and identity >= identity_threshold and bit_score > best_score:
                            best_score = bit_score
                            best_pathway = pathway_id
                            best_annotation = annotation

                if best_pathway:
                    id_to_pathway[query_id] = best_pathway
                    id_to_annotation[query_id] = best_annotation

        # Annotate the dataframe
    df['pathway'] = df['id'].map(id_to_pathway)
    df['annotation'] = df['id'].map(id_to_annotation)
    return df


def annotate_genomes_once(genome_dirs, kegg_db, temp_dir, annotated_dir, max_workers=32):
    if not os.path.exists(annotated_dir):
        os.makedirs(annotated_dir)

    genome_files = []
    for genome_dir in genome_dirs:
        print(f"Searching for genome files in {genome_dir}...", flush=True)
        for file in os.listdir(genome_dir):
            if file.endswith('.csv'):
                genome_files.append(os.path.join(genome_dir, file))

    print(f"Found {len(genome_files)} genome files to annotate.")

    def annotate_single_file(genome_file):
        print(f"Processing file: {genome_file}", flush=True)
        filename = os.path.basename(genome_file).replace('.csv', '_annotated.csv')
        output_path = os.path.join(annotated_dir, filename)

        if os.path.exists(output_path):
            print(f"Annotated file already exists: {output_path}, skipping.", flush=True)
            return

        try:
            print(f"Starting BLAST for {genome_file}...", flush=True)
            df_annotated = blast_and_map_to_kegg(genome_file, kegg_db, temp_dir)
            df_annotated.to_csv(output_path, index=False)
            print(f"Finished BLAST and saved {output_path}", flush=True)
        except subprocess.CalledProcessError as e:
            print(f"BLAST failed for {genome_file}: {e}", flush=True)
        except Exception as e:
            print(f"Unexpected error for {genome_file}: {e}", flush=True)

    with tqdm(total=len(genome_files), desc='Annotating genome files', unit='file') as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(annotate_single_file, genome_file) for genome_file in genome_files]
            for future in as_completed(futures):
                pbar.update(1)


# --- Sliding Window Functions ---

def clean_and_process_genes_and_pathways(genes_and_pathways):
    processed = {}
    for gene, values in genes_and_pathways.items():
        stripped_values = [''.join(c for c in v if c not in "[]'\"").split(",") for v in values]
        processed[gene] = stripped_values
    return processed


def find_first_common_element(genes_and_pathways, min_genes):
    from itertools import combinations, product

    processed = clean_and_process_genes_and_pathways(genes_and_pathways)
    gene_sets = {}
    for gene, paths in processed.items():
        flattened = set(item.strip() for sublist in paths for item in sublist)
        gene_sets[gene] = flattened

    gene_combinations = combinations(gene_sets.keys(), min_genes)
    for selected_genes in gene_combinations:
        selected_sets = [gene_sets[gene] for gene in selected_genes]
        all_combinations = product(*selected_sets)
        for combination in all_combinations:
            if len(set(combination)) == 1:
                return combination[0], list(selected_genes)
    return None, []


def process_annotated_file(annotated_file, output_file, file_lock, window_size, min_genes):
    print(f"Processing file: {annotated_file} with window size: {window_size}")
    total_matches = 0

    df = pd.read_csv(annotated_file)
    df["index"] = df.index
    filtered_df = df[df["pathway"].notna()]
    chromosomes = filtered_df['chromosome'].unique()
    prev_matches = set()

    for chromosome in chromosomes:
        chromosome_data = filtered_df[filtered_df['chromosome'] == chromosome]
        num_genes = len(chromosome_data)
        i = 0
        while i < num_genes:
            window = [chromosome_data.iloc[i]]
            start_index = chromosome_data.iloc[i]['index']
            for j in range(i+1, num_genes):
                end_index = chromosome_data.iloc[j]['index']
                if (end_index - start_index <= window_size):
                    window.append(chromosome_data.iloc[j])
                else:
                    break
            if len(window) >= min_genes:
                window_df = pd.DataFrame(window)
                genes_and_pathways = {row['id']: [row['pathway']] for idx, row in window_df.iterrows()}
                genes_and_annotations = {row['id']: [row['annotation']] for idx, row in window_df.iterrows()}
                pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_genes)
                metabolic_annotations = [genes_and_annotations[gene][0] for gene in metabolic_genes]
                if pathway and not tuple(metabolic_genes) in prev_matches:
                    prev_matches.add(tuple(metabolic_genes))
                    group = {
                        'pathway': pathway,
                        'genes': ','.join(window_df['id']),
                        'metabolic_genes': ','.join(metabolic_genes),
                        'metabolic_genes_annotations': ','.join(metabolic_annotations),
                        'start': window_df['start'].min(),
                        'end': window_df['end'].max(),
                        'source_file': annotated_file
                    }
                    with file_lock:
                        mode = 'a' if os.path.exists(output_file) else 'w'
                        header = not os.path.exists(output_file)
                        pd.DataFrame([group]).to_csv(output_file, mode=mode, header=header, index=False)
                        total_matches += 1
            i += 1
    print(f"Completed file: {annotated_file}, Matches Found: {total_matches}")
    return total_matches


def create_output_subdir(output_dir, min_genes):
    subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    return subdir


def main():
    full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
    # Define genome directories
    genome_dirs = [
        os.path.join(full_genome_dir, "ensembl/processed_annotations_test_no_chloroplast_with_sequences"),
        os.path.join(full_genome_dir, "phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"),
        os.path.join(full_genome_dir, "plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences")
    ]
    kegg_db = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/fasta/merged_metabolic_pathways/merged_metabolic_pathways"
    head_output_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output"
    output_dir = os.path.join(head_output_dir, "kegg_scanner_min_genes_based_metabolic")
    temp_dir = os.path.join(head_output_dir, "blast_temp_annotated_metabolic")
    annotated_dir = os.path.join(head_output_dir, "annotated_genomes_metabolic")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    if not os.path.exists(annotated_dir):
        os.makedirs(annotated_dir)

    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")

    # Stage 1: Annotate genomes (BLAST once)
    annotate_genomes_once(genome_dirs, kegg_db, temp_dir, annotated_dir)

    # Stage 2: Sliding window clustering
    annotated_files = [os.path.join(annotated_dir, f) for f in os.listdir(annotated_dir) if f.endswith('_annotated.csv')]
    print(annotated_files)

    for window_size in [10]:
        for min_genes in [3]:
            min_genes_subdir = create_output_subdir(output_dir, min_genes)
            output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}.csv")
            if os.path.exists(output_file):
                os.remove(output_file)
            file_lock = Lock()
            total_matches = 0
            with tqdm(total=len(annotated_files), desc=f"Sliding Window w{window_size}", unit='file') as pbar:
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = [
                        executor.submit(
                            process_annotated_file,
                            annotated_file,
                            output_file,
                            file_lock,
                            window_size,
                            min_genes
                        )
                        for annotated_file in annotated_files
                    ]
                    for future in as_completed(futures):
                        total_matches += future.result()
                        pbar.update(1)
            print(f"TOTAL MATCHES FOUND for window size {window_size} and min_genes {min_genes}: {total_matches}")
            print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()