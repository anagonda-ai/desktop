import os
import subprocess
import time
import pandas as pd
from threading import Lock
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from Bio.Blast import NCBIXML

# ------------------------------
# KEGG Annotation Function (Local or within SLURM)
# ------------------------------

def blast_and_map_to_kegg(genome_file, kegg_db, temp_dir, identity_threshold=90.0, evalue_threshold=1e-3, coverage_threshold=90.0):
    df = pd.read_csv(genome_file)
    fasta_query = os.path.join(temp_dir, os.path.basename(genome_file).replace('.csv', '.fasta'))

    # Write fasta
    with open(fasta_query, "w") as f:
        for idx, row in df.iterrows():
            f.write(f">{row['id']}\n{row['sequence']}\n")

    blast_output = os.path.join(temp_dir, f"blast_result_{os.path.basename(genome_file)}.xml")
    if not os.path.exists(blast_output):
        subprocess.run([
            "blastp",
            "-task", "blastp-fast",
            "-query", fasta_query,
            "-db", kegg_db,
            "-out", blast_output,
            "-outfmt", "5",
            "-evalue", str(evalue_threshold),
            "-num_threads", "4",
            "-max_target_seqs", "5"
        ], check=True)

    id_to_pathway, id_to_annotation = {}, {}
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
                    _, annotation, pathway = header.split('$')
                    annotation = annotation.strip()
                    pathway = pathway.strip()
                else:
                    continue

                for hsp in alignment.hsps:
                    coverage = (hsp.align_length / query_length) * 100
                    identity = (hsp.identities / hsp.align_length) * 100
                    bit_score = hsp.bits

                    if coverage >= coverage_threshold and identity >= identity_threshold and bit_score > best_score:
                        best_score = bit_score
                        best_pathway = pathway
                        best_annotation = annotation

            if best_pathway:
                id_to_pathway[query_id] = best_pathway
                id_to_annotation[query_id] = best_annotation

    df['pathway'] = df['id'].map(id_to_pathway)
    df['annotation'] = df['id'].map(id_to_annotation)
    return df

# ------------------------------
# Count Running Jobs
# ------------------------------

def count_running_jobs(user="alongonda"):
    """Count number of running jobs for a user."""
    result = subprocess.run(
        ["squeue", "-u", user, "-h"],
        capture_output=True, text=True
    )
    jobs = result.stdout.strip().splitlines()
    return len(jobs)

# ------------------------------
# SLURM Submission
# ------------------------------

def submit_annotation_jobs(genome_files, kegg_db, temp_dir, annotated_dir, max_jobs=100):
    slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/kegg_annotation.sh"
    os.makedirs(annotated_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    for genome_file in genome_files:
        job_name = os.path.basename(genome_file).replace('.csv', '')

        output_file = os.path.join(annotated_dir, os.path.basename(genome_file).replace('.csv', '_annotated.csv'))
        if os.path.exists(output_file):
            print(f"✔️ Already exists: {output_file}")
            continue
        
        # Throttle: wait if too many jobs
        while True:
            current_jobs = count_running_jobs(user="alongonda")
            if current_jobs < max_jobs:
                break
            print(f"⏳ {current_jobs} jobs running (limit {max_jobs}). Waiting 30 seconds...")
            time.sleep(30)
        
        sbatch_cmd = [
            "sbatch",
            "--job-name", f"kegg_{job_name}",
            "--output", os.path.join(annotated_dir, f"{job_name}.out"),
            "--error", os.path.join(annotated_dir, f"{job_name}.err"),
            slurm_script,
            genome_file,
            annotated_dir,
            temp_dir,
            kegg_db
        ]
        result = subprocess.run(sbatch_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✅ Submitted job for {job_name}: {result.stdout.strip()}")
        else:
            print(f"❌ Failed to submit job for {job_name}: {result.stderr.strip()}")

# ------------------------------
# Wait for SLURM Jobs
# ------------------------------

def wait_for_jobs(job_prefix="kegg_", user="alongonda", check_interval=60):
    print(f"⏳ Waiting for SLURM jobs starting with '{job_prefix}' by user '{user}'...")
    while True:
        result = subprocess.run(
            ["squeue", "-u", user, "-h", "-o", "%.18i %.9P %.100j %.8u %.2t %.10M %.6D %R"],
            capture_output=True, text=True
        )
        jobs = result.stdout.strip().splitlines()
        kegg_jobs = [j for j in jobs if job_prefix in j]

        if not kegg_jobs:
            print("✅ All KEGG annotation jobs are finished.")
            break
        else:
            print(f"🕒 {len(kegg_jobs)} jobs still running. Checking again in {check_interval} seconds...")
            time.sleep(check_interval)
            
# ------------------------------
# Sliding Window Clustering
# ------------------------------

def find_first_common_element(genes_and_pathways, min_genes):
    from itertools import combinations, product

    gene_sets = {gene: set([pathway for pathway in paths if pd.notna(pathway)])
                 for gene, paths in genes_and_pathways.items()}

    for group_size in range(len(gene_sets), min_genes - 1, -1):
        for selected_genes in combinations(gene_sets.keys(), group_size):
            selected_sets = [gene_sets[gene] for gene in selected_genes]
            for combination in product(*selected_sets):
                if len(set(combination)) == 1:
                    return combination[0], list(selected_genes)
    return None, []

def process_annotated_file(annotated_file, output_file, file_lock, window_size, min_genes):
    print(f"Processing: {annotated_file}")
    total_matches = 0
    df = pd.read_csv(annotated_file)
    df["index"] = df.index
    filtered_df = df[df["pathway"].notna()]
    chromosomes = filtered_df['chromosome'].unique()
    prev_matches = set()

    for chromosome in chromosomes:
        chr_data = filtered_df[filtered_df['chromosome'] == chromosome]
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
                genes_and_pathways = {row['id']: [row['pathway']] for _, row in window_df.iterrows()}
                genes_and_annotations = {row['id']: [row['annotation']] for _, row in window_df.iterrows()}
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

    print(f"✔️ {annotated_file} - Matches: {total_matches}")
    return total_matches

# ------------------------------
# Subset Filtering
# ------------------------------

def remove_subset_results(output_file):
    df = pd.read_csv(output_file)
    df['metabolic_genes_set'] = df['metabolic_genes'].apply(lambda x: set(x.split(',')))
    to_remove = set()
    for i, genes_i in enumerate(df['metabolic_genes_set']):
        for j, genes_j in enumerate(df['metabolic_genes_set']):
            if i != j and genes_i < genes_j:
                to_remove.add(i)
                break
    filtered_df = df.drop(list(to_remove)).drop(columns=['metabolic_genes_set'])
    filtered_output_file = output_file.replace(".csv", "_filtered.csv")
    filtered_df.to_csv(filtered_output_file, index=False)
    print(f"🧹 Removed {len(to_remove)} subset results: {filtered_output_file}")

# ------------------------------
# Main Function
# ------------------------------

def main():
    # Directories
    full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
    genome_dirs = [
        os.path.join(full_genome_dir, "ensembl/processed_annotations_test_no_chloroplast_with_sequences"),
        os.path.join(full_genome_dir, "phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences_with_strand"),
        os.path.join(full_genome_dir, "plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences")
    ]
    min_genes = 3
    kegg_db = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/fasta/merged_metabolic_pathways/merged_metabolic_pathways"
    head_output_dir = f"/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g{min_genes}_slurm"
    output_dir = os.path.join(head_output_dir, "kegg_scanner_min_genes_based_metabolic")
    temp_dir = os.path.join(head_output_dir, "blast_temp_annotated_metabolic")
    annotated_dir = os.path.join(head_output_dir, "annotated_genomes_metabolic")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(annotated_dir, exist_ok=True)

    # Genome files
    genome_files = []
    for genome_dir in genome_dirs:
        for file in os.listdir(genome_dir):
            if file.endswith('.csv'):
                genome_files.append(os.path.join(genome_dir, file))

    print(f"Found {len(genome_files)} genome files.")

    # Submit SLURM jobs for annotation
    submit_annotation_jobs(genome_files, kegg_db, temp_dir, annotated_dir, max_jobs=100)
    
    # Wait until all jobs finish
    wait_for_jobs(job_prefix="kegg_", user="alongonda")

    # Sliding Window (after annotation completes manually)
    annotated_files = [os.path.join(annotated_dir, f) for f in os.listdir(annotated_dir) if f.endswith('_annotated.csv')]

    min_genes_subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    os.makedirs(min_genes_subdir, exist_ok=True)

    for window_size in [10, 20]:
        output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}.csv")
        if os.path.exists(output_file):
            os.remove(output_file)

        file_lock = Lock()
        total_matches = 0

        with tqdm(total=len(annotated_files), desc=f"Sliding Window w{window_size}", unit='file') as pbar:
            with ThreadPoolExecutor(max_workers=32) as executor:
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

        print(f"✔️ Total Matches for w{window_size}, min_genes={min_genes}: {total_matches}")
        remove_subset_results(output_file)


if __name__ == "__main__":
    main()
