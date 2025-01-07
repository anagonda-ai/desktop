import os
import subprocess
from Bio import SeqIO
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

def read_fasta(fasta_file):
    """Reads a FASTA file and returns a DataFrame with identifiers and sequences."""
    print(f"Reading FASTA file: {fasta_file}")
    genes = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.id
        sequence = str(record.seq)
        genes.append({"qseqid": header, "sequence": sequence})
    print(f"Loaded {len(genes)} genes from FASTA file.")
    return pd.DataFrame(genes)

def read_blast_scores(blast_file):
    """Reads a BLAST score file and returns a DataFrame."""
    print(f"Reading BLAST scores file: {blast_file}")
    blast_df = pd.read_csv(blast_file, sep='\t')
    print(f"Loaded {len(blast_df)} rows of BLAST scores.")
    return blast_df

def create_sliding_windows(genes_df, window_size, blast_df):
    """Creates sliding windows with unique genes selected based on highest BLAST score."""
    print(f"Creating sliding windows with size {window_size}.")
    windows = []
    current_window = []
    seen_genes = set()

    genes_with_scores = pd.merge(genes_df, blast_df[['qseqid', 'normalized_composite_score']], on='qseqid', how='left')
    genes_with_scores['normalized_composite_score'] = genes_with_scores['normalized_composite_score'].fillna(0)

    for _, row in genes_with_scores.iterrows():
        gene_id = row['qseqid'].split('.')[0]

        if gene_id not in seen_genes:
            same_gene_transcripts = genes_with_scores[genes_with_scores['qseqid'].str.startswith(gene_id)]
            best_transcript = same_gene_transcripts.loc[same_gene_transcripts['normalized_composite_score'].idxmax()]
            current_window.append(best_transcript)
            seen_genes.add(gene_id)

        if len(current_window) == window_size:
            windows.append(pd.DataFrame(current_window).drop(columns=['normalized_composite_score'], errors='ignore'))
            current_window = []
            seen_genes = set()

    if current_window:
        windows.append(pd.DataFrame(current_window).drop(columns=['normalized_composite_score'], errors='ignore'))

    print(f"Created {len(windows)} sliding windows.")
    return windows

def filter_single_window(window, blast_df, score_threshold, min_genes_above_threshold):
    """Filters a single window based on score threshold and minimum gene count."""
    merged_window = pd.merge(window, blast_df, on="qseqid", how="left")
    genes_above_threshold = merged_window[merged_window["normalized_composite_score"] > score_threshold]
    return merged_window if len(genes_above_threshold) >= min_genes_above_threshold else None

def filter_windows_parallel(windows, blast_df, score_threshold, min_genes_above_threshold, num_workers):
    """Filters windows in parallel."""
    print(f"Filtering {len(windows)} windows in parallel with {num_workers} workers.")
    filtered = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(filter_single_window, window, blast_df, score_threshold, min_genes_above_threshold): idx
            for idx, window in enumerate(windows)
        }
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                filtered.append(result)
    print(f"Filtered down to {len(filtered)} windows.")
    return filtered

def blast_window_against_mgc(window, mgc_fasta_path, window_dir, window_idx, mgc_id):
    """Blast a window against a single MGC-specific FASTA file."""
    window_sequence = "".join(window["sequence"].values)
    window_fasta = os.path.join(window_dir, f"window_{window_idx}_query.fasta")
    blast_output = os.path.join(window_dir, f"window_{window_idx}_vs_{mgc_id}_blast_results.txt")

    # Save the window as a FASTA file
    with open(window_fasta, "w") as fasta_out:
        fasta_out.write(f">window_{window_idx}_query\n{window_sequence}\n")

    # Run BLAST
    blast_command = [
        "blastp",  # Use "blastn" if nucleotide sequences
        "-query", window_fasta,
        "-subject", mgc_fasta_path,
        "-out", blast_output,
        "-outfmt", "6",
        "-evalue", "1e-5",
        "-max_target_seqs", "10"
    ]
    subprocess.run(blast_command)
    print(f"BLAST completed for window {window_idx} against MGC {mgc_id}.")
    return blast_output

def split_mgc_fasta(mgc_fasta, output_dir):
    """Splits the MGC FASTA file by MGC IDs and saves each group as a separate FASTA file."""
    print(f"Splitting MGC FASTA file: {mgc_fasta}")
    os.makedirs(output_dir, exist_ok=True)
    mgc_groups = {}

    for record in SeqIO.parse(mgc_fasta, "fasta"):
        # Extract MGC ID (assuming it's part of the header)
        mgc_id = record.id.split("|")[0]  # Adjust based on actual header format
        if mgc_id not in mgc_groups:
            mgc_groups[mgc_id] = []
        mgc_groups[mgc_id].append(record)

    for mgc_id, sequences in mgc_groups.items():
        mgc_fasta_path = os.path.join(output_dir, f"{mgc_id}.fasta")
        with open(mgc_fasta_path, "w") as out_fasta:
            SeqIO.write(sequences, out_fasta, "fasta")
        print(f"Saved {len(sequences)} sequences for MGC {mgc_id}.")

    return mgc_groups.keys()

def blast_windows_against_mgcs(filtered_windows, mgc_dir, output_dir):
    """Blast each filtered window against each MGC-specific FASTA file, saving results in per-window directories."""
    print(f"Blasting {len(filtered_windows)} windows against MGC-specific FASTA files in {mgc_dir}.")
    os.makedirs(output_dir, exist_ok=True)
    mgc_fastas = [f for f in os.listdir(mgc_dir) if f.endswith(".fasta")]
    results_summary = []

    for window_idx, window in enumerate(filtered_windows):
        # Create a directory for the current window
        window_dir = os.path.join(output_dir, f"window_{window_idx}")
        os.makedirs(window_dir, exist_ok=True)

        for mgc_fasta in mgc_fastas:
            mgc_id = os.path.splitext(mgc_fasta)[0]
            mgc_fasta_path = os.path.join(mgc_dir, mgc_fasta)
            blast_output = blast_window_against_mgc(window, mgc_fasta_path, window_dir, window_idx, mgc_id)
            results_summary.append({"window_idx": window_idx, "mgc_id": mgc_id, "result_path": blast_output})

    # Save the results summary to a CSV file
    results_file = os.path.join(output_dir, "blast_results_summary.csv")
    results_df = pd.DataFrame(results_summary)
    results_df.to_csv(results_file, index=False)
    print(f"Saved results summary to {results_file}.")
    return results_summary

def main(fasta_file, blast_file, mgc_dir, window_size, score_threshold, min_genes_above_threshold, num_workers):
    print("Starting process...")
    genes_df = read_fasta(fasta_file)
    blast_scores = read_blast_scores(blast_file)
    windows = create_sliding_windows(genes_df, window_size, blast_scores)
    filtered_windows = filter_windows_parallel(windows, blast_scores, score_threshold, min_genes_above_threshold, num_workers)

    output_dir = os.path.splitext(fasta_file)[0] + "_windows"
    os.makedirs(output_dir, exist_ok=True)
    
    # Split the MGC FASTA file
    mgc_dir = os.path.join(output_dir, "MGC_FASTAs")
    split_mgc_fasta(mgc_fasta, mgc_dir)

    # Blast windows against each MGC-specific FASTA
    blast_results_dir = os.path.join(output_dir, "BLAST_results")
    blast_windows_against_mgcs(filtered_windows, mgc_dir, blast_results_dir)
    
    

    print("Process completed successfully!")


# פרמטרים
fasta_file = "/groups/itay_mayrose_nosnap/alongonda/desktop/arabidopsis/sorted_Arabidopsis_thaliana.TAIR10.pep.all.fa"
blast_file = "/groups/itay_mayrose_nosnap/alongonda/desktop/Arabidopsis_thaliana.TAIR10.pep.ordered.fa/best_normalized_blast_scores.csv"
mgc_fasta = "/groups/itay_mayrose_nosnap/alongonda/desktop/MGCs/all_genes_from_mibig/mibig_prot_seqs_4.0.fasta"  # Replace with the path to your MGCs FASTA file

window_size = 10  # מספר גנים בחלון
score_threshold = 0.8  # סף ציון
min_genes_above_threshold = 3  # מספר מינימלי של גנים עם ציון מעל הסף
num_workers = 8  # מספר תהליכים מקבילים (לפי מספר הליבות הזמינות)

# קריאה לפונקציית main
main(fasta_file, blast_file, mgc_fasta, window_size, score_threshold, min_genes_above_threshold, num_workers)