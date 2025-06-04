import os
import pandas as pd
import subprocess
from Bio.Blast import NCBIXML
from tqdm import tqdm
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional, Dict, List
import time
from pathlib import Path

# === Setup Logging ===
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(os.path.dirname(__file__), 'blast_pipeline.log'))
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class BlastHit:
    mibig_gene: str
    matched_gene: str
    source_file: str
    index_in_file: int
    bit_score: float
    identity: float
    coverage: float
    combined_score: float

    def to_dict(self):
        return self.__dict__

class BlastPipeline:
    def __init__(self):
        self.start_time = time.time()
        self._setup_directories()
        
    def _setup_directories(self):
        """Create necessary directories if they don't exist"""
        for dir_path in [OUTPUT_BASE, OUTPUT_DIR]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

# === CONFIG ===
GENOME_CSV_DIRS = [
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test_full_data",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations_with_chromosomes"
]
MIBIG_QUERY_FASTA = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files/merged_metabolic_pathways/merged_metabolic_pathways.fasta"
MGC_TO_ORGANISM_FILE = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/gbk_files/organisms.csv"

OUTPUT_BASE = "/groups/itay_mayrose/alongonda/Plant_MGC/mibig_metabolic_output/annotated_genomes_metabolic/merged_annotation_with_chloroplast"

MERGED_FASTA = os.path.join(OUTPUT_BASE, "merged_genomes.fasta")
BLAST_DB = os.path.join(OUTPUT_BASE, "merged_genomes_blastdb")
BLAST_OUT = os.path.join(OUTPUT_BASE, "blast_results.xml")
OUTPUT_DIR = os.path.join(OUTPUT_BASE, "split_by_mgc")
FINAL_OUTPUT = os.path.join(OUTPUT_BASE, "merged_annotation_id_eq_annotation.csv")
DISTANCES_OUTPUT = os.path.join(OUTPUT_BASE, "mgc_index_distances.csv")
# === Constants ===
EVALUE_THRESHOLD = 1e-3
NUM_THREADS = 32


# === STEP 1: Build merged FASTA with concurrency ===
def process_genome_csv(filepath):
    records = []
    try:
        df = pd.read_csv(filepath)
        for idx, row in df.iterrows():
            header = f"{row['id']}|{os.path.basename(filepath)}|{idx}"
            records.append((header, row["sequence"]))
    except Exception as e:
        print(f"❌ Failed on {filepath}: {e}")
    return records

def build_merged_fasta():
    if os.path.exists(MERGED_FASTA):
        print(f"✅ Merged FASTA already exists: {MERGED_FASTA}")
        return

    print("📦 Building merged FASTA with concurrency...")
    genome_csvs = []
    for genome_dir in GENOME_CSV_DIRS:
        for file in os.listdir(genome_dir):
            if file.endswith(".csv"):
                genome_csvs.append(os.path.join(genome_dir, file))

    records = []
    with ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = [executor.submit(process_genome_csv, f) for f in genome_csvs]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Merging FASTA"):
            records.extend(future.result())

    with open(MERGED_FASTA, 'w') as out:
        for header, seq in records:
            out.write(f">{header}\n{seq}\n")

    print(f"✅ Wrote merged FASTA: {MERGED_FASTA}")


# === STEP 2: BLAST DB ===
def build_blast_db():
    # Check for any BLAST DB files (including db.00.* for multi-volume DBs)
    db_files_exist = any(
        f.startswith(os.path.basename(BLAST_DB)) and (
            f.endswith('.pin') or f.endswith('.phr') or f.endswith('.psq') or
            f.endswith('.00.pin') or f.endswith('.00.phr') or f.endswith('.00.psq')
        )
        for f in os.listdir(os.path.dirname(BLAST_DB))
    )
    if db_files_exist:
        print(f"✅ BLAST DB already exists: {BLAST_DB}")
        return

    print("🛠️  Building BLAST DB...")
    subprocess.run([
        "makeblastdb",
        "-in", MERGED_FASTA,
        "-dbtype", "prot",
        "-out", BLAST_DB
    ], check=True)
    print("✅ BLAST DB created.")


# === STEP 3: Run BLAST ===
def run_blast():
    if os.path.exists(BLAST_OUT):
        print(f"✅ BLAST results already exist: {BLAST_OUT}")
        return

    print("🚀 Running BLASTP...")
    subprocess.run([
        "blastp",
        "-query", MIBIG_QUERY_FASTA,
        "-db", BLAST_DB,
        "-out", BLAST_OUT,
        "-outfmt", "5",
        "-evalue", str(EVALUE_THRESHOLD),
        "-num_threads", str(NUM_THREADS),
        "-max_target_seqs", "5"
    ], check=True)
    print("✅ BLASTP complete.")


# === STEP 4: Parse BLAST XML (parallel) ===
def extract_hits_per_file(record):
    query_id = record.query
    query_len = record.query_length
    hits_by_file = {}

    for alignment in record.alignments:
        try:
            gene_id, filepath, index = alignment.hit_def.split("|")
            complete_filepath = next(
                (os.path.join(dir_path, filepath) for dir_path in GENOME_CSV_DIRS
                 if os.path.exists(os.path.join(dir_path, filepath))),
                None
            )
            if complete_filepath is None:
                continue
        except ValueError:
            continue

        for hsp in alignment.hsps:
            identity = (hsp.identities / hsp.align_length) * 100
            coverage = (hsp.align_length / query_len) * 100
            combined_score = (
                (hsp.bits * 0.4) +
                (identity * 0.3) +
                (coverage * 0.3)
            )

            current = hits_by_file.get(complete_filepath)
            if current is None or combined_score > current["combined_score"]:
                hits_by_file[complete_filepath] = {
                    "mibig_gene": query_id,
                    "matched_gene": gene_id,
                    "source_file": complete_filepath,
                    "index_in_file": int(index),
                    "bit_score": hsp.bits,
                    "identity": identity,
                    "coverage": coverage,
                    "combined_score": combined_score
                }

    return list(hits_by_file.values())

def get_best_source_file(group, original_gene_count):
    """
    Get the best source file by:
    1. First selecting files with maximum number of matches
    2. Among those, choose the one with best combined score
    """
    source_stats = group.groupby('source_file').agg({
        'combined_score': ['count', 'sum', 'mean']
    })
    source_stats.columns = ['match_count', 'total_score', 'avg_score']
    
    # Get the maximum number of matches
    max_matches = source_stats['match_count'].max()
    
    # Filter to only sources with maximum matches
    best_match_sources = source_stats[source_stats['match_count'] == max_matches]
    
    # Among those with max matches, choose the one with highest average score
    return best_match_sources['avg_score'].idxmax()

def parse_blast_results():
    logger.info("Starting BLAST results parsing")
    pipeline = BlastPipeline()

    try:
        with open(BLAST_OUT) as handle:
            records = list(NCBIXML.parse(handle))
    except Exception as e:
        logger.error(f"Failed to parse BLAST XML: {str(e)}")
        raise

    all_hits: List[BlastHit] = []
    failed_records = 0

    with ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = [executor.submit(extract_hits_per_file, rec) for rec in records]
        with tqdm(total=len(futures), desc="Parsing BLAST hits") as pbar:
            for future in as_completed(futures):
                try:
                    hits = future.result()
                    all_hits.extend([BlastHit(**hit) for hit in hits])
                    pbar.update(1)
                except Exception as e:
                    failed_records += 1
                    logger.error(f"Failed to process record: {str(e)}")
                    pbar.update(1)

    if failed_records:
        logger.warning(f"Failed to process {failed_records} records")

    final_output = pd.DataFrame([hit.to_dict() for hit in all_hits])

    final_output["mgc_id"] = final_output["mibig_gene"].str.extract(r"\|\s*(BGC\d+)\s*\|")

    original_gene_counts = final_output.groupby('mgc_id')['mibig_gene'].nunique()
    all_mgc_genes = {mgc: set(group['mibig_gene'].unique()) 
                     for mgc, group in final_output.groupby('mgc_id')}

    mgc_best_sources = {}
    mgc_max_matches = {}

    for mgc_id, group in final_output.groupby("mgc_id"):
        source_counts = group.groupby('source_file').size()
        max_matches = source_counts.max()
        mgc_max_matches[mgc_id] = max_matches

        best_source = get_best_source_file(group, original_gene_counts[mgc_id])
        mgc_best_sources[mgc_id] = best_source

    final_output = final_output.apply(
        lambda row: row if row['source_file'] == mgc_best_sources[row['mgc_id']] else None,
        axis=1
    ).dropna()

    final_output.to_csv(FINAL_OUTPUT, index=False)
    logger.info(f"Saved {len(final_output)} best hits to: {FINAL_OUTPUT}")

    for mgc_id, group in final_output.groupby("mgc_id"):
        output_path = os.path.join(OUTPUT_DIR, f"{mgc_id}.csv")
        group.drop(columns=["mgc_id"]).to_csv(output_path, index=False)
        logger.debug(f"Saved individual MGC file: {output_path}")

    mgc_stats = []
    for mgc_id, group in final_output.groupby("mgc_id"):
        source_file = mgc_best_sources[mgc_id]
        original_genes = all_mgc_genes[mgc_id]
        matched_genes = set(group['mibig_gene'].unique())
        missing_genes = original_genes - matched_genes

        stats = {
            'mgc_id': mgc_id,
            'original_gene_count': original_gene_counts[mgc_id],
            'best_source_file': source_file,
            'num_matches': len(group),
            'max_matches_found': mgc_max_matches[mgc_id],
            'completeness': len(group) / original_gene_counts[mgc_id],
            'total_score': group['combined_score'].sum(),
            'avg_score': group['combined_score'].mean(),
            'max_score': group['combined_score'].max(),
            'missing_genes': '; '.join(sorted(missing_genes)) if missing_genes else 'None'
        }
        mgc_stats.append(stats)

    stats_df = pd.DataFrame(mgc_stats)
    stats_df = stats_df[['mgc_id', 'original_gene_count', 'best_source_file', 'num_matches', 'max_matches_found', 
                         'completeness', 'total_score', 'avg_score', 'max_score', 'missing_genes']]
    stats_df.to_csv(os.path.join(OUTPUT_BASE, "mgc_best_sources_stats.csv"), index=False)
    logger.info("Saved best source files summary with statistics")

    missing_genes_df = stats_df[['mgc_id', 'missing_genes']]
    missing_genes_df = missing_genes_df[missing_genes_df['missing_genes'] != 'None']
    missing_genes_df.to_csv(os.path.join(OUTPUT_BASE, "missing_genes.csv"), index=False)
    logger.info("Saved missing genes summary")

    elapsed_time = time.time() - pipeline.start_time
    logger.info(f"Pipeline completed in {elapsed_time:.2f} seconds")
    logger.info(f"Processed {len(records)} records with {len(final_output)} retained hits")
    return final_output


def distance_from_lowest_to_highest_index(group):
    indices = group["index_in_file"].astype(int)
    if not indices.empty:
        return indices.max() - indices.min()
    return 0

def save_distances_to_csv(distances, output_file):
    df = pd.DataFrame(distances)
    df.to_csv(output_file, index=False)
    print(f"✅ MGC index distances saved to: {output_file}")

def load_final_output_and_save_distances_to_csv():
    if not os.path.exists(FINAL_OUTPUT):
        print(f"❌ Final output file does not exist: {FINAL_OUTPUT}")
        return

    final_output = pd.read_csv(FINAL_OUTPUT)
    
    # Calculate distance from lowest to highest index gene in each MGC
    mgc_distances = []
    for mgc_id, group in final_output.groupby("mgc_id"):
        distance = distance_from_lowest_to_highest_index(group)
        if distance > 0:
            mgc_distances.append({"mgc_id": mgc_id, "distance": distance})
            
    # Save distances to CSV
    save_distances_to_csv(mgc_distances, DISTANCES_OUTPUT)

# === MAIN ===
if __name__ == "__main__":
    build_merged_fasta()
    build_blast_db()
    run_blast()
    parse_blast_results()
    load_final_output_and_save_distances_to_csv()
