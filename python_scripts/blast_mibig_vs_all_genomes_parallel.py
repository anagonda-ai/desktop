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
    source_organism: str
    query_organism: str

    def to_dict(self):
        return self.__dict__

class BlastPipeline:
    def __init__(self):
        self.start_time = time.time()
        self._setup_directories()
        self._load_mgc_mapping()
        
    def _setup_directories(self):
        """Create necessary directories if they don't exist"""
        for dir_path in [OUTPUT_BASE, OUTPUT_DIR]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
            
    def _load_mgc_mapping(self):
        """Load MGC to organism mapping with validation"""
        self.mgc_to_organism = {}
        if os.path.exists(MGC_TO_ORGANISM_FILE):
            try:
                df_org = pd.read_csv(MGC_TO_ORGANISM_FILE)
                if not all(col in df_org.columns for col in ["MGC", "Organism"]):
                    raise ValueError("Required columns 'MGC' and 'Organism' not found")
                self.mgc_to_organism = dict(zip(df_org["MGC"], df_org["Organism"]))
                logger.info(f"Loaded {len(self.mgc_to_organism)} MGC to organism mappings")
            except Exception as e:
                logger.error(f"Failed to load MGC mapping: {str(e)}")
                raise
        else:
            logger.error(f"Organism file not found: {MGC_TO_ORGANISM_FILE}")
            raise FileNotFoundError(f"Missing required file: {MGC_TO_ORGANISM_FILE}")

# === CONFIG ===
GENOME_CSV_DIRS = [
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"
]
MIBIG_QUERY_FASTA = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files/merged_metabolic_pathways/merged_metabolic_pathways.fasta"
MGC_TO_ORGANISM_FILE = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/gbk_files/organisms.csv"

OUTPUT_BASE = "/groups/itay_mayrose/alongonda/Plant_MGC/mibig_metabolic_output/annotated_genomes_metabolic/merged_annotation"

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
    if all(os.path.exists(f"{BLAST_DB}.{ext}") for ext in ['pin', 'phr', 'psq']):
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
def extract_best_hit(record):
    query_id = record.query
    query_mgc_id = query_id.split("|")[1].strip() if "|" in query_id else query_id
    query_organism = mgc_to_organism.get(query_mgc_id, "Unknown").lower()
    query_len = record.query_length
    
    # Initialize best hit tracking with sophisticated scoring
    best_score = {
        'combined_score': 0,  # Weighted combination of bit score, identity, and coverage
        'data': None
    }

    for alignment in record.alignments:
        hit_id = alignment.hit_def
        try:
            gene_id, filename, index = hit_id.split("|")
        except ValueError:
            continue  # malformed ID

        for hsp in alignment.hsps:
            identity = (hsp.identities / hsp.align_length) * 100
            coverage = (hsp.align_length / query_len) * 100
            source_organism = filename.replace("_"," ").replace("filtered.csv", "").lower().strip()
            
            # Sophisticated scoring system
            # Weight factors for different aspects (can be tuned)
            bit_weight = 0.4
            identity_weight = 0.3
            coverage_weight = 0.3
            organism_bonus = 1.5 if source_organism in query_organism else 1.0
            
            # Combined score calculation
            combined_score = (
                (hsp.bits * bit_weight) +
                (identity * identity_weight) +
                (coverage * coverage_weight)
            ) * organism_bonus

            # Update best hit if this is better
            if combined_score > best_score['combined_score']:
                best_score['combined_score'] = combined_score
                best_score['data'] = {
                    "mibig_gene": query_id,
                    "matched_gene": gene_id,
                    "source_file": filename,
                    "index_in_file": int(index),
                    "bit_score": hsp.bits,
                    "identity": identity,
                    "coverage": coverage,
                    "combined_score": combined_score,
                    "source_organism": source_organism,
                    "query_organism": query_organism
                }

    return best_score['data']


def parse_blast_results():
    logger.info("Starting BLAST results parsing")
    pipeline = BlastPipeline()
    
    try:
        with open(BLAST_OUT) as handle:
            records = list(NCBIXML.parse(handle))
    except Exception as e:
        logger.error(f"Failed to parse BLAST XML: {str(e)}")
        raise

    best_hits: List[BlastHit] = []
    failed_records = 0
    
    with ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = [executor.submit(extract_best_hit, rec) for rec in records]
        
        with tqdm(total=len(futures), desc="Parsing BLAST hits") as pbar:
            for future in as_completed(futures):
                try:
                    hit = future.result()
                    if hit:
                        best_hits.append(BlastHit(**hit))
                    pbar.update(1)
                except Exception as e:
                    failed_records += 1
                    logger.error(f"Failed to process record: {str(e)}")
                    pbar.update(1)

    if failed_records:
        logger.warning(f"Failed to process {failed_records} records")

    # Convert to DataFrame and validate
    final_output = pd.DataFrame([hit.to_dict() for hit in best_hits])
    
    # Ensure we have unique MGCs
    final_output.sort_values('combined_score', ascending=False, inplace=True)
    final_output.drop_duplicates(subset=['mibig_gene'], keep='first', inplace=True)
    
    # Extract MGC ID with validation
    final_output["mgc_id"] = final_output["mibig_gene"].str.extract(r"\|\s*(BGC\d+)\s*\|")
    if final_output["mgc_id"].isnull().any():
        logger.warning("Some records have invalid MGC IDs")
    
    # Save results
    final_output.to_csv(FINAL_OUTPUT, index=False)
    logger.info(f"Saved {len(final_output)} best hits to: {FINAL_OUTPUT}")
    
    # Save individual MGC files
    for mgc_id, group in final_output.groupby("mgc_id"):
        output_path = os.path.join(OUTPUT_DIR, f"{mgc_id}.csv")
        group.drop(columns=["mgc_id"]).to_csv(output_path, index=False)
        logger.debug(f"Saved individual MGC file: {output_path}")
    
    # Log completion statistics
    elapsed_time = time.time() - pipeline.start_time
    logger.info(f"Pipeline completed in {elapsed_time:.2f} seconds")
    logger.info(f"Processed {len(records)} records with {len(best_hits)} successful hits")
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
