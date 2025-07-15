
class ExtractPromotorMemeProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: extract_promotor_meme.py."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass


from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
import argparse
import pandas as pd
from Bio import SeqIO

# Function to extract dataset_id and organism name from mapping table
def get_organism_info(annotation_file, mapping_df):
    """Extract dataset_id and organism name from mapping table"""
    base = os.path.basename(annotation_file).split(".csv")[0].replace("_annotated", "")
    matches = mapping_df[mapping_df["Original Filename"].str.lower().str.contains(base.lower())]
    if matches.empty:
        raise ValueError(f"❌ No match found in mapping for {annotation_file}")
    row = matches.iloc[0]
    return str(row["Dataset ID"]), row["Organism"]

# Function to find fasta file matching dataset_id or organism name
def find_fasta_file(dataset_id, organism_name, fasta_dirs):
    """Find fasta file matching dataset_id or organism"""
    candidates = []

    for fasta_dir in fasta_dirs:
        if not os.path.isdir(fasta_dir):
            continue
        for file in os.listdir(fasta_dir):
            file_lower = file.lower()
            if (
                str(dataset_id) in file
                or organism_name.lower().replace(" ", "_") in file_lower
            ):
                candidates.append(os.path.join(fasta_dir, file))

    if not candidates:
        raise FileNotFoundError(
            f"❌ No matching fasta file found for {organism_name} (dataset_id={dataset_id})"
        )

    if len(candidates) > 1:
        logger.info(f"⚠️ Multiple fasta candidates found. Using the first:\n{candidates}")

    return candidates[0]

# Function to load genome sequences from a fasta file
def load_genome(fasta_file):
    """Load genome sequences into memory"""
    genome = {}
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            genome[record.id] = record.seq
    return genome

# Function to extract promoter sequences from a cluster dataframe and annotation dataframe
def extract_promoters(cluster_df, annotation_df, genome, upstream=1000):
    """Extract promoter sequences"""
    promoters = {}

    for _, row in cluster_df.iterrows():
        gene_id = row["gene_id"]

        hit = annotation_df[annotation_df['id'].str.strip() == gene_id.strip()]
        if hit.empty:
            logger.info(f"⚠️ Gene {gene_id} not found in annotation.")
            continue

        hit_row = hit.iloc[0]
        chrom = str(hit_row["chromosome"])
        strand = hit_row["strand"]
        start = int(hit_row["start"])
        end = int(hit_row["end"])

        if chrom not in genome:
            logger.info(f"⚠️ Chromosome {chrom} not found in genome fasta.")
            continue

        seq = genome[chrom]
        if strand == "+":
            promoter_start = max(0, start - upstream - 1)
            promoter_end = start - 1
            promoter_seq = seq[promoter_start:promoter_end]
        else:
            promoter_start = end
            promoter_end = min(len(seq), end + upstream)
            promoter_seq = seq[promoter_start:promoter_end].reverse_complement()

        header = f"{gene_id}|{chrom}:{promoter_start+1}-{promoter_end}|{strand}"
        promoters[header] = str(promoter_seq)

    return promoters

# Function to save promoter sequences to a fasta file
def save_promoters(promoters, output_path):
    """Save promoter sequences to fasta"""
    with open(output_path, "w") as f:
        for header, seq in promoters.items():
            f.write(f">{header}\n{seq}\n")
            
# Function to analyze promoter sequences for GC content, motif counts, and length
def analyze_promoters(fasta_file):
    """Analyze promoter sequences: GC content, motif counts, length."""
    results = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()

        gc = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0

        tata = seq.count('TATA')
        caat = seq.count('CAAT')
        at_rich = (seq.count('A') + seq.count('T')) / len(seq) if len(seq) > 0 else 0

        results.append({
            'gene': record.id,
            'length': len(seq),
            'GC_content': round(gc, 4),
            'AT_content': round(at_rich, 4),
            'TATA_count': tata,
            'CAAT_count': caat,
        })
    return pd.DataFrame(results)

# Main function to process the cluster CSV and extract promoter sequences
def main(cluster_csv, output_dir):
    mapping_file = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping.csv"
    annotation_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm_no_chloroplast/annotated_genomes_metabolic"
    fasta_dirs = [
        "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/phytozome_genomes/Phytozome",
        "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/ensembl_genomes/dna",
    ]

    os.makedirs(output_dir, exist_ok=True)

    cluster_df = pd.read_csv(cluster_csv)
    cluster_name = os.path.basename(cluster_csv).replace(".csv", "")
    logger.info(f"🚩 Processing cluster: {cluster_name}")

    mapping_df = pd.read_csv(mapping_file)

    # Match annotation file
    annotation_files = [
        os.path.join(annotation_dir, f) for f in os.listdir(annotation_dir) if f.endswith(".csv")
    ]

    best_file = None
    max_hits = 0
    for ann in annotation_files:
        try:
            ann_df = pd.read_csv(ann)
            if "id" not in ann_df.columns:
                continue
            hits = sum(
                cluster_df["gene_id"].str.contains("|".join(ann_df["id"].astype(str)), na=False)
            )
            if hits > max_hits:
                max_hits = hits
                best_file = ann
        except Exception:
            continue

    if not best_file:
        raise ValueError("❌ Could not match annotation file.")

    logger.info(f"✅ Best matching annotation file: {os.path.basename(best_file)} ({max_hits} gene hits)")

    ann_df = pd.read_csv(best_file)

    dataset_id, organism_name = get_organism_info(best_file, mapping_df)
    fasta_file = find_fasta_file(dataset_id, organism_name, fasta_dirs)
    logger.info(f"✅ Using fasta file: {os.path.basename(fasta_file)}")

    genome = load_genome(fasta_file)
    promoters = extract_promoters(cluster_df, ann_df, genome)

    output_fasta = os.path.join(output_dir, f"{cluster_name}_promoters.fasta")
    save_promoters(promoters, output_fasta)
    logger.info(f"✅ Saved {len(promoters)} promoters to {output_fasta}")
    analysis_df = analyze_promoters(output_fasta)
    output_csv = os.path.join(output_dir, f"{cluster_name}_promoters_analysis.csv")
    analysis_df.to_csv(output_csv, index=False)
    logger.info(f"✅ Saved promoter analysis to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract promoter sequences for cluster genes.")
    parser.add_argument("--cluster_csv", required=True, help="Cluster CSV file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    main(args.cluster_csv, args.output_dir)
