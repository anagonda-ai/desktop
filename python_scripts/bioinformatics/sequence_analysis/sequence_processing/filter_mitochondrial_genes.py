import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from tqdm import tqdm

# --------- CONFIG ----------
DATASET_DIR = "/groups/itay_mayrose/alongonda/datasets/full_genomes/final_dataset"
MITO_GBK = "/groups/itay_mayrose/alongonda/datasets/full_genomes/human_mitochondrion.gbk"
OUTPUT_DIR = os.path.join(DATASET_DIR, "filtered_no_mito")
EVALUE_THRESHOLD = 1e-5
COVERAGE_THRESHOLD = 70.0  # percent
NUM_THREADS = 16
TEMP_DIR = "temp_blast_mito"
# ---------------------------

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)

# Step 1: Convert GBK → protein FASTA
mito_fasta = os.path.join(TEMP_DIR, "human_mito_proteins.fasta")
if not os.path.exists(mito_fasta):
    mito_records = []
    with open(MITO_GBK) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    seq = feature.qualifiers["translation"][0]
                    gene_id = feature.qualifiers.get("gene", ["unknown"])[0]
                    mito_records.append(SeqRecord(Seq(seq), id=gene_id, description=""))
    SeqIO.write(mito_records, mito_fasta, "fasta")

# Step 2: Create BLAST DB
mito_db = os.path.join(TEMP_DIR, "mito_db")
if not os.path.exists(mito_db + ".pin"):
    subprocess.run(["makeblastdb", "-in", mito_fasta, "-dbtype", "prot", "-out", mito_db], check=True)

# Step 3: BLAST filter
def run_blastp_filter(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    if 'sequence' not in df.columns:
        print(f"⚠️ Skipping {input_csv} - no 'sequence' column")
        return

    # Write protein sequences to FASTA
    query_fasta = os.path.join(TEMP_DIR, "query.fasta")
    records = []
    id_to_idx = {}
    for idx, row in df.iterrows():
        seq = row['sequence'].replace("*", "")
        if seq and set(seq).issubset("ACDEFGHIKLMNPQRSTVWY"):  # check it's protein
            rec_id = f"seq_{idx}"
            records.append(SeqRecord(Seq(seq), id=rec_id, description=""))
            id_to_idx[rec_id] = idx
    SeqIO.write(records, query_fasta, "fasta")

    if not records:
        print(f"⚠️ No protein sequences in {input_csv}")
        df.to_csv(output_csv, index=False)
        return

    # Run BLASTP
    result_xml = os.path.join(TEMP_DIR, "blast_result.xml")
    NcbiblastpCommandline(
        query=query_fasta,
        db=mito_db,
        evalue=EVALUE_THRESHOLD,
        outfmt=5,
        out=result_xml,
        num_threads=NUM_THREADS
    )()

    # Parse BLAST results
    to_remove = set()
    with open(result_xml) as result_handle:
        for record in NCBIXML.parse(result_handle):
            if record.alignments:
                best_hsp = record.alignments[0].hsps[0]
                coverage = 100.0 * best_hsp.align_length / record.query_length
                if coverage >= COVERAGE_THRESHOLD:
                    orig_idx = id_to_idx.get(record.query)
                    if orig_idx is not None:
                        to_remove.add(orig_idx)

    # Remove rows and save
    df_filtered = df.drop(index=to_remove)
    df_filtered.to_csv(output_csv, index=False)

# Step 4: Loop through dataset
csv_files = [f for f in os.listdir(DATASET_DIR) if f.endswith(".csv")]
for file_name in tqdm(csv_files, desc="Filtering mitochondrial genes"):
    input_path = os.path.join(DATASET_DIR, file_name)
    output_path = os.path.join(OUTPUT_DIR, file_name)
    run_blastp_filter(input_path, output_path)

print("✅ All files filtered. Output saved to:", OUTPUT_DIR)
