
class ChloroplastFilterProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: chloroplast_filter.py."""
    
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
import glob
import subprocess
from pathlib import Path
from multiprocessing import Pool
from Bio import SeqIO
from argparse import ArgumentParser

def make_blast_db(fasta_file):
    """Make BLAST database from the chloroplast genes."""
    db_prefix = fasta_file.rsplit(".fasta", 1)[0]
    if not Path(f"{db_prefix}.pin").exists():
        subprocess.run(["makeblastdb", "-in", fasta_file, "-dbtype", "prot", "-out", db_prefix], check=True)
    return db_prefix

def run_blastp(query_fasta, db, evalue=1e-5):
    """Run BLASTP and return matched IDs."""
    output = subprocess.run([
        "blastp",
        "-query", query_fasta,
        "-db", db,
        "-evalue", str(evalue),
        "-outfmt", "6 qseqid"
    ], capture_output=True, text=True, check=True)
    
    matched_ids = set(output.stdout.strip().splitlines())
    return matched_ids

def filter_fasta(input_fasta, exclude_ids, output_fasta):
    """Write FASTA entries not in exclude_ids to output_fasta."""
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id not in exclude_ids:
                SeqIO.write(record, out_f, "fasta")

def process_organism(args):
    pep_fasta, output_dir, chloroplast_db, evalue = args
    org_name = Path(pep_fasta).parent.name
    filtered_output = Path(output_dir) / f"{org_name}.filtered.fasta"
    try:
        hits = run_blastp(pep_fasta, chloroplast_db, evalue)
        filter_fasta(pep_fasta, hits, filtered_output)
        return f"[✓] Processed {org_name} ({len(hits)} matches filtered)"
    except subprocess.CalledProcessError as e:
        return f"[✗] Error in {org_name}: {e}"

def main():
    parser = ArgumentParser(description="Filter chloroplast-like proteins from plant genomes.")
    parser.add_argument("--chloroplast_dir", required=True, help="Path to translated chloroplast genes (FASTA)")
    parser.add_argument("--organisms_dir", required=True, help="Path to directory with organism subdirectories")
    parser.add_argument("--output_dir", required=True, help="Where to save filtered FASTA files")
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--evalue", type=float, default=1e-5)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Merge chloroplast FASTAs into one
    combined_fasta = Path(args.output_dir) / "combined_chloroplast_genes.fasta"
    with open(combined_fasta, "w") as out_f:
        for fasta in glob.glob(f"{args.chloroplast_dir}/translated_*.fasta"):
            out_f.write(Path(fasta).read_text())

    # Create BLAST DB
    db_prefix = make_blast_db(str(combined_fasta))

    # Find all peptide FASTA files
    pep_fastas = glob.glob(f"{args.organisms_dir}/*/*.fasta")

    # Run filtering in parallel
    jobs = [(pep, args.output_dir, db_prefix, args.evalue) for pep in pep_fastas]
    with Pool(processes=args.threads) as pool:
        for result in pool.imap_unordered(process_organism, jobs):
            logger.info(result)

if __name__ == "__main__":
    main()
