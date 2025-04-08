import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import os

# NCBI email
Entrez.email = "alongonda@mail.tau.ac.il"
Entrez.api_key = "13b51c448a0ba246af8a501b77f7dc0fe309"

# Limit NCBI hits (officially 3/sec for anonymous users)
NCBI_DELAY = 0.1  # 3 queries/sec

# How many parallel threads?
MAX_WORKERS = 8

def fetch_gene_coordinates(gene_id):
    """Fetch both RefSeq and GenBank coordinates for a GeneID."""
    try:
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        docsum = record['DocumentSummarySet']['DocumentSummary'][0]
        genomic_info = docsum.get('GenomicInfo')

        refseq_data = None
        if genomic_info:
            chrom = genomic_info[0].get('ChrAccVer')
            start = int(genomic_info[0].get('ChrStart')) + 1
            end = int(genomic_info[0].get('ChrStop')) + 1
            if start > end:
                start, end = end, start
            refseq_data = (chrom, start, end)

        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        records = Entrez.read(handle)
        handle.close()

        return refseq_data

    except Exception as e:
        print(f"⚠️ Primary fetch failed for {gene_id}: {e}")
        print(f"🔄 Trying fallback search for locus tag...")
        
        # Now treat it as a Locus Tag
        try:
            search_handle = Entrez.esearch(db="gene", term=f"{gene_id}[Gene Name]", retmode="xml")
            search_results = Entrez.read(search_handle)
            search_handle.close()

            id_list = search_results.get("IdList", [])
            if id_list:
                real_gene_id = id_list[0]
                print(f"✅ Found real GeneID for {gene_id}: {real_gene_id}")
                return fetch_gene_coordinates(real_gene_id)  # recursive call with real ID
            else:
                print(f"❌ No match found for locus tag {gene_id}")
                return None, None, None

        except Exception as e2:
            print(f"❌ Failed fallback search for {gene_id}: {e2}")
            return None, None, None

def process_fasta_file(input_fasta, output_fasta):
    """Process one FASTA file concurrently."""
    updated_records = []

    gene_records = list(SeqIO.parse(input_fasta, "fasta"))

    def process_record(record):
        header = record.id
        organism_code, gene_id = header.split(":")
        refseq_info = fetch_gene_coordinates(gene_id)
        if refseq_info:
            chrom, gb_start, gb_end = refseq_info
            new_header = f"{header}|{chrom}|{gb_start}|{gb_end}"
        else:
            new_header = header

        updated_record = SeqRecord(
            Seq(str(record.seq)),
            id=new_header,
            description=""
        )
        time.sleep(NCBI_DELAY)  # respect API rate limit
        return updated_record

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(process_record, record) for record in gene_records]
        for future in as_completed(futures):
            try:
                updated_records.append(future.result())
            except Exception as e:
                print(f"Error processing a record: {e}")

    # Save updated FASTA
    with open(output_fasta, "w") as out_f:
        SeqIO.write(updated_records, out_f, "fasta-2line")

def process_all_kegg_fasta(root_input_folder, root_output_folder):
    """Process ALL .fasta files inside root_input_folder and save to root_output_folder."""
    for organism_dir in os.listdir(root_input_folder):
        organism_input_path = os.path.join(root_input_folder, organism_dir)

        if os.path.isdir(organism_input_path):
            print(f"\n📂 Processing organism: {organism_dir}")
            organism_output_path = os.path.join(root_output_folder, organism_dir)
            os.makedirs(organism_output_path, exist_ok=True)

            fasta_files = [f for f in os.listdir(organism_input_path) if f.endswith(".fasta")]

            for fasta_file in fasta_files:
                input_fasta = os.path.join(organism_input_path, fasta_file)
                output_fasta = os.path.join(organism_output_path, fasta_file)

                print(f"  ⚡ Processing {fasta_file}")

                process_fasta_file(input_fasta, output_fasta)

                print(f"  ✅ Done {fasta_file}")

# Example usage
input_root = "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta"
output_root = "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta_updated"

process_all_kegg_fasta(input_root, output_root)
