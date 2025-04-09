import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import pandas as pd

# Email and API key
Entrez.email = "alongonda@mail.tau.ac.il"
Entrez.api_key = "13b51c448a0ba246af8a501b77f7dc0fe309"

# Batch and speed settings
NCBI_DELAY = 0.1
BATCH_SIZE = 100
MAX_WORKERS = 8  # concurrency

def batch_fetch_gene_coordinates(gene_ids):
    """Batch fetch RefSeq coordinates for multiple GeneIDs."""
    gene_info = {}

    id_string = ",".join(gene_ids)

    try:
        handle = Entrez.esummary(db="gene", id=id_string, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        for docsum in record['DocumentSummarySet']['DocumentSummary']:
            real_gene_id = docsum.attributes.get('uid')
            genomic_info = docsum.get('GenomicInfo')
            if genomic_info:
                chrom_acc = genomic_info[0].get('ChrAccVer')
                chrom_num = genomic_info[0].get('ChrLoc')
                start = int(genomic_info[0].get('ChrStart')) + 1
                end = int(genomic_info[0].get('ChrStop')) + 1
                if start > end:
                    start, end = end, start
                gene_info[real_gene_id] = (chrom_num, chrom_acc, start, end)
            else:
                gene_info[real_gene_id] = None

    except Exception as e:
        print(f"⚠️ Batch fetch failed: {e}")

    return gene_info

def batch_resolve_locus_tags(locus_tags):
    """Batch resolve locus tags (gene names) to GeneIDs using OR search."""
    resolved = {}

    if not locus_tags:
        return resolved

    terms = [f"{tag}[Gene Name]" for tag in locus_tags]
    query = " OR ".join(terms)

    try:
        handle = Entrez.esearch(db="gene", term=query, retmax=5000, retmode="xml")
        search_results = Entrez.read(handle)
        handle.close()

        id_list = search_results.get("IdList", [])
        
        if id_list:
            id_string = ",".join(id_list)
            handle = Entrez.esummary(db="gene", id=id_string, retmode="xml")
            summaries = Entrez.read(handle)
            handle.close()

            for docsum in summaries['DocumentSummarySet']['DocumentSummary']:
                real_gene_id = docsum.attributes.get('uid')
                name = docsum.get('Name')
                if real_gene_id and name:
                    resolved[name] = real_gene_id

    except Exception as e:
        print(f"⚠️ Batch locus tag search failed: {e}")

    return resolved

def process_fasta_file(input_fasta, output_csv):
    """Process one FASTA file in batch + concurrent mode."""
    rows = []
    gene_records = list(SeqIO.parse(input_fasta, "fasta"))

    # Step 1: Collect IDs and locus tags separately
    id_map = {}
    gene_ids = []
    locus_tags = []

    for record in gene_records:
        header = record.id
        organism_code, gene_id = header.split(":")
        if gene_id.isdigit():
            id_map[header] = gene_id
            gene_ids.append(gene_id)
        else:
            id_map[header] = gene_id
            locus_tags.append(gene_id)

    # Step 2: Resolve locus tags (batch)
    resolved_locus_tags = batch_resolve_locus_tags(locus_tags)

    # Merge resolved locus tags
    for header, tag in id_map.items():
        if not tag.isdigit():
            real_id = resolved_locus_tags.get(tag)
            if real_id:
                id_map[header] = real_id
            else:
                id_map[header] = None  # unresolved

    # Step 3: Collect real GeneIDs
    all_gene_ids = [gid for gid in id_map.values() if gid is not None]

    # Step 4: Batch fetch gene coordinates concurrently
    batches = [all_gene_ids[i:i+BATCH_SIZE] for i in range(0, len(all_gene_ids), BATCH_SIZE)]
    gene_coordinates = {}

    def fetch_batch(batch):
        return batch_fetch_gene_coordinates(batch)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [executor.submit(fetch_batch, batch) for batch in batches]
        for future in as_completed(futures):
            try:
                batch_result = future.result()
                gene_coordinates.update(batch_result)
            except Exception as e:
                print(f"⚠️ Error fetching a batch: {e}")

    # Step 5: Update headers
    updated_count = 0
    missing_count = 0

    for record in gene_records:
        header = record.id
        organism_code, original_id = header.split(":")
        real_id = id_map.get(header)

        if real_id and real_id in gene_coordinates and gene_coordinates[real_id]:
            chrom_acc, chrom_num, start, end = gene_coordinates[real_id]
            row = {
                "Organism": organism_code,
                "GeneID": original_id,
                "Chromosome": chrom_num,
                "Chromosome_Accession": chrom_acc,
                "Start": start,
                "End": end,
                "Sequence": str(record.seq)
            }
            print(f"✅ Updated: {header}")
            updated_count += 1
        else:
            row = {
                "Organism": organism_code,
                "GeneID": original_id,
                "Chromosome": None,
                "Chromosome_Accession": None,
                "Start": None,
                "End": None,
                "Sequence": str(record.seq)
            }
            print(f"⚠️ Missing data for: {header}")
            missing_count += 1

        rows.append(row)

    # Step 5: Save CSV
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)

    # Summary
    print(f"\nSummary for {input_fasta}:")
    print(f"  Total: {len(gene_records)}")
    print(f"  Updated: {updated_count}")
    print(f"  Missing: {missing_count}\n")

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
                output_fasta = os.path.join(organism_output_path, fasta_file.replace(".fasta", "_updated.csv"))

                print(f"  ⚡ Processing {fasta_file}")
                process_fasta_file(input_fasta, output_fasta)
                print(f"  ✅ Done {fasta_file}")

# Example usage
input_root = "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta"
output_root = "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta_updated"

process_all_kegg_fasta(input_root, output_root)
