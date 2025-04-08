from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time

# Set your email address to identify your queries to NCBI
Entrez.email = "alongonda@mail.tau.ac.il"

def fetch_gene_coordinates(gene_id):
    """Fetch RefSeq and GenBank coordinates for a given NCBI GeneID."""
    try:
        # Step 1: Fetch gene summary
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()      
        
        docsum = record['DocumentSummarySet']['DocumentSummary'][0]
        genomic_info = docsum.get('GenomicInfo')
        # Fetch from RefSeq (NC_)
        refseq_data = None
        if genomic_info:
            chrom = genomic_info[0].get('ChrAccVer')    # RefSeq accession (e.g., NC_...)
            start = int(genomic_info[0].get('ChrStart')) + 1 # Convert 0-based to 1-based
            end = int(genomic_info[0].get('ChrStop')) + 1
            strand = genomic_info[0].get('Strand')
            if start > end:
                start, end = end, start
            refseq_data = (chrom, start, end, strand)
            
        # Step 2: Fetch full gene record (for Related Sequences GenBank)
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        records = Entrez.read(handle)
        handle.close()

        genbank_data = None
        related_sequences = None

        try:
            comments = records[0].get('Entrezgene_comments', [])
            for comment in comments:
                if comment.get('Gene-commentary_heading') == 'Related Sequences':
                    products = comment.get('Gene-commentary_products', [])
                    for prod in products:
                        if prod.get('Gene-commentary_heading') == 'Genomic':
                            acc = prod.get('Gene-commentary_accession')
                            version = prod.get('Gene-commentary_version')
                            full_acc = f"{acc}.{version}"

                            # Now get the start/end
                            seq_info = prod.get('Gene-commentary_seqs', [])
                            if seq_info:
                                interval = seq_info[0]['Seq-loc_int']['Seq-interval']
                                start_genbank = int(interval['Seq-interval_from']) + 1
                                end_genbank = int(interval['Seq-interval_to']) + 1
                                strand_info = interval['Seq-interval_strand']['Na-strand']
                                strand_value = strand_info.attributes.get('value', 'plus')
                                strand = 1 if strand_value == 'plus' else -1

                                if start_genbank > end_genbank:
                                    start_genbank, end_genbank = end_genbank, start_genbank

                                genbank_data = (full_acc, start_genbank, end_genbank, strand)
        except Exception as e:
            print(f"Error parsing related sequences: {e}")
            genbank_data = None

        return refseq_data, genbank_data

    except Exception as e:
        print(f"Error fetching data for GeneID {gene_id}: {e}")
        return None, None

def process_fasta(fasta_file, output_file):
    """Process a FASTA file to extract gene information and fetch coordinates."""
    updated_records = []
    
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            header = record.id
            organism_code, gene_id = header.split(":")
            
            refseq_info, genbank_info = fetch_gene_coordinates(gene_id)
            
            print(f"\nGene {gene_id} ({organism_code}):")
            
            if refseq_info:
                chrom, start, end, strand = refseq_info
                strand_str = "+" if strand == 1 else "-"
                print(f"  RefSeq: {chrom}:{start}-{end} ({strand_str} strand)")
            else:
                print("  RefSeq: Not available")

            if genbank_info:
                genbank_seq, gb_start, gb_end, gb_strand = genbank_info
                gb_strand_str = "+" if gb_strand == 1 else "-"
                print(f"  GenBank: {genbank_seq}:{gb_start}-{gb_end} ({gb_strand_str} strand)")
                new_header = f"{header}|{genbank_seq}|{gb_start}|{gb_end}"
            else:
                new_header = header
                print("  GenBank: Not available")
                
            # Create a new SeqRecord with updated header
            updated_record = SeqRecord(
                Seq(str(record.seq)),
                id=new_header,
                description=""
            )
            updated_records.append(updated_record)

            print(f"Updated {header} → {new_header}")
            # Be nice to NCBI servers
            time.sleep(0.35)
        
        # Write to new FASTA
    with open(output_file, "w") as out_f:
        SeqIO.write(updated_records, out_f, "fasta-2line")

# Example usage
process_fasta(
    "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta/aans/aans00071.fasta",
    "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta/aans/aans00071_UPDATED.fasta"
)
