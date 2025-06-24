import argparse
import os
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

PROMOTER_UPSTREAM = 1000
PROMOTER_DOWNSTREAM = 100


def parse_annotation_csv(annotation_csv):
    gene_coords = {}
    with open(annotation_csv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_id = row['id']
            gene_coords[gene_id] = {
                'chrom': row['chromosome'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': int(row['strand'])
            }
    return gene_coords


def load_genome_fasta(genome_fasta):
    return {rec.id: rec.seq for rec in SeqIO.parse(genome_fasta, "fasta")}


def extract_promoter(seq, start, end, strand):
    if strand == 1:
        promoter_start = max(0, start - PROMOTER_UPSTREAM)
        promoter_end = end + PROMOTER_DOWNSTREAM
        promoter_seq = seq[promoter_start:promoter_end]
    else:
        promoter_start = max(0, start - PROMOTER_DOWNSTREAM)
        promoter_end = end + PROMOTER_UPSTREAM
        promoter_seq = seq[promoter_start:promoter_end].reverse_complement()
    return promoter_seq


def parse_cluster_fasta(cluster_fasta):
    gene_entries = []
    for record in SeqIO.parse(cluster_fasta, "fasta"):
        header_parts = record.description.split("|")
        cluster_id = header_parts[0].strip()
        gene_id = header_parts[1].strip()
        genome_id = header_parts[3].strip()
        annotation_file = header_parts[3].strip()
        gene_entries.append((cluster_id, gene_id, genome_id, annotation_file))
    return gene_entries


def run_meme(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        "meme", input_fasta,
        "-oc", output_dir,
        "-dna",
        "-mod", "zoops",
        "-nmotifs", "5",
        "-minw", "6",
        "-maxw", "15"
    ]
    subprocess.run(cmd, check=True)


def main(cluster_fasta, genome_fasta, annotation_csv, output_dir):
    gene_coords = parse_annotation_csv(annotation_csv)
    genome_seqs = load_genome_fasta(genome_fasta)
    gene_entries = parse_cluster_fasta(cluster_fasta)

    cluster_promoters = {}
    for cluster_id, gene_id, genome_id, _ in gene_entries:
        if gene_id not in gene_coords:
            continue
        coords = gene_coords[gene_id]
        chrom = coords['chrom']
        if chrom not in genome_seqs:
            continue
        promoter_seq = extract_promoter(
            genome_seqs[chrom], coords['start'], coords['end'], coords['strand']
        )
        if cluster_id not in cluster_promoters:
            cluster_promoters[cluster_id] = []
        cluster_promoters[cluster_id].append(SeqRecord(
            promoter_seq,
            id=f"{gene_id}|{cluster_id}|{genome_id}",
            description=""
        ))

    for cluster_id, records in cluster_promoters.items():
        cluster_dir = Path(output_dir) / cluster_id
        cluster_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = cluster_dir / "promoters.fasta"
        SeqIO.write(records, fasta_path, "fasta")
        run_meme(str(fasta_path), str(cluster_dir / "meme"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract promoter sequences and run MEME")
    parser.add_argument("--cluster_fasta", required=True, help="FASTA file of protein sequences from MGC cluster")
    parser.add_argument("--genome_fasta", required=True, help="Genomic DNA FASTA file")
    parser.add_argument("--annotation_csv", required=True, help="Gene annotation CSV file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    main(
        cluster_fasta=args.cluster_fasta,
        genome_fasta=args.genome_fasta,
        annotation_csv=args.annotation_csv,
        output_dir=args.output_dir
    )
