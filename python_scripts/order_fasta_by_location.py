from Bio import SeqIO
import re

def parse_fasta_and_positions(fasta_file):
    gene_positions = []
    sequences = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract information from the description
        description = record.description
        match = re.search(r'chromosome:TAIR10:(\d+):(\d+):(\d+):(-?\d+)', description)
        if match:
            chrom, start, end, strand = match.groups()
            position = (int(start), int(end), int(strand))
            gene_positions.append((record.id, position))
            sequences[record.id] = str(record.seq)

    return gene_positions, sequences

def sort_sequences_by_position(gene_positions, sequences):
    # Sort genes by their start position
    sorted_genes = sorted(gene_positions, key=lambda x: x[1])
    sorted_sequences = {gene_id: sequences[gene_id] for gene_id, _ in sorted_genes}
    return sorted_sequences

def write_ordered_fasta(output_file, sorted_sequences):
    with open(output_file, 'w') as out_file:
        for gene_id, sequence in sorted_sequences.items():
            out_file.write(f">{gene_id}\n")
            out_file.write(f"{sequence}\n")

# Example usage
fasta_file = "/groups/itay_mayrose/alongonda/desktop/arabidopsis/Arabidopsis_thaliana.TAIR10.pep.all.fa"
output_file = "/groups/itay_mayrose/alongonda/desktop/arabidopsis/Arabidopsis_thaliana.TAIR10.pep.ordered.fa"

# Parse the file and get gene positions and sequences
gene_positions, sequences = parse_fasta_and_positions(fasta_file)

# Sort sequences by their positions
sorted_sequences = sort_sequences_by_position(gene_positions, sequences)

# Write the ordered sequences to a new FASTA file
write_ordered_fasta(output_file, sorted_sequences)

print(f"Ordered FASTA file written to {output_file}")
