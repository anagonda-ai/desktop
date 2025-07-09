#!/bin/bash
#SBATCH --job-name=kegg_annot
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=itaym
#SBATCH --time=04:00:00

# Arguments from Python
input_csv=$1
annotated_dir=$2
temp_dir=$3
kegg_db=$4

# File names
filename=$(basename "$input_csv")
output_file="${annotated_dir}/${filename/.csv/_annotated.csv}"
blast_output="${temp_dir}/blast_result_${filename}.xml"
fasta_query="${temp_dir}/${filename/.csv/.fasta}"

echo "=== KEGG Annotation Job ==="
echo "Input CSV: $input_csv"
echo "Output: $output_file"
echo "BLAST DB: $kegg_db"
echo "============================"

# Generate FASTA file from CSV
echo "Generating FASTA..."
awk -F',' 'BEGIN {OFS="\n"} NR==1 {for (i=1; i<=NF; i++) if ($i=="id") id=i; else if ($i=="sequence") seq=i} NR>1 {print ">"$id, $seq}' "$input_csv" > "$fasta_query"

# Run BLASTP
if [ ! -f "$blast_output" ]; then
    echo "Running BLASTP..."
    blastp -task blastp-fast \
        -query "$fasta_query" \
        -db "$kegg_db" \
        -out "$blast_output" \
        -outfmt 5 \
        -evalue 1e-3 \
        -num_threads 4 \
        -max_target_seqs 5
else
    echo "BLAST output already exists: $blast_output"
fi

# Parse BLAST and generate output CSV
echo "Parsing BLAST output..."

python3 - <<EOF
import os
import pandas as pd
from Bio.Blast import NCBIXML

input_csv = "$input_csv"
blast_output = "$blast_output"
output_file = "$output_file"

df = pd.read_csv(input_csv)
id_to_pathway = {}
id_to_annotation = {}

with open(blast_output) as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        query_id = record.query
        query_length = record.query_length
        best_score = 0
        best_pathway = None
        best_annotation = None

        for alignment in record.alignments:
            header = alignment.hit_def
            if '$' in header:
                _, annotation, pathway = header.split('$')
                annotation = annotation.strip()
                pathway = pathway.strip()
            else:
                continue

            for hsp in alignment.hsps:
                coverage = (hsp.align_length / query_length) * 100
                identity = (hsp.identities / hsp.align_length) * 100
                bit_score = hsp.bits

                if coverage >= 90.0 and identity >= 90.0 and bit_score > best_score:
                    best_score = bit_score
                    best_pathway = pathway
                    best_annotation = annotation

        if best_pathway:
            id_to_pathway[query_id] = best_pathway
            id_to_annotation[query_id] = best_annotation

df['pathway'] = df['id'].map(id_to_pathway)
df['annotation'] = df['id'].map(id_to_annotation)
df.to_csv(output_file, index=False)
print(f"✅ Saved annotated CSV to {output_file}")
EOF

echo "=== Job Complete ==="
