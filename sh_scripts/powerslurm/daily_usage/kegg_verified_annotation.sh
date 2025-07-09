#!/bin/bash
#SBATCH --job-name=kegg_annot
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=itaym
#SBATCH --time=04:00:00

# Arguments
input_csv=$1
annotated_dir=$2
temp_dir=$3
kegg_db=$4

filename=$(basename "$input_csv")
output_file="${annotated_dir}/${filename/.csv/_annotated.csv}"
blast_output="${temp_dir}/blast_result_${filename}.xml"
fasta_query="${temp_dir}/${filename/.csv/.fasta}"

echo "=== KEGG Annotation Job ==="
echo "Input CSV: $input_csv"
echo "Output: $output_file"
echo "BLAST DB: $kegg_db"
echo "============================"

# Generate FASTA
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
        -num_threads 4
else
    echo "BLAST output already exists: $blast_output"
fi

# Parse BLAST output and perform optimal bipartite matching
echo "Parsing BLAST output with optimal matching..."

python3 - <<EOF
import pandas as pd
from Bio.Blast import NCBIXML
import networkx as nx
from collections import defaultdict

df = pd.read_csv("$input_csv")
blast_output = "$blast_output"
query_ids = set(df['id'])

# Step 1: collect all valid edges
edges = []
hit_data = {}

with open(blast_output) as handle:
    for record in NCBIXML.parse(handle):
        query_id = record.query
        query_len = record.query_length
        for alignment in record.alignments:
            parts = alignment.hit_def.split('$')
            if len(parts) != 3:
                continue
            subject_id, annotation, pathway = parts
            annotation = annotation.strip()
            pathway = pathway.strip()
            for hsp in alignment.hsps:
                bit_score = hsp.bits
                # Save every hit — no strict filtering
                edges.append((query_id, subject_id, bit_score))
                hit_data[(query_id, subject_id)] = (annotation, pathway, bit_score)
                break  # take top HSP per alignment

# Step 2: Build bipartite graph
G = nx.Graph()
for query_id, subject_id, score in edges:
    G.add_edge(f"Q_{query_id}", f"S_{subject_id}", weight=score)

# Step 3: Max-weight matching
matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)

# Step 4: Extract best match per query gene
best_annotations = {}
best_pathways = {}

for u, v in matching:
    if u.startswith("Q_"):
        query_id, subject_id = u[2:], v[2:]
    else:
        query_id, subject_id = v[2:], u[2:]
    if (query_id, subject_id) in hit_data:
        annotation, pathway, _ = hit_data[(query_id, subject_id)]
        best_annotations[query_id] = annotation
        best_pathways[query_id] = pathway

# Step 5: Add to DataFrame
df['annotation'] = df['id'].map(best_annotations)
df['pathway'] = df['id'].map(best_pathways)
df.to_csv("$output_file", index=False)

print(f"✅ Saved annotated CSV to $output_file")
EOF

echo "=== Job Complete ==="