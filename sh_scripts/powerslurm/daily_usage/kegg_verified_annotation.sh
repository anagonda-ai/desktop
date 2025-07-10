#!/bin/bash
#SBATCH --job-name=kegg_annot
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=itaym
#SBATCH --time=04:00:00

input_csv=$1
annotated_dir=$2
temp_dir=$3
kegg_query_fasta=$4

filename=$(basename "$input_csv")
output_file="${annotated_dir}/${filename/.csv/_annotated.csv}"
blast_output="${temp_dir}/blast_result_${filename}.xml"
fasta_db="${temp_dir}/${filename/.csv/.fasta}"
id_map="${fasta_db/.fasta/_id_mapping.tsv}"

echo "=== KEGG Reverse Annotation Job ==="
echo "Input CSV: $input_csv"
echo "Output: $output_file"
echo "BLAST DB (genome FASTA): $fasta_db"
echo "KEGG Query: $kegg_query_fasta"
echo "============================"

# Step 1: Generate valid FASTA with short IDs and mapping
echo "Generating FASTA and ID mapping..."
python3 - <<EOF
import pandas as pd

csv_path = "$input_csv"
fasta_path = "$fasta_db"
map_path = "$id_map"

df = pd.read_csv(csv_path)
use_transcript = "transcript_name" in df.columns
assert "id" in df.columns and "sequence" in df.columns

with open(fasta_path, "w") as fasta, open(map_path, "w") as mapping:
    for i, row in df.iterrows():
        gene_id = str(row["id"])
        transcript = str(row["transcript_name"]) if use_transcript else ""
        seq = str(row["sequence"]).replace("*", "")  # remove stop codons if present
        short_id = f"gene_{i+1:05d}"
        long_id = f"{gene_id}|{transcript}" if transcript else gene_id
        fasta.write(f">{short_id}\n{seq}\n")
        mapping.write(f"{short_id}\t{gene_id}\n")
EOF

# Step 2: Create BLAST DB
echo "Creating BLAST DB..."
makeblastdb -in "$fasta_db" -dbtype prot -parse_seqids

# Step 3: Run BLAST
if [ ! -f "$blast_output" ]; then
    echo "Running BLASTP from KEGG → Genome..."
    blastp -task blastp-fast \
        -query "$kegg_query_fasta" \
        -db "$fasta_db" \
        -out "$blast_output" \
        -outfmt 5 \
        -evalue 1e-5 \
        -num_threads 4
else
    echo "BLAST output already exists: $blast_output"
fi

# Step 4: Annotate original CSV
echo "Parsing BLAST output and annotating..."

python3 - <<EOF
import pandas as pd
from Bio.Blast import NCBIXML
from collections import defaultdict

csv_path = "$input_csv"
blast_xml = "$blast_output"
map_path = "$id_map"
output_path = "$output_file"

df = pd.read_csv(csv_path)
id_map = pd.read_csv(map_path, sep="\t", header=None, names=["short", "orig"])
id_dict = dict(zip(id_map["short"], id_map["orig"]))

gene_to_kegg = defaultdict(list)

with open(blast_xml) as handle:
    for record in NCBIXML.parse(handle):
        query_def = record.query
        query_len = record.query_length
        if "$" not in query_def:
            continue
        parts = query_def.split("$")
        if len(parts) < 3:
            continue
        kegg_id, annotation, pathway = map(str.strip, parts[:3])

        best_hit = None
        best_score = -1

        for alignment in record.alignments:
            hit_id = alignment.hit_id.strip().split()[0]
            for hsp in alignment.hsps:
                identity = 100 * hsp.identities / hsp.align_length
                coverage = 100 * hsp.align_length / query_len
                if (
                    hsp.expect < 1e-5 and
                    identity >= 90 and
                    coverage >= 90
                ):
                    if hsp.bits > best_score:
                        best_hit = hit_id
                        best_score = hsp.bits
                break

        if best_hit and best_hit in id_dict:
            orig_id = id_dict[best_hit]
            gene_to_kegg[orig_id].append((kegg_id, annotation, pathway))

kegg_ids, annotations, pathways = [], [], []

for gene_id in df["id"]:
    hits = gene_to_kegg.get(gene_id, [])
    if hits:
        kegg_ids.append(",".join(k[0] for k in hits))
        annotations.append(",".join(k[1] for k in hits))
        pathways.append(",".join(k[2] for k in hits))
    else:
        kegg_ids.append(None)
        annotations.append(None)
        pathways.append(None)

df["kegg_ids"] = kegg_ids
df["annotation"] = annotations
df["pathway"] = pathways
df.to_csv(output_path, index=False)
print(f"✅ Annotation complete: {output_path}")
EOF

echo "=== Job Complete ==="
