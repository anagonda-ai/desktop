#!/bin/bash
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --partition=itaym

query_fasta="$1"
target_list_file="$2"
output_dir="$3"
blast_db_dir="$4"

mkdir -p "$output_dir"

while IFS= read -r target_fasta; do
    target_base=$(basename "$target_fasta" .fasta)
    db_path="$blast_db_dir/$target_base"

    if [ ! -f "${db_path}.pin" ]; then
        makeblastdb -in "$target_fasta" -dbtype prot -out "$db_path"
    fi

    out_file="$output_dir/$(basename "$query_fasta")_vs_$(basename "$target_fasta").txt"

    if [ ! -f "$out_file" ]; then
        blastp -query "$query_fasta" -db "$db_path" -evalue 0.001 -outfmt 6 -out "$out_file"
    fi
done < "$target_list_file"
