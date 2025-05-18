#!/bin/bash
#SBATCH --job-name=blast_with_params
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

# Get the input FASTA file from the command line argument
QUERY_FASTA=$1

# Define the BLAST database
DB=$2

# Define the output directory
BLAST_OUTPUT=$3



blastp -query "$QUERY_FASTA" -db "$DB" -out "$BLAST_OUTPUT" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 0.001 -num_threads $(nproc) -max_target_seqs 5