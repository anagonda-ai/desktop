#!/bin/bash
#SBATCH --job-name=blast_organism_on_fasta_and_normalize
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

# Check if a query file was provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta>"
    exit 1
fi

# Get the input FASTA file from the command line argument
INPUT_FASTA=$1

# Define the BLAST database
DB=$2

# Define the output directory
OUT_DIR=$(basename $INPUT_FASTA)

# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Define the output file
OUT="$OUT_DIR/blast_output.txt"


blastp -query "$INPUT_FASTA" -db "$DB" -out "$OUT" -outfmt 6 -evalue 1e-5 -num_threads $(nproc)

# Check if the BLAST run was successful
if [ $? -eq 0 ]; then
    echo "BLAST search completed successfully. Results are saved in $OUT"

    CSV_OUT="$OUT_DIR/blast_output.csv"
    # Convert the OUT to .csv
    python3 /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/sequence_analysis/blast_utils/blast_to_csv.py "$OUT" "$CSV_OUT"
    
    # Run the normalize_blast_results.py script with the BLAST output as input
    python3 /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/sequence_analysis/blast_utils/normalizer.py "$CSV_OUT"
    
    # Check if the normalization script ran successfully
    if [ $? -eq 0 ]; then
        echo "Normalization script completed successfully."
        
        # Now the normalized output file is generated with the name 'normalized_blast_scores.csv'
        NORMALIZED_FILE="$OUT_DIR/normalized_blast_scores.csv"
        
        # Run the score_filter.py script with the normalized results
        python3 /groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/sequence_analysis/blast_utils/score_filter.py "$NORMALIZED_FILE"

        # Remove the intermediate files
        rm "$$OUT" "$$CSV_OUT" "$$NORMALIZED_FILE"
        
        # Check if the filtering script ran successfully
        if [ $? -eq 0 ]; then
            echo "Filtering script completed successfully."
        else
            echo "Filtering script failed."
        fi
    else
        echo "Normalization script failed."
    fi
else
    echo "BLAST search failed."
fi
