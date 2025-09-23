#!/bin/bash
#SBATCH --job-name=colabfold_gpu_%j
#SBATCH --partition=gpu-itaymay
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/colabfold_gpu_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/colabfold_gpu_%j.err

# Get arguments
FASTA_FILE=$1
OUTPUT_DIR=$2

echo "Starting ColabFold GPU job"
echo "FASTA file: $FASTA_FILE"
echo "Output dir: $OUTPUT_DIR"

# Load conda environment
source ~/.bashrc
conda activate colabfold_gpu

# Create output directory
mkdir -p $OUTPUT_DIR

# Run ColabFold
echo "Running ColabFold batch..."
colabfold_batch \
    --templates \
    --num-recycle 3 \
    --num-models 1 \
    --model-type auto \
    $FASTA_FILE \
    $OUTPUT_DIR

echo "ColabFold GPU job completed"