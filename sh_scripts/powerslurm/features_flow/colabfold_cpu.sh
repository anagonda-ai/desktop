#!/bin/bash
#SBATCH --job-name=colabfold_cpu_%j
#SBATCH --partition=gpu-itaymay
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/colabfold_cpu_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/colabfold_cpu_%j.err

# Get arguments
FASTA_FILE=$1
OUTPUT_DIR=$2

echo "Starting ColabFold CPU job"
echo "FASTA file: $FASTA_FILE"
echo "Output dir: $OUTPUT_DIR"

# Load conda environment
source ~/.bashrc
conda activate colabfold_gpu

# Create output directory
mkdir -p $OUTPUT_DIR

# Run ColabFold (CPU mode)
echo "Running ColabFold batch (CPU)..."
colabfold_batch \
    --templates \
    --num-recycle 3 \
    --num-models 1 \
    --model-type auto \
    $FASTA_FILE \
    $OUTPUT_DIR

echo "ColabFold CPU job completed"