#!/bin/bash
#SBATCH --job-name=alon_interactive
#SBATCH --partition=gpu-itaymay
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1              # ✅ Request 1 GPU
#SBATCH --mem=32G
#SBATCH --time=24:00:00  # Adjust as needed

srun --partition=gpu-general --ntasks=1 --cpus-per-task=8 --gres=gpu:1 --mem=16G --time=04:00:00 --pty bash

