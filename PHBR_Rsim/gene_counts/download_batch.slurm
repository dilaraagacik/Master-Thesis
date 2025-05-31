#!/bin/bash
#SBATCH --job-name=gdc_r_dl
#SBATCH --output=logs/download_%A_%a.out
#SBATCH --error=logs/download_%A_%a.err
#SBATCH --array=0-115        
#SBATCH --time=02:00:00
#SBATCH --mem=2G

# Load conda and activate your environment
source /home/labgroups/ccgg/miniforge3/etc/profile.d/conda.sh
conda activate /home/dilara/projects/perimeta/perimeta

UUID_FILES=(uuid_batches/uuid_batch_*.txt)
UUID_FILE=${UUID_FILES[$SLURM_ARRAY_TASK_ID]}

Rscript download_batch.R "$UUID_FILE" "downloads_batch_$SLURM_ARRAY_TASK_ID"
