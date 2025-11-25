#!/bin/bash
#SBATCH --job-name=tiara_manual
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/tiara_manual_batch01_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/tiara_manual_batch01_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

# Activate conda env
source activate tiara

INPUT_DIR=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125
OUTPUT_DIR=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/tiara


mkdir -p "$OUTPUT_DIR"

# Process all .fa files in the directory
for file in "$INPUT_DIR"/*.fa; do
    if [[ -f "$file" ]]; then
        base_name=$(basename "$file" .fa)
        echo "Processing $base_name ..."
        tiara -m 1000 -i "$file" -o "$OUTPUT_DIR/tiara_${base_name}.txt" -t 8
    fi
done

# Deactivate conda env
conda deactivate
