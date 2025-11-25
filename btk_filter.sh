#!/bin/bash
#SBATCH --job-name=btk_filter
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_filter_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_filter_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044_best_assembly.fa
FILTER=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044_blobdir.current.json
OUTPUT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools/EGP017_25_044_blobdir_filtered
BLOBDIR=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools/EGP017_25_044_blobdir


#Align genome back to the reads to check coverage:

source /mnt/apps/users/ngarvey/conda/etc/profile.d/mamba.sh
mamba activate btk


blobtools filter --json "$FILTER" \
     --fasta "$SAMPLE" \
     "$BLOBDIR" \
     --output "$OUTPUT"
