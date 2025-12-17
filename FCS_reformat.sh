#!/bin/bash
#SBATCH --job-name=fcsx_gx
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/fcsx_gx_reformat_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/fcsx_gx_reformat_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=10g
#SBATCH --partition=short

RPT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/fcs/synthetic/EGP017_25_047_Com2.taxonomy.rpt
OUT_DIR=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/fcs/synthetic

mkdir -p "$OUT_DIR"

cut --complement -f 5,12,18,24,30 "$RPT" | tail -n +2 | sed '1s/^#//' > "$OUT_DIR"/EGP017_25_Com_2_FCS_GX.tsv


