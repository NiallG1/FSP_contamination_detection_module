#!/bin/bash
#SBATCH --job-name=tiara_fcs_R
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/tiara_fcs_R_10_003_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/tiara_fcs_R_10_003_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

# Activate conda env
#change source to your own profile source on CDC
source /mnt/apps/users/ngarvey/conda/etc/profile.d/conda.sh
conda activate R

#Run R script
Rscript /home/ngarvey/scratch/contamination_detection/manual_pipeline/scripts/tiara_vs_fcs_comparison_v4.R
