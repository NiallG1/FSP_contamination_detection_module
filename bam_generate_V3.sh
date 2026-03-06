#!/bin/bash
#SBATCH --job-name=bam_test
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/bam_v3_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/bam_v3_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/10_003_comparison/EGP017_25_B5_001/best_assembly_info_and_QC/samtools_pypolca/EGP017_25_B5_001_best_assembly_pypolca_sorted.bam
#R1=/mnt/shared/projects/rbgk/projects/FSP/00_RawData/01_SeqData/06_2025.07.09_EdGen_192_samples/EGP017_25_056/EGP017_25_056_R1.fastq.gz
#R2=/mnt/shared/projects/rbgk/projects/FSP/00_RawData/01_SeqData/06_2025.07.09_EdGen_192_samples/EGP017_25_056/EGP017_25_056_R2.fastq.gz
RESULT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/10_003_comparison/EGP017_25_B5_001/best_assembly_info_and_QC/samtools_pypolca

#blobtools requires a .bam.csi insread of .bam.bsi (reason unkown) so this script indexes bam files with the -c parameter

#Align genome back to the reads to check coverage:

source /mnt/apps/users/ngarvey/conda/etc/profile.d/conda.sh
conda activate bam_generation

samtools index -c "$RESULT/$(basename "$SAMPLE" .fa)_sorted_v3.bam"

conda deactivate




