#!/bin/bash
#SBATCH --job-name=bam_test
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/bam_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/bam_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_056_best_assembly.fa
R1=/mnt/shared/projects/rbgk/projects/FSP/00_RawData/01_SeqData/06_2025.07.09_EdGen_192_samples/EGP017_25_056/EGP017_25_056_S76_L001_R1.fastq.gz
R2=/mnt/shared/projects/rbgk/projects/FSP/00_RawData/01_SeqData/06_2025.07.09_EdGen_192_samples/EGP017_25_056/EGP017_25_056_S76_L001_R2.fastq.gz
RESULT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/bam

#blobtools requires a .bam.csi insread of .bam.bsi (reason unkown) so this script indexes bam files with the -c parameter

#Align genome back to the reads to check coverage:

source /mnt/apps/users/ngarvey/conda/etc/profile.d/conda.sh
conda activate bwa-mem2

bwa-mem2 index "$SAMPLE"

bwa-mem2 mem -t 16 "$SAMPLE" "$R1" "$R2" > alignment.sam

conda activate samtools

samtools view -bS alignment.sam -@ 16 > alignment.bam

samtools sort -@ 16 -o "$RESULT/$(basename "$SAMPLE" .fa)_sorted.bam" alignment.bam

samtools index -c "$RESULT/$(basename "$SAMPLE" .fa)_sorted.bam"

# Optional: clean up intermediate files
rm -f alignment.sam alignment.bam

conda deactivate



