#!/bin/bash
#SBATCH --job-name=btk_filter
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_filter_url_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_filter_url_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_003_best_assembly.fa
#FILTER=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044_blobdir.current.json
OUTPUT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools
BLOBDIR=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools/EGP017_25_003_blobdir


#Align genome back to the reads to check coverage:

source /mnt/apps/users/ngarvey/conda/etc/profile.d/mamba.sh
mamba activate btk



blobtools filter \
    --query-string "http://localhost:8080/view/all/dataset/EGP017_25_003_blobdir/blob?plotShape=circle&plotGraphics=svg&taxonomy--Keys=25%2C14%2C17%2C32%2C13%2C10%2C1%2C7%2C0#Filters" \
    --fasta "$SAMPLE" \
    "$BLOBDIR"


#blobtools filter \
#    --query-string "http://localhost:8080/view/all/dataset/Myblobdir/blob?gc--Max=0.500&hifi.aln.sorted_cov--Active=true&hifi.aln.sorted_cov--Max=5180#Filters" \
#    --fasta assembly/original_assembly.fa \
#    ./Myblobdir



blobtools filter \
    --query-string "http://localhost:8080/view/all/dataset/EGP017_25_003_blobdir/blob?plotShape=circle&plotGraphics=svg&taxonomy--Keys=25%2C14%2C17%2C32%2C13%2C10%2C1%2C7%2C0#Filters" \
    --fasta /home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_003_best_assembly.fa \
    /home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools/EGP017_25_003_blobdir

blobtools filter \
     --json /home/ngarvey/scratch/contamination_detection/manual_pipeline/blobtools/EGP017_25_003_blobdir.current.json \
     --fasta /home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_003_best_assembly.fa \
     /home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools/EGP017_25_003_blobdir
    
