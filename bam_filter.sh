#!/bin/bash
#SBATCH --job-name=bam_filter
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/bam_filter__%A_%a.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/bam_filter_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --partition=short

set -euo pipefail

# Input arguments
BAM_DIR=/home/ngarvey/scratch/contamination_detection/manual_pipeline/bam
MINLEN=250

# Load samtools 
source activate samtools

echo "Starting BAM filtering job on $(date)"
echo "Directory: $BAM_DIR"
echo "Minimum contig length: $MINLEN bp"
echo

for bam in "$BAM_DIR"/*.bam; do
    [ -e "$bam" ] || continue
    echo "Processing: $bam"

    outbam="${bam%.bam}_filtered.bam"

    samtools idxstats "$bam" | awk -v m=$MINLEN '$2 >= m {print $1}' > keep_contigs.txt

    if [[ -s keep_contigs.txt ]]; then
        samtools view -h "$bam" -L <(awk '{print $1 "\t0\t999999999"}' keep_contigs.txt) -o "$outbam"
        echo "  → Filtered BAM written: $outbam"
    else
        echo "  ⚠️ No contigs ≥ $MINLEN bp found; skipping $bam."
    fi

    echo
done
