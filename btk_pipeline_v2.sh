#!/bin/bash
#SBATCH --job-name=btk_generate
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_v2_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_v2_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044_best_assembly.fa 
BAM=/home/ngarvey/scratch/contamination_detection/manual_pipeline/bam/EGP017_25_044_best_assembly.bam
YAML=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_025_044.yaml
RESULT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools
BUSCO=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044_busco_table.tsv
TAXONOMY=/mnt/shared/projects/rbgk/projects/FSP/03_Output/01_QC/03_Decontamination/manual_pipeline/results/comparison/graphs/EGP017_25_044_blobtools_taxonomy_underscore.tsv

source /mnt/apps/users/ngarvey/conda/etc/profile.d/mamba.sh

mamba activate btk

# ================================
# Step 1: Create BlobToolKit dataset
# ================================

cd "$RESULT"


blobtools create \
    --fasta "$SAMPLE" \
    --meta "$YAML" \
    ./EGP017_25_044_blobdir



#================================
# Step 3: Add BUSCO completeness data
# ================================

blobtools add \
    --busco "$BUSCO" \
    ./EGP017_25_044_blobdir



# ================================
# Step 4: Add taxonomy assignments
# ================================

#add taxonomy assignments from tiara and fcs
blobtools add \
    --text "$TAXONOMY" \
    --text-delimiter "\t" \
    --text-cols "seq_id=identifiers,taxonomy=taxonomy" \
    --text-header \
    --key plot.cat=taxonomy \
    ./EGP017_25_044_blobdir


#add coverage data
blobtools add \
    --cov "$BAM" \
    --threads 24 \
    ./EGP017_25_044_blobdir


