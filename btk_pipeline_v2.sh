#!/bin/bash
#SBATCH --job-name=btk_generate
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_v2_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_v2_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_003_best_assembly.fa 
BAM=/home/ngarvey/scratch/contamination_detection/manual_pipeline/bam/EGP017_25_003_best_assembly_sorted.bam
YAML=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_025_003.yaml
RESULT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools
BUSCO=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_003_busco.tsv
TAXONOMY=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/comparison/EGP017_25_003_blobtools_taxonomy.tsv

source /mnt/apps/users/ngarvey/conda/etc/profile.d/mamba.sh

mamba activate btk

# ================================
# Step 1: Create BlobToolKit dataset
# ================================

cd "$RESULT"


blobtools create \
    --fasta "$SAMPLE" \
    --meta "$YAML" \
    ./EGP017_25_003_blobdir



#================================
# Step 3: Add BUSCO completeness data
# ================================

blobtools add \
    --busco "$BUSCO" \
    ./EGP017_25_003_blobdir



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
    ./EGP017_25_003_blobdir


#add coverage data
blobtools add \
    --cov "$BAM" \
    --threads 24 \
    ./EGP017_25_003_blobdir


