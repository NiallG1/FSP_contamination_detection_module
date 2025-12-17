#!/bin/bash
#SBATCH --job-name=btk_generate
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_v2_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_v2_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G


SAMPLE=/mnt/shared/projects/rbgk/projects/FSP/03_Output/01_QC/03_Decontamination/02_synthetic_genomes/Com_1.fa 
BAM=/home/ngarvey/scratch/contamination_detection/manual_pipeline/bam/EGP017_25_056_best_assembly_sorted.bam
YAML=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_025_Com1.yaml
RESULT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools
BUSCO=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_056_busco.tsv
TAXONOMY=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/comparison/EGP017_Com1_blobtools_taxonomy.tsv



#you have to use mamba for this instead of conda, as conda has many issues running btk.

source /mnt/apps/users/ngarvey/conda/etc/profile.d/mamba.sh

mamba activate btk

# ================================
# Step 1: Create BlobToolKit dataset
# ================================

#move into the directory where you want to store your blobdirectory.
#you have to change the name to match your desired name.


cd "$RESULT"


blobtools create \
    --fasta "$SAMPLE" \
    --meta "$YAML" \
    ./EGP017_25_Com1_blobdir



#================================
# Step 2: Add BUSCO completeness data
# ================================

blobtools add \
    --busco "$BUSCO" \
    ./EGP017_25_Com1_blobdir



# ================================
# Step 3: Add taxonomy assignments
# ================================

#add taxonomy assignments from tiara and fcs
blobtools add \
    --text "$TAXONOMY" \
    --text-delimiter "\t" \
    --text-cols "seq_id=identifiers,taxonomy=taxonomy" \
    --text-header \
    --key plot.cat=taxonomy \
    ./EGP017_25_Com1_blobdir


# ================================
# Step 4: Add coverage
# ================================

#add coverage data
blobtools add \
    --cov "$BAM" \
    --threads 24 \
    ./EGP017_25_Com1_blobdir


