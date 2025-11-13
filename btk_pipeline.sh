#!/bin/bash
#SBATCH --job-name=btk_generate
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/btk_pipeline_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G


SAMPLE=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/genomes_for_decontamination_Niall/EGP017_25_044.fa
BAM=/mnt/shared/projects/rbgk/projects/FSP/03_Output/02_Genome_Assemblies/02_Final_genome_assemblies/06_EG/EGP017_25_044/best_assembly_info_and_QC/bwa_mem2_samtools_pilon/EGP017_25_044_best_assembly_pilon_sorted.bam
YAML=/home/ngarvey/scratch/contamination_detection/manual_pipeline/EG_025_044.yml
RESULT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/blobtools
BUSCO=/mnt/shared/projects/rbgk/projects/FSP/03_Output/02_Genome_Assemblies/02_Final_genome_assemblies/06_EG/EGP017_25_044/best_assembly_info_and_QC/busco_specific_pilon/full_table.tsv
TAXONOMY=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044__blobtools_taxonomy.tsv

source activate btk

# ================================
# Step 1: Create BlobToolKit dataset
# ================================
BLOB_DIR="${RESULT}/${SAMPLE_NAME}_blob_dir"

echo "Creating BlobToolKit dataset for ${SAMPLE_NAME} ..."
blobtools create \
    --fasta "$SAMPLE" \
    --meta "$YAML" \
    --out "$BLOB_DIR"

echo "✅ BlobToolKit dataset created at: $BLOB_DIR"


# ================================
# Step 2: Add coverage information (BAM)
# ================================
echo "Adding coverage information from BAM file..."
blobtools add \
    --cov "$BAM" \
    --threads 24 \
    ./"$BLOB_DIR"

# ================================
# Step 3: Add BUSCO completeness data
# ================================
echo "Adding BUSCO completeness results..."
blobtools add \
    --busco "$BUSCO" \
    "$BLOB_DIR"

echo "✅ BUSCO results added."


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
    "$BLOB_DIR"

