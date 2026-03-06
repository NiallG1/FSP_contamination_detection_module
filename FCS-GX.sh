#!/bin/bash
#SBATCH --job-name=fcsx_gx
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/Com_5_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/Com_5_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=520G
#SBATCH --partition=himem


#input query fasta and taxid
FASTA=/mnt/shared/projects/rbgk/projects/FSP/03_Output/01_QC/03_Decontamination/02_synthetic_genomes/Com_5.fa
OUT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/fcs/synthetic
TAXID=5061
GXDB=/home/ngarvey/scratch/contamination_detection/FCS/gxdb/gxdb/all.gxi
OUTBASENAME="EGP017_Com_5.fa"

mkdir -p "$OUT"

source /mnt/apps/users/ngarvey/conda/etc/profile.d/conda.sh

conda activate ncbi_fcsgx

# Match threads to SLURM allocation
export GX_NUM_CORES=$SLURM_CPUS_PER_TASK

run_gx.py --fasta "$FASTA" \
--out-dir "$OUT" \
--gx-db "$GXDB" \
--tax-id "$TAXID" \
--out-basename "$OUTBASENAME"





