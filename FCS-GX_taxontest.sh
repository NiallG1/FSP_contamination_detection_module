#!/bin/bash
#SBATCH --job-name=fcsx_gx
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/fcsx_gx_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/fcsx_gx_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=520G
#SBATCH --partition=himem


#input query fasta and taxid
FASTA=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/batch01_191125/EGP017_25_047_best_assembly.fa
OUT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/fcs/taxon_test
TAXID=60108
GXDB=/home/ngarvey/scratch/contamination_detection/FCS/gxdb/gxdb/all.gxi
OUTBASENAME="EGP017_25_047_family_level"

mkdir -p "$OUT"

source activate ncbi_fcsgx

# Match threads to SLURM allocation
export GX_NUM_CORES=$SLURM_CPUS_PER_TASK

run_gx.py --fasta "$FASTA" \
--out-dir "$OUT" \
--gx-db "$GXDB" \
--tax-id "$TAXID" \
--out-basename "$OUTBASENAME"




