#!/bin/bash
#SBATCH --job-name=fcsx_gx
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/fcsx_gx_%j.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/fcsx_gx_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=520G
#SBATCH --partition=himem


#input query fasta and taxid
FASTA=/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/EGP017_25_044_best_assembly.fa
OUT=/home/ngarvey/scratch/contamination_detection/manual_pipeline/results/fcs
TAXID=1354922
GXDB=/home/ngarvey/scratch/contamination_detection/FCS/gxdb/gxdb/all.gxi

mkdir -p "$OUT"

source activate ncbi_fcsgx

# Match threads to SLURM allocation
export GX_NUM_CORES=$SLURM_CPUS_PER_TASK

run_gx.py --fasta "$FASTA" \
--out-dir "$OUT" \
--gx-db "$GXDB" \
--tax-id "$TAXID"


# -----------------------------
# Postprocess GX output report
# -----------------------------
REPORT="$OUT/$(basename "$FASTA" .fasta).FCS-GX.taxonomy.rpt"
CLEAN_REPORT="$OUT/$(basename "$FASTA" .fasta).FCS-GX.taxonomy.clean.tsv"

# Skip first line, remove leading # from header, convert | to tabs, collapse extra tabs
tail -n +2 "$REPORT" \
    | sed '1s/^#//' \
    | sed 's/|/\t/g' \
    | tr -s '\t' \
    > "$CLEAN_REPORT"

echo "Cleaned FCS-GX report saved to: $CLEAN_REPORT"
