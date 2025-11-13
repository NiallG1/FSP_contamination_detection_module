#!/bin/bash
#SBATCH --job-name=contam_preprocess
#SBATCH --output=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/contam_preprocess_%A_%a.out
#SBATCH --error=/home/ngarvey/scratch/contamination_detection/manual_pipeline/error_out/contam_preprocess_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --partition=short
#SBATCH --array=1-5     # <-- adjust based on number of FASTA files

# ---- Activate seqkit environment ----
source activate seqkit

# ---- Define input and output directories ----
INPUT_DIR="/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/genomes_for_decontamination_Niall"
OUTDIR="/home/ngarvey/scratch/contamination_detection/manual_pipeline/genomes/processed"
mkdir -p "$OUTDIR"

# ---- Determine which FASTA file to process ----
input_fasta=$(ls "$INPUT_DIR"/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")

# ---- Validate input ----
if [ -z "$input_fasta" ]; then
    echo "❌ Error: No input FASTA found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

# ---- Get base name ----
base_name=$(basename "$input_fasta" .fa)

# ---- Define output file paths ----
above_1kbp="${OUTDIR}/${base_name}_gt1kbp.fasta"
below_1kbp="${OUTDIR}/${base_name}_lt1kbp.fasta"

# ---- Filter and split ----
echo "🚀 Processing: $input_fasta"

seqkit seq -g -m 1000 "$input_fasta" > "$above_1kbp"
seqkit seq -g -M 999 "$input_fasta" > "$below_1kbp"

# ---- Report ----
echo "✅ Done for $base_name"
echo "📂 Output files:"
echo "   ≥1kbp:     $above_1kbp"
echo "   250–999bp: $below_1kbp"

conda deactivate
