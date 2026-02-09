#!/bin/bash
#SBATCH -p workq
#SBATCH -w proxima2
#SBATCH --job-name=mmseqs-tax
#SBATCH --cpus-per-task=20
#SBATCH --gpus=2
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=logs/mmseqs_tax_%j.log
#SBATCH --hint=nomultithread

set -euo pipefail

############################
# CONFIG
############################
THREADS=${SLURM_CPUS_PER_TASK:-20}

TAXDB=/home/agorbunov/shared/databases/mmseqs-ictv-nr/ictv_nr_db_pad

TMPBASE=${SLURM_TMPDIR:-$PWD/tmp_mmseqs}
mkdir -p "$TMPBASE"

############################
# ARGS CHECK
############################
if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch $0 <fasta1> [fasta2 ... fastaN]"
    exit 1
fi

echo "Using tmp dir: $TMPBASE"
echo "Threads: $THREADS"
echo "Start time: $(date)"
echo

############################
# MAIN LOOP
############################
for FASTA in "$@"; do
    if [ ! -f "$FASTA" ]; then
        echo "ERROR: file not found: $FASTA"
        exit 1
    fi

    BASENAME=$(basename "$FASTA")
    NAME=${BASENAME%.*}

    OUTDIR="${NAME}/vir_tax_${NAME}"
    TMPDIR="$TMPBASE/$NAME"

    echo "----------------------------------------"
    echo "Processing: $FASTA"
    echo "Output dir: $OUTDIR"
    echo "Tmp dir: $TMPDIR"

    mmseqs easy-taxonomy \
        "$FASTA" \
        "$TAXDB" \
        "$OUTDIR" \
        "$TMPDIR" \
        --threads "$THREADS" \
        --gpu 1 \
        --blacklist "" \
        --tax-lineage 1
done

echo
echo "Done!"
echo "End time: $(date)"