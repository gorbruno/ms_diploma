#!/bin/bash
#SBATCH -p workq
#SBATCH -w proxima2
#SBATCH --job-name=repeatmasker
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=logs/repeatmasker_%j.log
#SBATCH --hint=nomultithread

set -euo pipefail

############################
# CONFIG
############################
THREADS=${SLURM_CPUS_PER_TASK:-20}

# Проверяем наличие RepeatMasker в PATH
if ! command -v RepeatMasker &> /dev/null
then
    echo "RepeatMasker не найден в PATH. Добавьте его в PATH перед запуском."
    exit 1
fi

############################
# ARGS CHECK
############################
if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch $0 <repeatmasker_flags>"
    echo "Example: sbatch $0 -pa 20 -lib Saguaro-families.fa -dir results -s -a -gff Cgig_v2_SGP5p_21.fa"
    exit 1
fi

FLAGS="$@"

############################
# RUN REPEATMASKER
############################
echo "----------------------------------------"
echo "Running RepeatMasker with the following command:"
echo "RepeatMasker $FLAGS"
echo "Threads allocated by SLURM: $THREADS"
echo "Start time: $(date)"
echo

RepeatMasker $FLAGS

echo
echo "Done!"
echo "End time: $(date)"