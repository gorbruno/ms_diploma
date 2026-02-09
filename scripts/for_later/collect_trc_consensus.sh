#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <root_dir> [output.fasta]"
    exit 1
fi

ROOT_DIR="$1"
OUT="${2:-all_TRC_consensus.fasta}"

> "$OUT"

find "$ROOT_DIR" -maxdepth 1 -type d -name 'TRC_*.fasta_tarean' | sort | while read -r dir; do
    consensus="$dir/consensus.fasta"
    [ -f "$consensus" ] || continue

    trc_name=$(basename "$dir" | sed 's/\.fasta_tarean//')

    seqkit head -n 1 "$consensus" \
      | seqkit replace -p '^(.+)$' -r ">$trc_name \$1" \
      >> "$OUT"
done
