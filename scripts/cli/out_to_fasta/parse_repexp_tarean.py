#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
from pathlib import Path
import re
import argparse
import sys


CSV_HEADER_SIZE = 5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate TAREAN consensus FASTA with cluster metadata"
    )

    parser.add_argument(
        "-i", "--input-dir",
        type=Path,
        required=True,
        help="Directory with TAREAN*.fasta files"
    )

    parser.add_argument(
        "-s", "--sample-name",
        type=str,
        required=True,
        default="sample",
        help="Sample name that would be placed in all sequence names"
    )

    parser.add_argument(
        "-c", "--cluster-table",
        type=Path,
        required=True,
        help="Path to CLUSTER_TABLE.csv"
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=False,
        help="Output FASTA file (default: all_tarean_consensus.fasta)"
    )

    return parser.parse_args()


def main(args):
    records = []
    sample_name = args.sample_name.lower().strip().replace(" ", "_")

    # --- Parse FASTA ---
    for fasta_path in args.input_dir.glob("TAREAN*.fasta"):
        for record in SeqIO.parse(fasta_path, "fasta"):
            match_rank = re.search(
                r"tarean_consensus_rank_(\d+)",
                fasta_path.name,
                flags=re.IGNORECASE
            )
            tarean_consensus_rank = (
                int(match_rank.group(1)) if match_rank else "undefined"
            )

            match_id = re.search(
                r"cl(\d+)_.*_(\d+)nt$",
                record.id,
                flags=re.IGNORECASE
            )
            cluster, length = (
                int(match_id.group(1)),
                int(match_id.group(2))
            ) if match_id else ("undefined", "undefined")

            record.annotations["cluster"] = cluster
            record.annotations["length"] = length
            record.annotations["tarean_consensus_rank"] = tarean_consensus_rank

            records.append(record)

    if not records:
        raise RuntimeError("No TAREAN FASTA records found")

    # --- Read cluster table ---
    df = pd.read_table(
        args.cluster_table,
        header=CSV_HEADER_SIZE,
    )

    cluster_to_keep = {rec.annotations["cluster"] for rec in records}

    df = df[df["Cluster"].isin(cluster_to_keep)]

    cluster_meta = (
        df.set_index("Cluster")
          .to_dict(orient="index")
    )

    # --- Read summary table ---
    table_meta = (
        pd.read_table(
            args.cluster_table,
            skiprows=lambda x: x not in range(6),
            names=["statement", "value"]
        )
        .set_index("statement")["value"]
        .to_dict()
    )

    # --- Annotate records ---
    for rec in records:
        cluster = rec.annotations.get("cluster")
        meta = {**cluster_meta.get(cluster, {}), **table_meta}

        annotation = meta.get("TAREAN_annotation", "").lower()

        match annotation:
            case s if "high" in s:
                suffix = "/tarean_high_conf"
            case s if "low" in s:
                suffix = "/tarean_low_conf"
            case _:
                suffix = "/tarean_unkn_conf"

        rec.id += ("#repeatexplorer_" + sample_name + suffix)

        for k, v in meta.items():
            rec.annotations[k] = v

        rec.description = ",".join(
            f"{k}={v}" for k, v in rec.annotations.items()
        )
    
    # --- Sort records by cluster ---
    records = sorted(
        records,
        key=lambda r: r.annotations.get("cluster", float("inf"))
    )
    # --- Write output ---
    if args.output:
        SeqIO.write(records, args.output, "fasta")
    else:
        SeqIO.write(records, sys.stdout, "fasta")

if __name__ == "__main__":
    args = parse_args()
    main(args)