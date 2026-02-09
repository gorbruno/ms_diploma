#!/usr/bin/env python3

from pathlib import Path
import re
import argparse
from Bio import SeqIO
import sys


PATTERN = (
    r"\(\s*"
    r".*?RepeatScout\s+Family\s+Size\s*=\s*(\d+)\s*,"
    r".*?Final\s+Multiple\s+Alignment\s+Size\s*=\s*(\d+)\s*,"
    r".*?Localized\s+to\s*(\d+)\s+out\s+of\s*(\d+)"
    r".*?\)"
)

def annotate_record(record, pattern):
    match = re.search(pattern, record.description, flags=re.IGNORECASE)

    if match:
        record.annotations["repeatscout_family_size"] = int(match.group(1))
        record.annotations["final_mult_aln_size"] = int(match.group(2))
        record.annotations["contigs_localized"] = int(match.group(3))
        record.annotations["contigs_total"] = int(match.group(4))
    else:
        record.annotations["repeatscout_family_size"] = "unknown"
        record.annotations["final_mult_aln_size"] = "unknown"
        record.annotations["contigs_localized"] = "unknown"
        record.annotations["contigs_total"] = "unknown"

    record.description = ",".join(
        f"{k}={v}" for k, v in record.annotations.items()
    )

    return record


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate RepeatModeler FASTA headers with parsed RepeatScout statistics"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        type=Path,
        help="Input FASTA file from RepeatModeler"
    )

    parser.add_argument(
        "-o", "--output",
        required=False,
        type=Path,
        help="Output FASTA with normalized descriptions"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    records = []

    for record in SeqIO.parse(args.input, "fasta"):
        records.append(annotate_record(record, PATTERN))

    if args.output:
        SeqIO.write(records, args.output, "fasta")
    else:
        SeqIO.write(records, sys.stdout, "fasta")

if __name__ == "__main__":
    main()