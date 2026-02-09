#!/usr/bin/env python3
import re
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def main():
    parser = argparse.ArgumentParser(description="Convert tidecluster TAREAN report TSV to FASTA")
    parser.add_argument("-i", "--input", required=True, help="Path to input TSV file")
    parser.add_argument("-o", "--output", default=None, help="Path to output FASTA file")
    args = parser.parse_args()

    df = pd.read_table(args.input)

    records = []

    for row in df.itertuples():
        record = SeqRecord(
            Seq(re.sub(r"<pre>|\r|\n", "", str(row.Consensus))),
            id=str(row.TRC).strip() + "#tidecluster/tarean",
            description=""
        )

        record.annotations["monomer_length"] = row.monomer_length
        record.annotations["kmer"] = row.kmer
        record.annotations["total_score"] = row.total_score
        record.annotations["n_gap50"] = row.n_gap50
        record.annotations["number_of_arrays"] = row.number_of_arrays
        record.annotations["min_array_length"] = row.min_array_length
        record.annotations["max_array_length"] = row.max_array_length
        record.annotations["median_array_length"] = row.median_array_length
        record.annotations["type"] = row.type
        # row.size_of_arrays может быть строкой вроде "min:5011;median:16586;max:1442660"
        size_str = str(row.size_of_arrays).replace(" ", "").replace("<br>", ";")
        size_parts = dict(item.split(":") for item in size_str.split(";") if ":" in item)

        record.annotations["size_of_arrays_min"] = int(size_parts.get("min", 0))
        record.annotations["size_of_arrays_median"] = int(size_parts.get("median", 0))
        record.annotations["size_of_arrays_max"] = int(size_parts.get("max", 0))
        
        record.annotations["total_size"] = row.Total_size

        record.description = ",".join(
            f"{k}={v}" for k, v in record.annotations.items()
        )

        records.append(record)

    if args.output:
        SeqIO.write(records, args.output, "fasta")
    else:
        SeqIO.write(records, sys.stdout, "fasta")

if __name__ == "__main__":
    main()