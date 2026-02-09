#!/usr/bin/env python3

import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Join mmseqs2 alignment with taxonomy report.")
    parser.add_argument("-a", "--aln_file", type=Path, required=True,
                        help="Path to the mmseqs2 alignment file")
    parser.add_argument("-r", "--report_file", type=Path, required=True,
                        help="Path to the mmseqss2 taxonomy report file")
    parser.add_argument("-o", "--out_file", type=Path, required=True,
                        help="Path to the output TSV file")

    args = parser.parse_args()

    aln_file = args.aln_file
    report_file = args.report_file
    out_file = args.out_file

    # ---------- read taxonomy report ----------
    tax = {}
    with report_file.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            accession = fields[0]
            tax[accession] = fields[1:]

    # ---------- write joined file with header ----------
    header_aln = [
        "query_id", "accession", "pident", "aln_len", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bitscore"
    ]

    header_tax = [
        "n_hits", "score", "sum_score", "coverage", "taxid", "rank", "taxon_name", "lineage"
    ]

    with aln_file.open() as fin, out_file.open("w") as fout:
        # write header
        fout.write("\t".join(header_aln + header_tax) + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            accession = fields[1]

            if accession in tax:
                fout.write(line + "\t" + "\t".join(tax[accession]) + "\n")
            else:
                fout.write(line + "\t" + "\t".join(["NA"] * len(header_tax)) + "\n")

    print(f"[OK] written: {out_file}")

if __name__ == "__main__":
    main()
