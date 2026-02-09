#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Clusterize and annotate FASTA sequences based on MMseqs cluster output.")
    parser.add_argument("-c", "--cluster_tsv", required=True, help="Path to MMseqs cluster TSV file")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file with sequences")
    parser.add_argument("-o", "--output", required=False, help="Path to write updated FASTA")
    return parser.parse_args()

def parse_annotation(seq_id):
    if "#" in seq_id:
        name = seq_id.split("#", 1)[0]
        annot = seq_id.split("#", 1)[1]
        annot_type = annot.split("/")[0]
        if "tidecluster" in annot_type.lower():
            source = "tidecluster"
        elif "repeatexplorer" in annot_type.lower():
            source = "repeatexplorer"
        elif "rnd" in annot_type.lower():
            source = "repeatmodeler"
        else:
            source = "unknown"
        return {"name": name, "annot": annot, "source": source, "full_name": seq_id}
    return None

def main():
    args = parse_args()

    # читаем кластеры
    clusters = defaultdict(list)
    with open(args.cluster_tsv) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            clusters[rep].append(member)

    # создаём аннотации для реплик
    rep_annotations = {}
    for rep, members in clusters.items():
        meta = defaultdict(list)
        a_rep = parse_annotation(rep)
        meta["self"] = a_rep

        for m in members:
            a = parse_annotation(m)
            if a:
                meta["other"].append(a)
        rep_annotations[rep] = meta

    # проверяем наличие TideCluster/RepeatExplorer и меняем full_name
    for rep, meta in rep_annotations.items():
        self_annot = meta["self"]
        has_tarean = any(o["source"] in ("tidecluster", "repeatexplorer") for o in meta["other"])
        if has_tarean:
            tarean_annot = "Satellite"
            self_annot["full_name"] = f"{self_annot['name']}#{tarean_annot}"

    # обновляем FASTA
    records = list(SeqIO.parse(args.fasta, "fasta"))

    for record in records:
        if record.id in rep_annotations:
            meta = rep_annotations[record.id]

            # флаг кластеризации
            clusterized = True
            if len(meta.get("other")) == 1 and meta.get("self").get("full_name") == meta.get("other")[0].get("full_name"):
                clusterized = False
            clusterized_flag = str(clusterized).lower()

            self_annot = meta.get("self")
            other = meta.get("other")
            record.id = self_annot["full_name"]

            # объединяем мета-информацию
            other_contigs = ";".join(o["full_name"] for o in other) if other else ""
            meta_merged = f"merged_contigs={other_contigs},"
            meta_info = meta_merged + f"clusterized={clusterized_flag}"

            # обновляем описание
            desc_parts = record.description.split(maxsplit=1)
            if len(desc_parts) > 1:
                old_desc = desc_parts[1]
                record.description = old_desc.strip() + "," + meta_info
            else:
                record.description = meta_info

    # сохраняем обновлённый FASTA
    if args.output:
        SeqIO.write(records, args.output, "fasta")
    else:
        SeqIO.write(records, sys.stdout, "fasta")
        

if __name__ == "__main__":
    main()