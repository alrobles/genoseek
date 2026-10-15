#!/usr/bin/env python3
"""Find the best reference genome for the top genominer candidates.

Usage:
    python3 scripts/ingest_reference_genomes.py \
        --candidates data/candidates.tsv \
        --top-n 20 \
        --out data/reference_genome_report.tsv
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import reference_genomes
from src.spa_scheduler import candidate_priority


def load_candidates(tsv_path: str, top_n: int):
    with open(tsv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        r["priority"] = candidate_priority(r)
    rows.sort(key=lambda r: r["priority"], reverse=True)
    return rows[:top_n]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", default="data/candidates.tsv")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--out", default="data/reference_genome_report.tsv")
    parser.add_argument("--sleep", type=float, default=0.5)
    args = parser.parse_args()

    selected = load_candidates(args.candidates, args.top_n)
    print(f"Looking up reference genomes for {len(selected)} candidates")

    fieldnames = [
        "binomial", "dataset", "ref_genome_accession", "ref_genome_name",
        "ref_organism_name", "ref_tax_id", "ref_common_name",
        "ref_assembly_level", "ref_refseq_category", "ref_bioproject",
        "ref_total_sequence_length", "ref_scaffold_n50",
        "ref_number_of_chromosomes", "ref_gc_percent", "ref_genome_coverage",
        "ref_ftp_genomic_fasta", "ref_score", "ref_n_alternative_assemblies",
        "status"
    ]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    # Merge with the previous report so reference lookups stay cumulative
    # across cycles; rows for species re-queried this run are refreshed.
    merged: dict[str, dict] = {}
    out_path = Path(args.out)
    if out_path.is_file():
        with open(out_path, newline="", encoding="utf-8") as f:
            for old in csv.DictReader(f, delimiter="\t"):
                if old.get("binomial"):
                    merged[old["binomial"]] = old
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in selected:
            binomial = row["binomial"]
            dataset = binomial.replace(" ", "_").lower()
            try:
                best = reference_genomes.select_reference(binomial, sleep=args.sleep)
                if best:
                    merged[binomial] = {
                        "binomial": binomial,
                        "dataset": dataset,
                        "ref_genome_accession": best["accession"],
                        "ref_genome_name": best["assembly_name"],
                        "ref_organism_name": best["organism_name"],
                        "ref_tax_id": best.get("tax_id", ""),
                        "ref_common_name": best.get("common_name", ""),
                        "ref_assembly_level": best["assembly_level"],
                        "ref_refseq_category": best["refseq_category"],
                        "ref_bioproject": best.get("bioproject_accession", ""),
                        "ref_total_sequence_length": best["total_sequence_length"],
                        "ref_scaffold_n50": best["scaffold_n50"],
                        "ref_number_of_chromosomes": best["number_of_chromosomes"],
                        "ref_gc_percent": best.get("gc_percent", ""),
                        "ref_genome_coverage": best.get("genome_coverage", ""),
                        "ref_ftp_genomic_fasta": best["ftp_genomic_fasta"],
                        "ref_score": best["score"],
                        "ref_n_alternative_assemblies": best["n_alternative_assemblies"],
                        "status": "ok",
                    }
                    print(f"  {binomial}: {best['accession']} ({best['assembly_level']}, {best['total_sequence_length']} bp)")
                else:
                    merged[binomial] = {"binomial": binomial, "dataset": dataset, "status": "no_reference_found"}
                    print(f"  {binomial}: no reference found")
            except Exception as e:
                merged[binomial] = {"binomial": binomial, "dataset": dataset, "status": f"error: {e}"}
                print(f"  {binomial}: ERROR {e}")

        w.writerows(merged.values())
    print(f"Wrote {args.out} ({len(merged)} species total)")


if __name__ == "__main__":
    main()
