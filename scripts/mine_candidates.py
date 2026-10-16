#!/usr/bin/env python3
"""Mine PubMed and GBIF Literature for species with population-genetic data.

Usage:
    python3 scripts/mine_candidates.py --seeds config/seeds.csv \
        --groups config/target_groups.json --out data/mine.db \
        --tsv data/candidates.tsv --limit-per-query 100 \
        --enrich-gbif --enrich-ncbi --enrich-places --top-n 50
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import query_builder
from src import search
from src import entities
from src import export
from src import geoparse
from src import gbif_check
from src import ncbi_check


def _process_rows(rows, source, allowed_genera, seed_geo, enrich_places):
    records = []
    for row in rows:
        text = f"{row.get('title', '')} {row.get('abstract', '')}"
        markers = entities.detect_markers(text)
        if not markers:
            continue
        binomials = entities.extract_binomials(text, allowed_genera)
        if not binomials:
            continue
        score = entities.evidence_score(text, markers)
        places = ",".join(sorted(geoparse.extract_places(text))) if enrich_places else ""
        if source == "pubmed":
            record = {
                "source": "pubmed",
                "external_id": row.get("pmid", ""),
                "doi": row.get("doi", ""),
                "title": row.get("title", ""),
                "abstract": row.get("abstract", ""),
                "authors": row.get("authors", ""),
                "pub_year": row.get("pub_year", ""),
                "journal": row.get("journal", ""),
                "keywords": row.get("keywords", ""),
            }
        else:
            record = {
                "source": "gbif_literature",
                "external_id": row.get("gbif_id", ""),
                "doi": row.get("doi", ""),
                "title": row.get("title", ""),
                "abstract": row.get("abstract", ""),
                "authors": row.get("authors", ""),
                "pub_year": row.get("year", ""),
                "journal": row.get("source", ""),
                "keywords": ",".join(row.get("keywords") or []),
            }
        record.update({
            "binomials": binomials,
            "markers": markers,
            "geography": seed_geo,
            "places": places,
            "score": score,
        })
        records.append(record)
    return records


def main():
    parser = argparse.ArgumentParser(description="Mine literature for population-genetic candidates.")
    parser.add_argument("--seeds", default="config/seeds.csv")
    parser.add_argument("--groups", default="config/target_groups.json")
    parser.add_argument("--out", default="data/mine.db")
    parser.add_argument("--tsv", default="data/candidates.tsv")
    parser.add_argument("--limit-per-query", type=int, default=100)
    parser.add_argument("--year-min", type=int, default=1990)
    parser.add_argument("--year-max", type=int, default=None)
    parser.add_argument("--sources", default="pubmed,gbif", help="Comma-separated sources to query")
    parser.add_argument("--enrich-gbif", action="store_true", help="Query GBIF occurrence counts/extent for top candidates")
    parser.add_argument("--enrich-ncbi", action="store_true", help="Query NCBI SRA/BioSample/BioProject counts for top candidates")
    parser.add_argument("--enrich-places", action="store_true", help="Extract geographic place names from abstracts")
    parser.add_argument("--top-n", type=int, default=50, help="Number of top candidates to enrich with external APIs")
    args = parser.parse_args()

    seeds = query_builder.load_seeds(args.seeds)
    groups = query_builder.load_target_groups(args.groups)
    allowed_genera = entities.load_genera(args.groups)

    sources = set(args.sources.split(","))
    records = []

    n_queries = 0
    for seed in seeds:
        queries = query_builder.build_queries(seed, groups, max_per_seed=60)
        print(f"[{seed.get('seed_id','?')}] {len(queries)} queries", flush=True)
        for q, query_geo in queries:
            n_queries += 1
            if n_queries % 50 == 0:
                print(f"  ...processed {n_queries} queries, {len(records)} records", flush=True)
            if "pubmed" in sources:
                rows = search.search_pubmed(
                    q, limit=args.limit_per_query,
                    year_min=args.year_min, year_max=args.year_max,
                )
                records.extend(_process_rows(rows, "pubmed", allowed_genera, query_geo, args.enrich_places))
            if "gbif" in sources:
                rows = search.search_gbif_literature(
                    q, limit=args.limit_per_query,
                    year_min=args.year_min, year_max=args.year_max,
                )
                records.extend(_process_rows(rows, "gbif_literature", allowed_genera, query_geo, args.enrich_places))

    print(f"Collected {len(records)} paper×marker records with candidate species.", flush=True)

    # Enrich top candidates with GBIF / NCBI
    enrichment: dict[str, dict] = {}
    if records and (args.enrich_gbif or args.enrich_ncbi):
        # Rank candidates by summed score
        cand_score: dict[str, float] = defaultdict(float)
        for r in records:
            for b in r["binomials"]:
                cand_score[b] += r["score"]
        top = sorted(cand_score.items(), key=lambda x: x[1], reverse=True)[:args.top_n]
        print(f"Enriching top {len(top)} candidates via GBIF/NCBI...", flush=True)
        for i, (binomial, _) in enumerate(top, 1):
            if i % 10 == 0:
                print(f"  ...enriched {i}/{len(top)} candidates", flush=True)
            extra: dict = {}
            if args.enrich_gbif:
                extra.update(gbif_check.check(binomial))
            if args.enrich_ncbi:
                extra.update(ncbi_check.check(binomial))
            enrichment[binomial] = extra

    export.save(args.out, args.tsv, records, enrichment=enrichment)
    print(f"Wrote {args.out} and {args.tsv}.", flush=True)


if __name__ == "__main__":
    main()
