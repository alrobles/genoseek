#!/usr/bin/env python3
"""Ingest SRA/BioSample metadata for the top genominer candidates and build popmaps.

Usage:
    python3 scripts/ingest_genetic_candidates.py \
        --candidates data/candidates.tsv \
        --top-n 10 \
        --out data/genetic_ingest_report.tsv \
        --genetic-dir data/genetic
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import sra_popmap
from src.spa_scheduler import candidate_priority


def load_candidates(tsv_path: str, top_n: int, min_sra: int = 5,
                    prefer_missing_dir: str | None = None):
    with open(tsv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        r["priority"] = candidate_priority(r)
    # penalize species with very low spatial/genetic scores
    rows.sort(key=lambda r: r["priority"], reverse=True)
    rows = [r for r in rows if int(r.get("ncbi_sra_count") or 0) >= min_sra]
    if prefer_missing_dir:
        gdir = Path(prefer_missing_dir)

        def _has_popmap(row) -> bool:
            ds = row["binomial"].replace(" ", "_").lower()
            d = gdir / ds
            return d.is_dir() and any(d.glob("*_popmap.tsv"))

        missing = [r for r in rows if not _has_popmap(r)]
        present = [r for r in rows if _has_popmap(r)]
        rows = missing + present
        print(f"{len(missing)} candidates have no popmap yet (prioritized)")
    return rows[:top_n]


def main():
    parser = argparse.ArgumentParser(description="Build popmaps from SRA metadata for top candidates.")
    parser.add_argument("--candidates", default="data/candidates.tsv")
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--min-sra", type=int, default=5)
    parser.add_argument("--prefer-missing", action="store_true",
                        help="Prioritize species with no existing popmap in --genetic-dir")
    parser.add_argument("--out", default="data/genetic_ingest_report.tsv")
    parser.add_argument("--genetic-dir", default="data/genetic")
    parser.add_argument("--max-runs", type=int, default=1000)
    parser.add_argument("--coord-res-deg", type=float, default=0.5)
    parser.add_argument("--sleep", type=float, default=0.4)
    args = parser.parse_args()

    selected = load_candidates(
        args.candidates, args.top_n, args.min_sra,
        prefer_missing_dir=args.genetic_dir if args.prefer_missing else None,
    )
    print(f"Selected {len(selected)} candidates for genetic metadata ingestion")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.genetic_dir).mkdir(parents=True, exist_ok=True)

    report_rows = []
    for row in selected:
        binomial = row["binomial"]
        dataset = binomial.replace(" ", "_").lower()
        out_dir = Path(args.genetic_dir) / dataset
        try:
            summary = sra_popmap.build_popmap(
                binomial=binomial,
                out_dir=str(out_dir),
                max_runs=args.max_runs,
                coord_res_deg=args.coord_res_deg,
                sleep=args.sleep,
            )
            report_rows.append({
                "binomial": binomial,
                "dataset": dataset,
                "n_sra_runs": summary["n_sra_runs"],
                "n_biosamples": summary["n_biosamples"],
                "n_with_coords": summary["n_with_coords"],
                "n_populations": summary["n_populations"],
                "popmap": summary["popmap"],
                "status": "ok",
            })
            print(f"  {binomial}: {summary['n_sra_runs']} runs, {summary['n_with_coords']} with coords, {summary['n_populations']} pops")
        except Exception as e:
            report_rows.append({
                "binomial": binomial,
                "dataset": dataset,
                "n_sra_runs": 0,
                "n_biosamples": 0,
                "n_with_coords": 0,
                "n_populations": 0,
                "popmap": "",
                "status": f"error: {e}",
            })
            print(f"  {binomial}: ERROR {e}")

    fieldnames = [
        "binomial", "dataset", "n_sra_runs", "n_biosamples", "n_with_coords",
        "n_populations", "popmap", "status"
    ]
    # Merge with the previous report so cycles stay cumulative: species
    # ingested in earlier runs keep their rows unless re-processed now.
    merged: dict[str, dict] = {}
    out_path = Path(args.out)
    if out_path.is_file():
        with open(out_path, newline="", encoding="utf-8") as f:
            for old in csv.DictReader(f, delimiter="\t"):
                if old.get("binomial"):
                    merged[old["binomial"]] = old
    for r in report_rows:
        merged[r["binomial"]] = r

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(merged.values())
    print(f"Wrote report: {args.out} ({len(merged)} species total, {len(report_rows)} this run)")


if __name__ == "__main__":
    main()
