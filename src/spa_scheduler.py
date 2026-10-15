"""Prepare a SPA run queue from genominer candidates and submit Slurm jobs."""
from __future__ import annotations
import csv
import json
import math
import os
import shutil
from pathlib import Path
from collections import defaultdict

from src import build_occurrences


def load_groups(groups_path: str) -> dict:
    with open(groups_path, encoding="utf-8") as f:
        return json.load(f)


def group_for_genus(genus: str, groups: dict) -> str:
    for gk, g in groups.items():
        if genus in g.get("genera", []):
            return g.get("group", gk.title())
    return ""


def candidate_priority(row: dict) -> float:
    """Higher = better for SPA pilot scheduling."""
    score = float(row.get("score", 0) or 0)
    occ = int(row.get("gbif_occurrence_count") or 0)
    sra = int(row.get("ncbi_sra_count") or 0)
    # Reward score, log spatial coverage, log genetic evidence
    return score + 2.0 * math.log10(occ + 1) + 3.0 * math.log10(sra + 1)


def build_queue(
    candidates_tsv: str,
    groups_path: str,
    out_csv: str,
    min_occurrences: int = 5000,
    min_sra: int = 5,
    top_n: int = 20,
):
    """Write a prioritized SPA queue CSV."""
    groups = load_groups(groups_path)
    with open(candidates_tsv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    filtered = []
    for r in rows:
        occ = int(r.get("gbif_occurrence_count") or 0)
        sra = int(r.get("ncbi_sra_count") or 0)
        if occ >= min_occurrences and sra >= min_sra:
            r["priority"] = candidate_priority(r)
            r["group"] = group_for_genus(r.get("genus", ""), groups)
            filtered.append(r)

    filtered.sort(key=lambda x: x["priority"], reverse=True)
    top = filtered[:top_n]

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dataset", "canonical_name", "group", "gbif_taxon_key",
        "gbif_occurrence_count", "ncbi_sra_count", "ncbi_biosample_count",
        "ncbi_bioproject_count", "markers", "priority", "notes"
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in top:
            w.writerow({
                "dataset": r["binomial"].replace(" ", "_").lower(),
                "canonical_name": r["binomial"],
                "group": r["group"],
                "gbif_taxon_key": r.get("gbif_taxon_key", ""),
                "gbif_occurrence_count": r.get("gbif_occurrence_count", ""),
                "ncbi_sra_count": r.get("ncbi_sra_count", ""),
                "ncbi_biosample_count": r.get("ncbi_biosample_count", ""),
                "ncbi_bioproject_count": r.get("ncbi_bioproject_count", ""),
                "markers": r.get("markers", ""),
                "priority": round(float(r["priority"]), 2),
                "notes": f"genominer top; n_papers={r.get('n_papers','')}; geographies={r.get('geographies','')}",
            })
    return top


def _genetic_coords_path(dataset: str, repo_dir: Path) -> Path | None:
    """Return the population (or sample) coordinate file from genominer's genetic ingestion."""
    base = repo_dir / "data" / "genetic" / dataset
    for suffix in ("_population_coordinates.csv", "_sample_coordinates.csv"):
        p = base / f"{dataset}{suffix}"
        if p.is_file():
            return p
    return None


def prepare_occurrences(
    queue_csv: str,
    spa_results_root: str,
    max_points: int = 120,
    res_deg: float = 0.5,
    skip_existing: bool = False,
    repo_dir: Path | None = None,
):
    """Prepare population coordinates and a GBIF occurrence backdrop for queued species.

    If genetic coordinates were ingested for a dataset (in data/genetic/<dataset>/),
    those are copied to spa_1000_sp/results/<dataset>/<dataset>_all_markers_coordinates.csv
    and a stale *_occurrences.csv is removed so that 39d_sample_gbif_duckdb.R can rebuild
    the background from the genetic populations. Otherwise GBIF occurrences are thinned and
    used directly as a fallback (only appropriate for species without genetic data).
    """
    if repo_dir is None:
        repo_dir = Path(__file__).resolve().parents[2]
    with open(queue_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    prepared = []
    for r in rows:
        taxon_key = r.get("gbif_taxon_key", "")
        if not taxon_key:
            print(f"Skipping {r['canonical_name']}: no GBIF taxon key")
            prepared.append((r, None))
            continue
        out_dir = Path(spa_results_root) / r["dataset"]
        coords_path = out_dir / f"{r['canonical_name'].replace(' ', '_').lower()}_all_markers_coordinates.csv"
        occurrences_path = out_dir / f"{r['dataset']}_occurrences.csv"

        genetic = _genetic_coords_path(r["dataset"], repo_dir)
        if genetic:
            out_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(genetic, coords_path)
            if occurrences_path.is_file():
                occurrences_path.unlink()
            print(f"Copied genetic coords for {r['canonical_name']}: {genetic} -> {coords_path}")
            prepared.append((r, coords_path))
            continue

        if skip_existing and coords_path.is_file():
            print(f"Using existing occurrence coords for {r['canonical_name']}: {coords_path}")
            prepared.append((r, coords_path))
            continue

        try:
            key = int(taxon_key)
        except ValueError:
            print(f"Skipping {r['canonical_name']}: invalid GBIF taxon key {taxon_key!r}")
            prepared.append((r, None))
            continue
        coords_path = build_occurrences.build_for_species(
            canonical_name=r["canonical_name"],
            taxon_key=key,
            out_dir=str(out_dir),
            max_points=max_points,
            res_deg=res_deg,
        )
        prepared.append((r, coords_path))
    return prepared


def write_species_manifest_rows(queue_csv: str, species_csv: str, mode: str = "append"):
    """Merge queued species into spa_1000_sp config/species.csv, updating GBIF taxon keys."""
    with open(queue_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    fieldnames = [
        "dataset", "canonical_name", "group", "gbif_taxon_key", "status",
        "has_genetic_data", "popmap", "snp_vcf", "haps_vcf", "parent_dataset",
        "species_filter", "notes"
    ]

    # Read existing manifest
    manifest_rows: list[dict] = []
    if os.path.exists(species_csv):
        with open(species_csv, newline="", encoding="utf-8") as f:
            manifest_rows = list(csv.DictReader(f))

    # Index by dataset, prefer populated gbif_taxon_key
    by_dataset: dict[str, dict] = {r["dataset"]: r for r in manifest_rows if r.get("dataset")}
    for r in rows:
        ds = r["dataset"]
        if ds in by_dataset:
            existing = by_dataset[ds]
            if not existing.get("gbif_taxon_key") and r.get("gbif_taxon_key"):
                existing["gbif_taxon_key"] = r["gbif_taxon_key"]
                existing["notes"] = f"genominer queue; {r['notes']}"
        else:
            by_dataset[ds] = {
                "dataset": ds,
                "canonical_name": r["canonical_name"],
                "group": r["group"],
                "gbif_taxon_key": r["gbif_taxon_key"],
                "status": "genominer",
                "has_genetic_data": "no",
                "popmap": "",
                "snp_vcf": "",
                "haps_vcf": "",
                "parent_dataset": "",
                "species_filter": r["canonical_name"],
                "notes": f"genominer queue; {r['notes']}",
            }

    with open(species_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in by_dataset.values():
            w.writerow({k: r.get(k, "") for k in fieldnames})
