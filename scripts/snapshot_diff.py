#!/usr/bin/env python3
"""Archive a mining-cycle output set into data/snapshots/<date>/ and diff it
against the previous snapshot.

Detects what a new mining cycle contributed:
  - new_species.tsv: binomials absent from the previous snapshot
  - grown_species.tsv: binomials whose NCBI SRA/BioSample/BioProject counts grew
  - gaps.tsv: candidates with SRA data but no genetic ingest yet
    ("sequences we don't have") and ingested species without a VCF flag
  - snapshot_report.json: counts/deltas for capacity tracking
  - history.tsv (under snapshots/): one appended row per cycle

Usage:
    python3 scripts/snapshot_diff.py --data-dir data --snapshots-dir data/snapshots
    python3 scripts/snapshot_diff.py --date 2026-10-14 --min-sra 5
"""
from __future__ import annotations
import argparse
import csv
import json
import shutil
import sys
from datetime import date
from pathlib import Path

CYCLE_FILES = [
    "candidates.tsv",
    "genetic_ingest_report.tsv",
    "reference_genome_report.tsv",
    "spa_ready_candidates.tsv",
]

GROWTH_COLS = ["ncbi_sra_count", "ncbi_biosample_count", "ncbi_bioproject_count"]


def _load_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _by_binomial(rows: list[dict]) -> dict[str, dict]:
    return {r["binomial"]: r for r in rows if r.get("binomial")}


def _int(row: dict, col: str) -> int:
    try:
        return int(float(row.get(col) or 0))
    except (ValueError, TypeError):
        return 0


def _previous_snapshot(snapshots_dir: Path, today: str) -> Path | None:
    """Most recent snapshot directory older than today's."""
    dirs = sorted(
        d for d in snapshots_dir.iterdir()
        if d.is_dir() and d.name != "latest" and not d.name.startswith(".")
    )
    prev = [d for d in dirs if d.name < today]
    return prev[-1] if prev else None


def _append_history(history_path: Path, row: dict) -> None:
    history_path.parent.mkdir(parents=True, exist_ok=True)
    exists = history_path.exists()
    with open(history_path, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()), delimiter="\t")
        if not exists:
            w.writeheader()
        w.writerow(row)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", default="data")
    p.add_argument("--snapshots-dir", default="data/snapshots")
    p.add_argument("--date", default=None, help="Snapshot label (default: today)")
    p.add_argument("--min-sra", type=int, default=5,
                   help="Min SRA count for a species to count as a usable gap")
    p.add_argument("--vcf-manifest", default=None,
                   help="Optional vcf_manifest.tsv to flag species with a published VCF")
    args = p.parse_args()

    data_dir = Path(args.data_dir)
    snapshots_dir = Path(args.snapshots_dir)
    snap_date = args.date or date.today().isoformat()
    snap_dir = snapshots_dir / snap_date

    cur_candidates = data_dir / "candidates.tsv"
    if not cur_candidates.is_file():
        print(f"ERROR: {cur_candidates} not found", file=sys.stderr)
        return 1

    # --- Archive current cycle outputs -------------------------------------
    snap_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for name in CYCLE_FILES:
        src = data_dir / name
        if src.is_file():
            shutil.copy2(src, snap_dir / name)
            copied.append(name)
        else:
            print(f"[warn] {src} missing, not archived")

    # --- Load current and previous ------------------------------------------
    cur = _by_binomial(_load_tsv(cur_candidates))
    prev_dir = _previous_snapshot(snapshots_dir, snap_date)
    prev = _by_binomial(_load_tsv(prev_dir / "candidates.tsv")) if prev_dir else {}

    # --- Diff ----------------------------------------------------------------
    new_species = [b for b in sorted(cur) if b not in prev]
    grown = []
    for b, row in cur.items():
        if b not in prev:
            continue
        deltas = {c: _int(row, c) - _int(prev[b], c) for c in GROWTH_COLS}
        if any(d > 0 for d in deltas.values()):
            grown.append((b, deltas))

    ingest = _by_binomial(_load_tsv(data_dir / "genetic_ingest_report.tsv"))
    spa_ready = _by_binomial(_load_tsv(data_dir / "spa_ready_candidates.tsv"))

    # Datasets present in the published VCF manifest (e.g. "quercus_rubra_wgs"
    # normalizes to the "quercus_rubra" species dataset name).
    import re
    vcf_datasets: set[str] = set()
    if args.vcf_manifest:
        for mrow in _load_tsv(Path(args.vcf_manifest)):
            ds = (mrow.get("dataset") or "").strip()
            if ds:
                vcf_datasets.add(ds)
                vcf_datasets.add(
                    re.sub(r"_(wgs|full|reduced|pilot|reference)$", "", ds)
                )

    gaps = []
    for b, row in cur.items():
        sra = _int(row, "ncbi_sra_count")
        if sra < args.min_sra:
            continue
        dataset = b.replace(" ", "_").lower()
        ing = ingest.get(b, {})
        has_popmap = bool((ing.get("popmap") or "").strip()) and ing.get("status") == "ok"
        n_coords = _int(ing, "n_with_coords")
        has_vcf = dataset in vcf_datasets if vcf_datasets else None
        reasons = []
        if not has_popmap:
            reasons.append("no_ingest")
        elif n_coords < 3:
            reasons.append("no_coords")
        if has_vcf is False:
            reasons.append("no_vcf")
        if reasons:
            gaps.append({
                "binomial": b,
                "genus": row.get("genus", ""),
                "score": row.get("score", ""),
                "markers": row.get("markers", ""),
                "ncbi_sra_count": row.get("ncbi_sra_count", ""),
                "gbif_occurrence_count": row.get("gbif_occurrence_count", ""),
                "has_popmap": int(has_popmap),
                "n_with_coords": n_coords,
                "has_vcf": "" if has_vcf is None else int(has_vcf),
                "ready_for_spa": spa_ready.get(b, {}).get("ready_for_spa", ""),
                "gap_reason": ";".join(reasons),
            })

    # --- Write outputs --------------------------------------------------------
    for fname, rows in (
        ("new_species.tsv", [cur[b] for b in new_species]),
        ("grown_species.tsv", [
            {**cur[b], **{f"delta_{c}": d[c] for c in GROWTH_COLS},
             **{f"prev_{c}": prev[b].get(c, "") for c in GROWTH_COLS}}
            for b, d in grown
        ]),
        ("gaps.tsv", gaps),
    ):
        if not rows:
            Path(snap_dir / fname).write_text("")
            continue
        fieldnames = list(rows[0].keys())
        with open(snap_dir / fname, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    report = {
        "date": snap_date,
        "previous_snapshot": prev_dir.name if prev_dir else None,
        "archived_files": copied,
        "candidates_total": len(cur),
        "candidates_previous": len(prev),
        "new_species": len(new_species),
        "grown_species": len(grown),
        "gap_species": len(gaps),
        "gap_no_ingest": sum(1 for g in gaps if "no_ingest" in g["gap_reason"]),
        "gap_no_vcf": sum(1 for g in gaps if "no_vcf" in g["gap_reason"]),
        "ingested": sum(1 for r in ingest.values() if r.get("status") == "ok"),
        "spa_ready": sum(
            1 for r in spa_ready.values() if r.get("ready_for_spa") == "1"
        ),
        "new_species_list": new_species,
        "grown_species_list": [b for b, _ in grown],
    }
    with open(snap_dir / "snapshot_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    _append_history(snapshots_dir / "history.tsv", {
        "date": snap_date,
        "candidates_total": report["candidates_total"],
        "new_species": report["new_species"],
        "grown_species": report["grown_species"],
        "gap_species": report["gap_species"],
        "ingested": report["ingested"],
        "spa_ready": report["spa_ready"],
        "previous_snapshot": report["previous_snapshot"] or "",
    })

    print(f"[snapshot] archived {len(copied)} files -> {snap_dir}")
    print(f"[diff] vs {report['previous_snapshot'] or 'NONE'}: "
          f"{len(new_species)} new species, {len(grown)} grown, {len(gaps)} gaps")
    print(f"[history] appended to {snapshots_dir / 'history.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
