#!/usr/bin/env python3
"""Daily progress and scratch-usage report for genoseek/spa_1000_sp.

Run on the HPC head node or from a Slurm job that can see
$HOME/genoseek/scratch. The report is printed to stdout
and appended to a JSON/TSV log file.

Example:
    python3 scripts/progress_report.py \
        --db $HOME/genoseek/scratch/publish/genoseek/genoseek.db \
        --scratch $HOME/genoseek/scratch \
        --out-tsv $HOME/genoseek/scratch/publish/progress.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sqlite3
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _human_bytes(num: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(num) < 1024.0:
            return f"{num:.2f} {unit}"
        num /= 1024.0
    return f"{num:.2f} EB"


def _run(cmd: list[str]) -> str:
    try:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=600,
        ).stdout.strip()
    except subprocess.TimeoutExpired:
        return ""


def _du_summary(scratch: Path) -> dict[str, Any]:
    """Return a quick summary of scratch usage."""
    summary: dict[str, Any] = {
        "scratch_root": str(scratch),
        "total_bytes": None,
        "total_human": None,
        "wgs_bytes": None,
        "wgs_human": None,
        "results_bytes": None,
        "results_human": None,
        "publish_bytes": None,
        "publish_human": None,
        "raw_hotspots": [],
        "top_consumers": [],
    }

    # Total via du -sbL (follow symlinks so publish VCFs are counted)
    total_out = _run(["du", "-sbL", str(scratch)])
    if total_out:
        summary["total_bytes"] = int(total_out.split()[0])
        summary["total_human"] = _human_bytes(summary["total_bytes"])

    # Categorised sizes for known top-level directories
    for label, path in (
        ("scratch", scratch),
        ("results", scratch / "results"),
        ("publish", scratch / "publish"),
    ):
        if path.exists():
            out = _run(["du", "-sbL", str(path)])
            if out:
                try:
                    b = int(out.split()[0])
                    summary[f"{label}_bytes"] = b
                    summary[f"{label}_human"] = _human_bytes(b)
                except ValueError:
                    pass

    # Find raw-data hotspots (temp/cache/bam/fastq directories)
    raw_cmd = [
        "find",
        str(scratch),
        "-maxdepth",
        "2",
        "-type",
        "d",
        "(",
        "-name", "_sra_tmp",
        "-o",
        "-name", "sra-cache",
        "-o",
        "-name", "fastq_tmp",
        "-o",
        "-name", "bams",
        "-o",
        "-name", "fastq",
        "-o",
        "-name", "vcf_parts",
        ")",
        "-print",
    ]
    raw_dirs = _run(raw_cmd).splitlines()
    if raw_dirs:
        du_lines = _run(["du", "-sb"] + raw_dirs).splitlines()
        sizes = []
        for line in du_lines:
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                try:
                    sizes.append((int(parts[0]), parts[1]))
                except ValueError:
                    continue
        sizes.sort(reverse=True)
        summary["raw_hotspots"] = [
            {"path": p, "bytes": b, "human": _human_bytes(b)}
            for b, p in sizes[:20]
        ]

    # Top-level consumers (use shell glob via Python)
    top = []
    for entry in sorted(scratch.glob("*")):
        if entry.is_dir():
            out = _run(["du", "-sbL", str(entry)])
            if out:
                parts = out.split(maxsplit=1)
                try:
                    top.append((int(parts[0]), str(entry)))
                except (ValueError, IndexError):
                    continue
    top.sort(reverse=True)
    summary["top_consumers"] = [
        {"path": p, "bytes": b, "human": _human_bytes(b)}
        for b, p in top[:10]
    ]

    return summary


def _db_progress(db_path: Path) -> dict[str, Any]:
    """Count publishable products in genoseek.db."""
    if not db_path.exists():
        return {"error": f"{db_path} not found"}

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    counts: dict[str, Any] = {}

    queries = {
        "candidates_total": "SELECT COUNT(*) FROM candidates",
        "spa_ready": "SELECT COUNT(*) FROM spa_ready_candidates",
        "reference_genomes": "SELECT COUNT(*) FROM reference_genome_report",
        "with_vcf": "SELECT COUNT(DISTINCT binomial) FROM vcf_outputs WHERE relative_path IS NOT NULL AND relative_path != ''",
        "with_popmap": "SELECT COUNT(DISTINCT binomial) FROM vcf_outputs WHERE popmap_relative_path IS NOT NULL AND popmap_relative_path != ''",
        "with_gst": "SELECT COUNT(DISTINCT binomial) FROM vcf_outputs WHERE gst_relative_path IS NOT NULL AND gst_relative_path != ''",
        "with_fst": "SELECT COUNT(DISTINCT binomial) FROM vcf_outputs WHERE fst_relative_path IS NOT NULL AND fst_relative_path != ''",
        "vcf_files": "SELECT COUNT(*) FROM vcf_outputs WHERE relative_path IS NOT NULL AND relative_path != ''",
    }

    for key, sql in queries.items():
        try:
            counts[key] = cur.execute(sql).fetchone()[0]
        except sqlite3.OperationalError:
            counts[key] = None

    try:
        row = cur.execute(
            "SELECT SUM(CAST(size_bytes AS INTEGER)) FROM vcf_outputs WHERE relative_path IS NOT NULL AND relative_path != ''"
        ).fetchone()
        vcf_bytes = row[0] or 0
        counts["vcf_size_bytes"] = vcf_bytes
        counts["vcf_size_human"] = _human_bytes(vcf_bytes)
    except sqlite3.OperationalError:
        counts["vcf_size_bytes"] = None
        counts["vcf_size_human"] = None

    conn.close()
    return counts


def _report(db_path: Path, scratch: Path) -> dict[str, Any]:
    now = datetime.now(timezone.utc).isoformat()
    report = {
        "timestamp": now,
        "db": _db_progress(db_path),
        "scratch": _du_summary(scratch),
    }
    return report


def _append_tsv(report: dict[str, Any], tsv_path: Path) -> None:
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    flat = {
        "timestamp": report["timestamp"],
        **report["db"],
        "scratch_total_bytes": report["scratch"]["total_bytes"],
        "scratch_total_human": report["scratch"]["total_human"],
        "scratch_results_bytes": report["scratch"]["results_bytes"],
        "scratch_results_human": report["scratch"]["results_human"],
        "scratch_publish_bytes": report["scratch"]["publish_bytes"],
        "scratch_publish_human": report["scratch"]["publish_human"],
    }
    fieldnames = list(flat.keys())
    file_exists = tsv_path.exists()
    with open(tsv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerow(flat)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    default_scratch = os.environ.get("GENOSEEK_SCRATCH", str(Path.home() / "genoseek" / "scratch"))
    default_db = os.environ.get(
        "GENOSEEK_DB",
        default_scratch + "/publish/genoseek/genoseek.db",
    )
    default_out_tsv = os.environ.get(
        "GENOSEEK_PROGRESS_TSV",
        default_scratch + "/publish/progress.tsv",
    )
    p.add_argument("--db", default=default_db)
    p.add_argument("--scratch", default=default_scratch)
    p.add_argument("--out-tsv", default=default_out_tsv)
    p.add_argument("--out-json")
    p.add_argument("--markdown", action="store_true", default=True)
    args = p.parse_args()

    db_path = Path(args.db)
    scratch = Path(args.scratch)
    report = _report(db_path, scratch)

    if args.out_json:
        Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_json, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)

    if args.out_tsv:
        _append_tsv(report, Path(args.out_tsv))

    if args.markdown:
        d = report["db"]
        s = report["scratch"]
        print("# genoseek daily progress")
        print(f"**timestamp:** {report['timestamp']}")
        print()
        print("## catalog")
        print(f"- candidates: {d.get('candidates_total')}")
        print(f"- spa_ready: {d.get('spa_ready')}")
        print(f"- reference_genomes: {d.get('reference_genomes')}")
        print(f"- with_vcf: {d.get('with_vcf')}")
        print(f"- with_popmap: {d.get('with_popmap')}")
        print(f"- with_gst: {d.get('with_gst')}")
        print(f"- with_fst: {d.get('with_fst')}")
        print(f"- vcf_files: {d.get('vcf_files')} ({d.get('vcf_size_human')})")
        print()
        print("## scratch")
        print(f"- total: {s.get('total_human')} ({s.get('total_bytes')} bytes)")
        print(f"- results: {s.get('results_human')}")
        print(f"- publish: {s.get('publish_human')}")
        print()
        print("### raw-data hotspots (top 20)")
        for h in s.get("raw_hotspots", [])[:20]:
            print(f"- `{h['path']}`: {h['human']}")
        print()
        print("### top consumers")
        for c in s.get("top_consumers", [])[:10]:
            print(f"- `{c['path']}`: {c['human']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
