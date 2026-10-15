#!/usr/bin/env python3
"""Safe scratch cleanup for genominer/spa_1000_sp WGS pipelines.

Removes raw/temporary data while keeping final VCFs, popmaps, GST/FST CSVs,
reference genomes and the published catalog. Run with `--dry-run` first.

Example:
    python3 scripts/clean_scratch.py \
        --db $HOME/genominer/scratch/publish/genominer/genominer.db \
        --scratch $HOME/genominer/scratch \
        --dry-run
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import subprocess
from pathlib import Path

RAW_SUBDIRS = (
    "fastq",
    "_sra_tmp",
    "sra-cache",
    "fastq_tmp",
    "bams",
    "vcf",
    "vcf_parts",
    "tmp",
    "bed_chunks",
)

ROOT_INTERMEDIATE_GLOBS = (
    "merged.bam*",
    "raw_variants.vcf.gz*",
    "merged_1kb.*",
    "covered_windows_1kb.bed.gz*",
    "*_haplotypecaller*.log",
)


def _run(cmd: list[str]) -> str:
    try:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=120,
        ).stdout.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return ""


def _active_datasets() -> set[str]:
    """Return dataset names currently in squeue for this user."""
    out = _run(["squeue", "-u", os.environ.get("USER", ""), "--format=%j", "--noheader"])
    active: set[str] = set()
    for name in out.splitlines():
        name = name.strip()
        if not name:
            continue
        # job names look like dataset_stage02_download or dataset_stage04_merge
        parts = name.split("_stage")
        if len(parts) > 1:
            active.add(parts[0])
    return active


def _published_datasets(db_path: Path) -> set[str]:
    """Dataset keys that have a published VCF in genominer.db."""
    if not db_path.exists():
        return set()
    conn = sqlite3.connect(db_path)
    rows = conn.execute(
        "SELECT DISTINCT dataset FROM vcf_outputs WHERE relative_path IS NOT NULL AND relative_path != ''"
    ).fetchall()
    conn.close()
    return {r[0] for r in rows}


def _source_dataset_dirs(db_path: Path) -> set[Path]:
    """Unique source directories referenced by the published VCF manifest."""
    if not db_path.exists():
        return set()
    conn = sqlite3.connect(db_path)
    rows = conn.execute(
        "SELECT DISTINCT source_path FROM vcf_outputs WHERE source_path IS NOT NULL AND source_path != ''"
    ).fetchall()
    conn.close()
    dirs: set[Path] = set()
    for r in rows:
        p = Path(r[0])
        if p.exists():
            dirs.add(p if p.is_dir() else p.parent)
    return dirs


def _size(path: Path) -> int:
    out = _run(["du", "-sbL", str(path)])
    try:
        return int(out.split()[0])
    except (ValueError, IndexError):
        return 0


def _remove(path: Path, dry_run: bool) -> None:
    if dry_run:
        return
    if path.is_dir():
        _run(["rm", "-rf", str(path)])
    elif path.is_file():
        _run(["rm", "-f", str(path)])


def _maybe_remove(path: Path, reason: str, actions: list[tuple[str, int]], dry_run: bool) -> int:
    if not path.exists():
        return 0
    size = _size(path)
    if size == 0:
        return 0
    actions.append((f"{path} [{reason}]", size))
    _remove(path, dry_run)
    return size


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    default_db = os.environ.get(
        "GENOMINER_DB",
        os.environ.get("GENOMINER_SCRATCH", str(Path.home() / "genominer" / "scratch")) + "/publish/genominer/genominer.db",
    )
    default_scratch = os.environ.get("GENOMINER_SCRATCH", str(Path.home() / "genominer" / "scratch"))
    p.add_argument("--db", default=default_db)
    p.add_argument("--scratch", default=default_scratch)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force-stale", nargs="*", default=[], help="Additional dataset dirs to treat as stale/cancelled (remove raw subdirs)")
    args = p.parse_args()

    scratch = Path(args.scratch)
    db_path = Path(args.db)
    published = _published_datasets(db_path)
    active = _active_datasets()
    published_dirs = _source_dataset_dirs(db_path)

    # Discover candidate dataset dirs: *_wgs plus any dir named like a published dataset
    candidate_dirs: set[Path] = set(scratch.glob("*_wgs"))
    for pub_dir in published_dirs:
        if pub_dir.is_relative_to(scratch):
            candidate_dirs.add(pub_dir)

    total_saved = 0
    actions: list[tuple[str, int]] = []

    for dataset_dir in sorted(candidate_dirs):
        dataset = dataset_dir.name
        # normalize *_wgs to match vcf_outputs dataset key
        dataset_key = dataset.removesuffix("_wgs")
        is_published = dataset_key in published or dataset in published
        is_active = (dataset_key in active) or (dataset in active)
        force_stale = dataset_key in args.force_stale or dataset in args.force_stale

        # Published/stale: safe to remove raw subdirs and root intermediates.
        # Active: only safe to remove temp cache subdirs.
        if is_published or force_stale:
            can_clean_raw = True
            can_clean_root_intermediate = True
        elif is_active:
            can_clean_raw = False
            can_clean_root_intermediate = False
        else:
            can_clean_raw = False
            can_clean_root_intermediate = False

        if can_clean_raw:
            for sub in RAW_SUBDIRS:
                subdir = dataset_dir / sub
                if subdir.exists():
                    total_saved += _maybe_remove(subdir, "published" if is_published else "stale", actions, args.dry_run)

        if can_clean_root_intermediate:
            for pattern in ROOT_INTERMEDIATE_GLOBS:
                for item in dataset_dir.glob(pattern):
                    total_saved += _maybe_remove(item, "published" if is_published else "stale", actions, args.dry_run)

        if is_active:
            for sub in ("_sra_tmp", "sra-cache", "fastq_tmp"):
                subdir = dataset_dir / sub
                if subdir.exists():
                    total_saved += _maybe_remove(subdir, "active temp", actions, args.dry_run)

    # Clean any leftover top-level temp/cache directories
    for top in ("sra-cache", "_sra_tmp", "fasterq.tmp"):
        for pth in scratch.glob(f"**/{top}"):
            if pth.is_dir():
                size = _size(pth)
                if size:
                    actions.append((str(pth), size))
                    _remove(pth, args.dry_run)
                    total_saved += size

    mode = "DRY-RUN" if args.dry_run else "DELETED"
    print(f"# scratch cleanup ({mode})")
    print(f"published datasets: {len(published)}")
    print(f"active datasets: {len(active)}")
    if args.force_stale:
        print(f"force-stale datasets: {args.force_stale}")
    print(f"total actions: {len(actions)}")
    print(f"total space: {_human(total_saved)}")
    for label, size in sorted(actions, key=lambda x: x[1], reverse=True)[:30]:
        print(f"- {label}: {_human(size)}")
    return 0


def _human(num: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(num) < 1024.0:
            return f"{num:.2f} {unit}"
        num /= 1024.0
    return f"{num:.2f} EB"


if __name__ == "__main__":
    raise SystemExit(main())
