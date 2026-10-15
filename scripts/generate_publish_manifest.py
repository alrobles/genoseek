#!/usr/bin/env python3
"""Generate a list of relative paths that must be mirrored from a publish route.

The list is derived from genominer.db (vcf_outputs and sidecars) so rsync can
copy only the files the catalog actually references, skipping scratch temp data.
"""
from __future__ import annotations
import argparse
import sqlite3
import sys
from pathlib import Path


SIDE_CAR_COLUMNS = [
    "relative_path",
    "popmap_relative_path",
    "gst_relative_path",
    "fst_relative_path",
    "metrics_relative_path",
]


def main() -> int:
    p = argparse.ArgumentParser(
        description="Write a manifest of relative paths under a publish route."
    )
    p.add_argument("--db", default="data/genominer.db",
                   help="Path to genominer.db")
    p.add_argument("--publish-root", default=".",
                   help="Root directory mirrored by rsync (contains genominer/ and vcf_database/)")
    p.add_argument("--out", default="-",
                   help="Output file ('-' for stdout)")
    p.add_argument("--existing-only", action="store_true",
                   help="Only include files that actually exist under publish-root")
    args = p.parse_args()

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"ERROR: {db_path} not found", file=sys.stderr)
        return 1

    root = Path(args.publish_root).resolve()
    paths: set[str] = set()

    # Core catalog files at the root of the publish tree.
    paths.add("genominer/genominer.db")
    paths.add("genominer/genominer_catalog.tsv")
    paths.add("vcf_database/vcf_manifest.tsv")

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            "SELECT " + ", ".join(SIDE_CAR_COLUMNS) + " FROM vcf_outputs"
        ).fetchall()
    except sqlite3.OperationalError as e:
        print(f"ERROR: cannot query vcf_outputs: {e}", file=sys.stderr)
        return 1
    finally:
        conn.close()

    for row in rows:
        for col in SIDE_CAR_COLUMNS:
            try:
                rel = row[col]
            except IndexError:
                continue
            if rel:
                paths.add(str(rel))

    if args.existing_only:
        existing = {rel for rel in paths if (root / rel).is_file()}
        missing = paths - existing
        if missing:
            print(f"[warn] {len(missing)} referenced files not found under {root}",
                  file=sys.stderr)
        paths = existing

    out_lines = sorted(paths)

    if args.out == "-":
        print("\n".join(out_lines))
    else:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text("\n".join(out_lines) + "\n")
        print(f"[done] {len(out_lines)} paths -> {args.out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
