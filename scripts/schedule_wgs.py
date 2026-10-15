#!/usr/bin/env python3
"""Schedule staged WGS -> G'ST pipelines for genominer WGS-ready candidates.

Reads data/spa_ready_candidates.tsv and submits the staged pipeline in
spa_1000_sp/hpc/submit_wgs_staged.py for the top candidates.

Usage:
    python3 scripts/schedule_wgs.py \
        --spa-repo /home/$USER/work/spa_1000_sp \
        --candidates data/spa_ready_candidates.tsv \
        --top-n 2
"""
from __future__ import annotations
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def is_running(dataset: str, spa_repo: str) -> bool:
    """Best-effort check that the dataset already has a pending/running Slurm job."""
    try:
        result = subprocess.run(
            ["squeue", "-u", os.environ.get("USER", ""), "--format=%j", "--noheader"],
            capture_output=True, text=True, timeout=10
        )
        jobs = result.stdout.splitlines()
        prefix = f"{dataset}_stage"
        return any(j.strip().startswith(prefix) or dataset in j for j in jobs)
    except Exception:
        return False


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", default="data/spa_ready_candidates.tsv")
    default_spa_repo = os.environ.get(
        "SPA_REPO",
        os.environ.get("GENOMINER_WORK", str(Path.home() / "work")) + "/spa_1000_sp",
    )
    p.add_argument("--spa-repo", default=default_spa_repo)
    p.add_argument("--top-n", type=int, default=2, help="Number of WGS-ready candidates to submit")
    p.add_argument("--align-max-jobs", type=int, default=8)
    p.add_argument("--call-n-jobs", type=int, default=6)
    p.add_argument("--gpu", action="store_true", help="Use GPU-accelerated stage03/04 (Parabricks)")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    submitter = Path(args.spa_repo) / "hpc" / "submit_wgs_staged.py"
    if not submitter.exists():
        print(f"ERROR: submit script not found: {submitter}", file=sys.stderr)
        return 1

    rows: list[dict] = []
    with open(args.candidates, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    ready = [r for r in rows if int(r.get("ready_for_wgs_gst", "0") or 0) == 1]
    ready.sort(key=lambda r: float(r.get("priority", "0") or 0), reverse=True)
    top = ready[:args.top_n]

    if not top:
        print("No WGS-ready candidates found.")
        return 0

    submitted: list[tuple[str, str, str]] = []
    for r in top:
        dataset = r["dataset"]
        binomial = r["binomial"]
        if is_running(dataset, args.spa_repo):
            print(f"Skipping {dataset}: already has a running/pending Slurm job")
            continue
        cmd = [
            sys.executable, str(submitter),
            "--dataset", dataset,
            "--binomial", binomial,
            "--repo", args.spa_repo,
            "--align-max-jobs", str(args.align_max_jobs),
            "--call-n-jobs", str(args.call_n_jobs),
        ]
        if args.gpu:
            cmd.append("--gpu")
        print("Submitting:", " ".join(cmd))
        if args.dry_run:
            submitted.append((dataset, "DRY", "DRY"))
            continue
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR submitting {dataset}: {result.stderr}", file=sys.stderr)
            continue
        print(result.stdout)
        submitted.append((dataset, binomial, result.stdout.strip()))

    print(f"\nSubmitted {len(submitted)} staged WGS pipelines:")
    for dataset, binomial, jid in submitted:
        print(f"  {dataset} ({binomial}): {jid}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
