#!/usr/bin/env python3
"""Build the SPA queue, prepare GBIF occurrences, update spa_1000_sp manifest, and submit Slurm jobs.

Usage:
    python3 scripts/schedule_spa.py \
        --candidates data/candidates.tsv \
        --groups config/target_groups.json \
        --queue data/spa_queue.csv \
        --spa-repo /home/$USER/work/spa_1000_sp \
        --submit
"""
from __future__ import annotations
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import spa_scheduler


def main():
    parser = argparse.ArgumentParser(description="Queue and submit SPA jobs from genoseek candidates.")
    parser.add_argument("--candidates", default="data/candidates.tsv")
    parser.add_argument("--groups", default="config/target_groups.json")
    parser.add_argument("--queue", default="data/spa_queue.csv")
    default_spa_repo = os.environ.get(
        "SPA_REPO",
        os.environ.get("GENOSEEK_WORK", str(Path.home() / "work")) + "/spa_1000_sp",
    )
    parser.add_argument("--spa-repo", default=default_spa_repo)
    parser.add_argument("--min-occurrences", type=int, default=5000)
    parser.add_argument("--min-sra", type=int, default=5)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--max-points", type=int, default=120)
    parser.add_argument("--res-deg", type=float, default=0.5)
    parser.add_argument("--submit", action="store_true", help="Submit Slurm jobs for queued species")
    parser.add_argument("--partition", default="sixhour")
    parser.add_argument("--clim-dataset", default="Krapp2021")
    parser.add_argument("--marker", default="snps", help="Marker name for G'ST regression if a pairwise GST file exists")
    parser.add_argument("--skip-existing-occurrences", action="store_true", help="Do not overwrite existing _all_markers_coordinates.csv")
    args = parser.parse_args()

    # 1. Build queue
    queue = spa_scheduler.build_queue(
        args.candidates, args.groups, args.queue,
        min_occurrences=args.min_occurrences,
        min_sra=args.min_sra,
        top_n=args.top_n,
    )
    print(f"Selected {len(queue)} species for SPA queue: {args.queue}")
    if not queue:
        return

    # 2. Append manifest rows to spa_1000_sp config/species.csv
    species_csv = os.path.join(args.spa_repo, "config", "species.csv")
    spa_scheduler.write_species_manifest_rows(args.queue, species_csv, mode="append")
    print(f"Appended {len(queue)} rows to {species_csv}")

    # 3. Prepare occurrence files (prefer genetic coords ingested by genoseek)
    repo_dir = Path(__file__).resolve().parents[1]
    prepared = spa_scheduler.prepare_occurrences(
        args.queue,
        os.path.join(args.spa_repo, "results"),
        max_points=args.max_points,
        res_deg=args.res_deg,
        skip_existing=args.skip_existing_occurrences,
        repo_dir=repo_dir,
    )

    # 4. Optionally submit Slurm jobs
    if args.submit:
        for row, coords_path in prepared:
            if coords_path is None:
                print(f"Skipping {row['canonical_name']}: no occurrence file")
                continue
            dataset = row["dataset"]
            coords = str(coords_path)
            occurrences = coords.replace("_all_markers_coordinates.csv", "_occurrences.csv")
            spa_results = Path(args.spa_repo) / "results" / dataset
            pairwise_gst = spa_results / f"{dataset}_{args.marker}_pairwise_gst.csv"
            if pairwise_gst.is_file():
                stop_before = "0"
                pairwise_arg = f"SPA_PAIRWISE_GST={pairwise_gst},"
            else:
                stop_before = "1"
                pairwise_arg = ""
            export = (
                f"--export=ALL,SPA_DATASET={dataset},"
                f"SPA_POP_COORDS={coords},"
                f"SPA_SPA_COORDS={coords},"
                f"SPA_OCCURRENCES={occurrences},"
                f"{pairwise_arg}"
                f"SPA_STOP_BEFORE_GST={stop_before},"
                f"SPA_CLIM_DATASET={args.clim_dataset},"
                f"TEMPORAL_MAX_STEP_KM=75,"
                f"THIN_N_MAX={args.max_points},"
                f"THIN_RES_DEG={args.res_deg}"
            )
            cmd = [
                "sbatch",
                export,
                "--partition", args.partition,
                f"{args.spa_repo}/hpc/generic_nicher_optimization_pipeline.sh",
            ]
            print("Submitting:", " ".join(cmd))
            subprocess.run(cmd, check=False)


if __name__ == "__main__":
    main()
