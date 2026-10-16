#!/bin/bash
set -e

# Configurable work paths. Defaults assume a standard
# $HOME-based layout; set the env vars to override for your cluster.
export GENOSEEK_WORK="${GENOSEEK_WORK:-$HOME/work}"
export SPA_REPO="${SPA_REPO:-$GENOSEEK_WORK/spa_1000_sp}"
export GENOSEEK_SCRATCH="${GENOSEEK_SCRATCH:-$HOME/genoseek/scratch}"
export PUBMED_INDEX_DIR="${PUBMED_INDEX_DIR:-$GENOSEEK_WORK/pubmed/index}"
export GBIF_LITERATURE_DIR="${GBIF_LITERATURE_DIR:-$GENOSEEK_WORK/gbif_literature}"

REPO_ROOT="${REPO_ROOT:-$GENOSEEK_WORK/genoseek}"
cd "$REPO_ROOT"
mkdir -p logs

mine=$(sbatch --parsable --partition=sixhour --chdir="$REPO_ROOT" --time=06:00:00 hpc/mine_candidates.job)
ingest=$(sbatch --parsable --partition=sixhour --dependency=afterok:$mine --chdir="$REPO_ROOT" hpc/ingest_and_combine.job)
build=$(sbatch --parsable --partition=sixhour --dependency=afterok:$ingest \
    --chdir="$REPO_ROOT" --time=01:00:00 hpc/build_database.job)
snapshot=$(sbatch --parsable --partition=sixhour --dependency=afterok:$build \
    --chdir="$REPO_ROOT" hpc/snapshot_report.job)
echo "mine=$mine"
echo "ingest=$ingest"
echo "build=$build"
echo "snapshot=$snapshot"
