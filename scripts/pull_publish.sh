#!/usr/bin/env bash
# Efficient rsync of a genominer publish route from an HPC cluster to a local path.
#
# The script generates a manifest on the HPC side from genominer.db, then
# rsyncs only the listed files. This avoids scanning the entire scratch tree
# and lets us resume interrupted transfers of large VCFs.
#
# Example:
#   HPC_USER=you \
#   HPC_HOST=hpc.example.edu \
#   HPC_PUBLISH_ROOT=/scratch/you/genominer/publish \
#   HPC_GENOMINER_REPO=/home/you/work/genominer \
#   LOCAL_TARGET=$HOME/genominer/data \
#   bash scripts/pull_publish.sh

set -euo pipefail

HPC_USER="${HPC_USER:-$USER}"
HPC_HOST="${HPC_HOST:-hpc.example.edu}"
HPC_PUBLISH_ROOT="${HPC_PUBLISH_ROOT:-/scratch/$HPC_USER/genominer/publish}"
HPC_GENOMINER_REPO="${HPC_GENOMINER_REPO:-/home/$HPC_USER/work/genominer}"
LOCAL_TARGET="${LOCAL_TARGET:-$HOME/genominer/data}"

LOCK_FILE="${LOCAL_TARGET}/.pull_publish.lock"

# Optional flags
DRY_RUN="${DRY_RUN:-0}"
NO_DELETE="${NO_DELETE:-0}"

# -L dereferences symlinks: the HPC publish tree stores VCFs as links into
# the results dirs, and we want real files on this side.
RSYNC_FLAGS=(-avLP --inplace --partial)
if [[ "$DRY_RUN" == "1" ]]; then
    RSYNC_FLAGS+=("-n")
fi
if [[ "$NO_DELETE" != "1" ]]; then
    RSYNC_FLAGS+=(--delete --delete-excluded)
fi

mkdir -p "$LOCAL_TARGET"

# Generate manifest on the HPC side.
# This list contains relative paths like:
#   genominer/genominer.db
#   genominer/genominer_catalog.tsv
#   vcf_database/vcf_manifest.tsv
#   vcf_database/files/<dataset>/<vcf_type>/...
MANIFEST_REMOTE="/tmp/publish_manifest_${HPC_USER}.txt"
echo "=== Generating manifest on $HPC_HOST ==="
ssh "${HPC_USER}@${HPC_HOST}" \
    "python3 ${HPC_GENOMINER_REPO}/scripts/generate_publish_manifest.py \
        --db ${HPC_PUBLISH_ROOT}/genominer/genominer.db \
        --publish-root ${HPC_PUBLISH_ROOT} \
        --existing-only \
        --out ${MANIFEST_REMOTE}"

# Copy manifest to local for inspection / rsync --files-from.
MANIFEST_LOCAL="${LOCAL_TARGET}/.publish_manifest.txt"
scp -q "${HPC_USER}@${HPC_HOST}:${MANIFEST_REMOTE}" "$MANIFEST_LOCAL"

# Count what we are about to transfer.
N_FILES=$(wc -l < "$MANIFEST_LOCAL" | tr -d ' ')
echo "=== Manifest has $N_FILES files ==="

# Pull only the files the catalog references.
echo "=== Starting rsync to $LOCAL_TARGET ==="
# Use a file lock so a second run started by cron does not race.
exec 200>"$LOCK_FILE"
if ! flock -n 200; then
    echo "ERROR: another pull_publish.sh is already running" >&2
    exit 1
fi

rsync "${RSYNC_FLAGS[@]}" \
    --files-from="$MANIFEST_LOCAL" \
    "${HPC_USER}@${HPC_HOST}:${HPC_PUBLISH_ROOT}/" \
    "$LOCAL_TARGET/"

# Keep the manifest as a record.
echo "=== Done. Copied manifest to $MANIFEST_LOCAL ==="
