#!/usr/bin/env bash
# Dated, hardlink-deduplicated snapshot of the live genominer data dir onto an
# archive disk (e.g. an external USB drive on the serving host).
#
# Source is the LOCAL serving copy (kept fresh by a periodic sync cron), so
# each snapshot freezes whatever the cluster published last. Unchanged files
# are hardlinked to the previous snapshot -> each snapshot only costs the delta.
# Snapshots are independent: deleting an old one never breaks newer ones.
#
# Suggested cron (daily 03:30):
#   30 3 * * *  $HOME/genominer/scripts/snapshot_disk.sh >> $HOME/logs/snapshot_disk.log 2>&1
#
set -euo pipefail

LIVE_DIR="${LIVE_DIR:-$HOME/genominer/data}"
# Point SNAP_ROOT at an external/archive disk.
SNAP_ROOT="${SNAP_ROOT:-$HOME/genominer_snapshots}"
KEEP="${KEEP:-8}"

TODAY="$(date +%F)"
DEST="$SNAP_ROOT/$TODAY"
LOCK="$SNAP_ROOT/.snapshot.lock"

mkdir -p "$SNAP_ROOT"
exec 200>"$LOCK"
flock -n 200 || { echo "$(date -Iseconds) another snapshot_disk.sh is running"; exit 1; }

# Most recent previous snapshot = hardlink basis for unchanged files.
PREV="$(ls -1d "$SNAP_ROOT"/2*/ 2>/dev/null | sort | tail -1 || true)"

LINK_ARGS=()
if [[ -n "$PREV" && "${PREV%/}" != "$DEST" ]]; then
    LINK_ARGS=(--link-dest="$PREV")
    echo "$(date -Iseconds) link-dest basis: $PREV"
fi

echo "$(date -Iseconds) snapshot $LIVE_DIR -> $DEST"
rsync -a --delete "${LINK_ARGS[@]}" "$LIVE_DIR/" "$DEST/"
ln -sfn "$DEST" "$SNAP_ROOT/latest"

SIZE=$(du -sh --apparent-size "$DEST" 2>/dev/null | cut -f1)
REAL=$(du -sh "$DEST" 2>/dev/null | cut -f1)
echo "$(date -Iseconds) done: logical=$SIZE on-disk(delta)=$REAL"

# Retention: keep the newest KEEP snapshots.
ls -1d "$SNAP_ROOT"/2*/ 2>/dev/null | sort | head -n -"$KEEP" | while read -r old; do
    echo "$(date -Iseconds) pruning $old"
    rm -rf "$old"
done
