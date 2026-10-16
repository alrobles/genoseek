#!/usr/bin/env bash
# genoseek mining-cycle driver.
#
# Submits the full genoseek Slurm chain on the HPC cluster:
#   mine_candidates -> ingest_and_combine -> build_database -> snapshot_report
#
# Designed to be invoked from cron on any host with ssh access to the cluster
# (e.g. a sync/serving host), or run directly on the login node.
#
# Suggested cron (biweekly, 1st and 15th at 06:00):
#   0 6 1,15 * *  $HOME/genoseek/scripts/snapshot_cycle.sh >> $HOME/logs/genoseek_cycle.log 2>&1
#
set -euo pipefail

HPC_USER="${HPC_USER:-$USER}"
HPC_HOST="${HPC_HOST:-hpc.example.edu}"
HPC_KEY="${HPC_KEY:-$HOME/.ssh/id_ed25519}"
HPC_REPO="${HPC_REPO:-/home/$HPC_USER/work/genoseek}"

echo "=== $(date -Iseconds) genoseek snapshot cycle ==="

# Local mode only where sbatch exists (the HPC login node); everywhere else
# drive the submission over ssh (e.g. from a sync-host cron).
if command -v sbatch >/dev/null 2>&1 && [[ -d "$HPC_REPO" ]]; then
    cd "$HPC_REPO"
else
    SSH=(ssh -i "$HPC_KEY" -o ConnectTimeout=30 -o StrictHostKeyChecking=no
         "$HPC_USER@$HPC_HOST")
fi

run_remote() {
    if [[ -n "${SSH[*]:-}" ]]; then
        "${SSH[@]}" "$1"
    else
        bash -lc "$1"
    fi
}

echo "--- refreshing repo ---"
run_remote "cd $HPC_REPO && git fetch origin && git merge --ff-only origin/main 2>&1 | tail -3" \
    || echo "[warn] repo update failed; continuing with current checkout"

echo "--- submitting slurm chain ---"
run_remote "cd $HPC_REPO && mkdir -p logs && bash genoseek_submit.sh"

echo "=== cycle submitted ==="
