#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/work/genominer}"
SPA_REPO="${SPA_REPO:-$HOME/work/spa_1000_sp}"
QUEUE="${QUEUE:-data/spa_ready_queue.csv}"
LOG="${LOG:-${REPO_DIR}/schedule_spa_ready.log}"

cd "${REPO_DIR}"

# Ensure latest scheduler scripts are in place
git fetch origin
git checkout origin/devin/schedule-spa-ready -- scripts/schedule_spa.py src/spa_scheduler.py

# Pull latest candidate tables via relay repo (not tracked in genominer to avoid large data commits)
RELAY_DIR="${RELAY_DIR:-/tmp/genrelay_$$}"
rm -rf "${RELAY_DIR}"
git clone https://github.com/alrobles/test.git "${RELAY_DIR}" --depth 1
cp "${RELAY_DIR}/genominer_data/candidates.tsv" data/candidates.tsv
cp "${RELAY_DIR}/genominer_data/spa_ready_candidates.tsv" data/spa_ready_candidates.tsv
rm -rf "${RELAY_DIR}"

# Pull the fixed generic nicher script from spa_1000_sp main (it removes stale
# --exclude nodes and respects SPA_STOP_BEFORE_GST=1 for SPA-only runs).
if [[ -d "${SPA_REPO}/.git" ]]; then
  git -C "${SPA_REPO}" fetch origin
  git -C "${SPA_REPO}" checkout origin/main -- hpc/generic_nicher_optimization_pipeline.sh
fi

# Subset of species flagged as ready_for_spa=1 in spa_ready_candidates.tsv
SPECIES=(
  "Desmognathus wrighti"
  "Lepus americanus"
  "Lepus californicus"
  "Lepus europaeus"
  "Lepus timidus"
  "Ursus americanus"
  "Ischnura elegans"
  "Picea glauca"
  "Plethodon cinereus"
  "Plethodon teyahalee"
  "Quercus rubra"
)

# Build a TSV with only the SPA-ready species from candidates.tsv
# (filter by the species names listed in spa_ready_candidates.tsv with ready_for_spa=1)
awk -F'\t' 'NR==1 {print; next}' data/candidates.tsv > data/spa_ready_only.tsv
for sp in "${SPECIES[@]}"; do
  awk -F'\t' -v sp="${sp}" '$1==sp' data/candidates.tsv >> data/spa_ready_only.tsv
done

# Ensure all 11 species were matched
n_selected=$(($(wc -l < data/spa_ready_only.tsv) - 1))
if [[ "${n_selected}" -ne "${#SPECIES[@]}" ]]; then
  echo "WARNING: only ${n_selected}/${#SPECIES[@]} species matched in candidates.tsv" >&2
fi

echo "Selected $(($(wc -l < data/spa_ready_only.tsv) - 1)) species from candidates.tsv"

# Launch scheduler in the background; it will fetch occurrences for species without coords
# and submit sbatch jobs for all 11 species.
nohup python3 scripts/schedule_spa.py \
  --candidates data/spa_ready_only.tsv \
  --groups config/target_groups.json \
  --queue "${QUEUE}" \
  --spa-repo "${SPA_REPO}" \
  --min-occurrences 0 \
  --min-sra 0 \
  --top-n 20 \
  --max-points 120 \
  --res-deg 0.5 \
  --partition sixhour \
  --clim-dataset Krapp2021 \
  --skip-existing-occurrences \
  --submit > "${LOG}" 2>&1 < /dev/null &
disown

echo "schedule_spa.py started in background; log: ${LOG}"
sleep 2
tail -n 5 "${LOG}" 2>/dev/null || true
