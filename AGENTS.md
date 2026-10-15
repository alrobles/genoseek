# Agent Guide: genoseek

## Purpose
Build a reproducible pipeline that mines PubMed and GBIF Literature for species with population-genetic data, producing a publishable candidate database for landscape-genetic / refugia studies.

## Dev setup

No Python dependencies beyond the standard library are required for the core pipeline (`sqlite3`, `csv`, `re`, `json`, `argparse`). Optional analysis can use `pandas` if available.

## How to run

1. Set the full-text index paths (defaults point to your HPC copies):
   - `PUBMED_INDEX_DIR` → directory containing `pubmed_fts.db`
   - `GBIF_LITERATURE_DIR` → directory containing `gbif_literature_fts.db`

2. Build/update a seed list in `config/seeds.csv`.

3. Run the miner:
   ```bash
   python3 scripts/mine_candidates.py
   ```

4. Inspect results:
   ```bash
   sqlite3 data/genominer.db ".schema"
   sqlite3 data/genominer.db "SELECT * FROM candidates ORDER BY score DESC LIMIT 20;"
   ```

## Cluster workflow

For large searches, submit the Slurm wrapper:
```bash
sbatch hpc/mine_candidates.job
```

## Expected files

- `config/seeds.csv` — seed papers / taxon-marker-geography combinations.
- `config/taxa/salamander_genera.txt` — example taxonomic filter.
- `data/candidates.tsv` — flat candidate summary.
- `data/genominer.db` — relational output (papers, candidates, links).

## Snapshot system

Each mining cycle freezes its outputs and diffs them against the previous
cycle so new discoveries are tracked:

- `genominer_submit.sh` chains four Slurm jobs: `mine_candidates` →
  `ingest_and_combine` → `build_database` → `snapshot_report`.
- `scripts/snapshot_diff.py` archives `data/*.tsv` into
  `data/snapshots/<YYYY-MM-DD>/` and writes `new_species.tsv` (binomials not
  seen before), `grown_species.tsv` (SRA/BioSample count deltas),
  `gaps.tsv` (species with SRA data but no popmap/VCF yet), and appends a row
  to `data/snapshots/history.tsv`.
- `scripts/snapshot_cycle.sh` triggers the chain over ssh — run from cron on
  a sync host biweekly (`0 6 1,15 * *`).
- `scripts/snapshot_disk.sh` makes a dated hardlink snapshot (`rsync
  --link-dest`) of the live data dir to an archive disk — runs daily on the
  serving host, keeps the newest 8.
- `scripts/pull_publish.sh` mirrors only DB-referenced files (rsync `-L`
  dereferences the publish-tree symlinks).

## Notes

- The indexes are built by `alrobles/ecoseek-litdump`; do not re-download them here.
- Re-mining only finds new papers when the FTS indexes are refreshed or
  seeds/genera expand — snapshot diffs make that visible in `history.tsv`.
- Avoid destructive operations against `pubmed_fts.db` (read-only access).
- Never commit the SQLite index or large TSV outputs.
