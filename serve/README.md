# GenoSeek Web Server

FastAPI + SQLite file server for the `genominer` population-genetic candidate database.

Mirrors the [EcoSeek Bioclim](https://bioclim.ecoseek.org) deployment pattern:

- `genominer.db` (SQLite) and `vcf_database/` live in the repo's `data/` dir
  (`$HOME/genominer/data` on the serving host), kept in sync with the HPC
  publish route by a periodic `rsync` cron (e.g. every 5 min).
- Dated hardlink snapshots are archived to an external drive by
  `scripts/snapshot_disk.sh` (`$SNAP_ROOT`, default `$HOME/genominer_snapshots`).
- Docker container reads the data read-only and exposes the UI/JSON API on port `8651`.
- A tunnel (e.g. Cloudflare Tunnel) maps the public hostname → `http://127.0.0.1:8651`.

## Build & run locally

```bash
cd genominer/serve
GENOMINER_DATA_DIR=/path/to/genominer_data docker compose up -d
# open http://127.0.0.1:8651
```

## Required data layout

```
/path/to/genominer_data/
├── genominer.db                 # built by scripts/build_genominer_database.py
└── vcf_database/
    ├── vcf_manifest.tsv
    └── files/
        └── <dataset>/
            └── <vcf_type>/
                ├── *.vcf.gz
                ├── popmap.tsv          # optional sidecar
                ├── pairwise_gst.csv    # optional sidecar
                └── pairwise_fst.csv    # optional sidecar
```

## Build database on HPC and publish to the serving host

On the cluster the build writes to a single **publish route** under
`GENOMINER_SCRATCH`:

```bash
export GENOMINER_WORK=/home/$USER/work
export SPA_REPO=$GENOMINER_WORK/spa_1000_sp
export GENOMINER_SCRATCH=$HOME/genominer/scratch

cd $GENOMINER_WORK/genominer
python3 scripts/build_genominer_database.py \
  --data-dir data \
  --repo $SPA_REPO \
  --vcf-manifest $GENOMINER_SCRATCH/publish/vcf_database/vcf_manifest.tsv \
  --genominer-data data \
  --out-db $GENOMINER_SCRATCH/publish/genominer/genominer.db \
  --publish-tsv $GENOMINER_SCRATCH/publish/genominer/genominer_catalog.tsv
```

### Efficient sync to the serving / archive host

`scripts/pull_publish.sh` pulls only the files that `genominer.db` actually
references (DB, catalog, VCF manifest, VCFs and sidecars). It uses
`--inplace --partial` so large VCFs can resume and does not re-scan the whole
scratch tree.

On the serving host:

```bash
export HPC_USER=you
export HPC_HOST=hpc.example.edu
export HPC_PUBLISH_ROOT=/scratch/$HPC_USER/genominer/publish
export HPC_GENOMINER_REPO=/home/$HPC_USER/work/genominer
export LOCAL_TARGET=$HOME/genominer/data

# Dry run first
DRY_RUN=1 bash /path/to/genominer/scripts/pull_publish.sh

# Real pull
bash /path/to/genominer/scripts/pull_publish.sh
```

If the HPC publish route is not yet complete, run the build step above first.

## Cloudflare Tunnel

Run the `genoseek` tunnel as a systemd user service on the serving host:

```bash
systemctl --user status genoseek-tunnel
```

To set up manually:

```bash
cloudflared tunnel create genoseek
cloudflared tunnel route dns <id> genominer.example.org
cloudflared tunnel --config ~/.cloudflared/genoseek.yml run
```

## API

| Endpoint | Description |
|---|---|
| `GET /` | HTML catalog with search, facets, and downloads |
| `GET /api/health` | Health + data-dir status |
| `GET /api/catalog` | Flat catalog (species + VCF join); supports `q`, `genus`, `has_vcf`, `limit`, `offset`, `sort` |
| `GET /api/candidates` | Candidate list with same filters |
| `GET /api/candidate/{binomial}` | Single candidate detail |
| `GET /api/facets` | Datasette-style facet counts (genus, markers, vcf_status, geographies) |
| `GET /api/vcfs` | All VCF files; supports `dataset`, `binomial`, `vcf_type`, `has_sidecar` |
| `GET /api/vcf/{dataset}/{vcf_type}/{filename}` | Download a VCF |
| `GET /api/download/{relative_path}` | Download any sidecar (popmap, GST, FST, index) by relative path |
| `GET /api/vcf_manifest.tsv` | Raw VCF manifest TSV |
