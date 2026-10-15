# genominer

**genominer** is a reproducible pipeline and database that mines the primary literature to discover species with population-genetic data, geographic coverage, and a usable reference genome for landscape-genetic and glacial-refugia studies.

Live catalog: [genominer.ecoseek.org](https://genominer.ecoseek.org) · Web frontend: [alrobles/genominer-web-devel](https://github.com/alrobles/genominer-web-devel)

It is intended to be a **publishable, queryable database of candidate species** linking:
- papers (PubMed / GBIF Literature),
- genetic markers and sequence archives (NCBI SRA / BioSample / BioProject),
- sample coordinates (when available in BioSample attributes),
- spatial occurrence data (GBIF),
- a recommended reference genome (NCBI Datasets / Assembly).

For each candidate, genominer produces a `candidates` record, a `popmap`, and a flag (`ready_for_spa` / `ready_for_wgs_gst`) that lets downstream pipelines (`spa_1000_sp`) run niche/SPA, temporal connectivity, and `G'ST ~ SPA` regressions automatically.

## What it does

The pipeline currently has three stages:

1. **Literature mining** (`scripts/mine_candidates.py`)
   - Starts from `config/seeds.csv` (taxon × marker × geography combinations).
   - Expands each seed into a family of FTS5 queries.
   - Searches a local PubMed index (~36 M articles) and GBIF Literature (~61 k papers).
   - Detects population-genetic markers: SNP, WGS, RADseq, GBS, microsatellite, mtDNA, allozyme, STR.
   - Extracts binomial species names and filters false positives.
   - Scores and deduplicates candidate species × paper pairs.
   - Enriches each candidate with GBIF occurrence counts/extent and NCBI SRA / BioSample / BioProject counts.
   - Extracts geographic place names from abstracts.
   - Exports `data/genominer.db` and `data/candidates.tsv`.

2. **Genetic metadata ingestion** (`scripts/ingest_genetic_candidates.py`)
   - For each top candidate, queries NCBI SRA and BioSample.
   - Parses `lat_lon` / `geographic location` attributes.
   - Builds a `popmap.tsv` and `population_coordinates.csv` by grouping samples into geographic cells.
   - Produces `data/genetic_ingest_report.tsv`.

3. **Reference genome recommendation** (`scripts/ingest_reference_genomes.py`)
   - Uses the NCBI Datasets API to find the best assembly for the species or its genus.
   - Reports assembly level (Chromosome / Scaffold / Contig), size, N50, chromosome count, and an FTP URL for the genomic FASTA.
   - Produces `data/reference_genome_report.tsv`.

Finally, `scripts/combine_spa_ready_report.py` joins the three reports into `data/spa_ready_candidates.tsv`, flagging species that have:
- `ready_for_spa` = enough coordinates + GBIF coverage + a reference genome + SRA data.
- `ready_for_wgs_gst` = enough WGS samples with coordinates + a reference genome (suitable for `spa_1000_sp` `G'ST` pipeline).

## Quick start

```bash
# Set index paths (defaults fall back to $HOME/work/...)
export PUBMED_INDEX_DIR=$HOME/work/pubmed/index
export GBIF_LITERATURE_DIR=$HOME/work/gbif_literature

# 1. Mine literature
python3 scripts/mine_candidates.py \
  --seeds config/seeds.csv \
  --out data/genominer.db \
  --tsv data/candidates.tsv \
  --enrich-gbif --enrich-ncbi --enrich-places

# 2. Ingest SRA/BioSample popmaps for the top 20 candidates
python3 scripts/ingest_genetic_candidates.py \
  --candidates data/candidates.tsv \
  --top-n 20 \
  --out data/genetic_ingest_report.tsv \
  --genetic-dir data/genetic

# 3. Find a reference genome for each
python3 scripts/ingest_reference_genomes.py \
  --candidates data/candidates.tsv \
  --top-n 20 \
  --out data/reference_genome_report.tsv

# 4. Combine into an actionable SPA-ready table
python3 scripts/combine_spa_ready_report.py \
  --out data/spa_ready_candidates.tsv
```

## Repository layout

- `config/` — seed lists and taxonomic dictionaries (`seeds.csv`, `target_groups.json`).
- `src/` — reusable modules.
  - `search.py`, `query_builder.py` — FTS5 search over PubMed and GBIF Literature.
  - `entities.py` — genetic marker detection and species-name extraction.
  - `gbif_check.py`, `ncbi_check.py` — GBIF and NCBI enrichment.
  - `geoparse.py` — geographic place-name extraction from abstracts.
  - `sra_popmap.py` — SRA/BioSample metadata download and popmap construction.
  - `reference_genomes.py` — NCBI Datasets API reference-genome search and ranking.
  - `export.py`, `spa_scheduler.py` — database export and SPA job scheduling.
- `scripts/` — runnable entry points for each stage.
- `hpc/` — Slurm job wrappers (`mine_candidates.job`, `ingest_genetic_candidates.job`, `schedule_spa.job`).
- `data/genetic/` — per-dataset SRA runinfo, BioSample attributes, sample/population coordinates and popmaps for ~30 candidate species (public NCBI-derived metadata, included as sample data). Large generated outputs (SQLite db, candidate TSVs) are not committed; see `.gitignore`.

## Main outputs

| File | Description |
|------|-------------|
| `data/candidates.tsv` | Species with paper counts, markers, geographies, GBIF occurrence counts, NCBI SRA/BioSample/BioProject counts. |
| `data/genetic_ingest_report.tsv` | SRA run counts, samples with coordinates, number of populations, and popmap paths. |
| `data/reference_genome_report.tsv` | Recommended reference genome per species: accession, assembly level, size, FTP URL. |
| `data/spa_ready_candidates.tsv` | Combined, ranked table with `ready_for_spa` and `ready_for_wgs_gst` flags. |
| `data/genetic/<dataset>/` | Per-candidate SRA runinfo, BioSample attributes, sample coordinates, popmap, and population coordinates. |
| `data/genominer.db` | SQLite database with `papers`, `candidates`, and `paper_candidate` tables. |

## License

MIT
