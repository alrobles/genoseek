#!/usr/bin/env python3
"""Build a single genoseek SQLite database from existing TSVs and VCF outputs.

Usage:
    python3 scripts/build_genoseek_database.py \
        --data-dir data \
        --vcf-manifest /path/to/vcf_manifest.tsv \
        --out-db data/genoseek.db \
        --publish-tsv publish/genoseek_catalog.tsv

The database is meant to be regenerated after each mining/pipeline pass. Do not
commit the SQLite file or the generated catalog TSV; both can be large.
"""
from __future__ import annotations
import argparse
import csv
import os
import shutil
import sqlite3
from pathlib import Path


def _safe_col(name: str) -> str:
    """Make a SQLite-safe column name from a TSV header."""
    return "".join(c if c.isalnum() or c == "_" else "_" for c in name).strip("_")


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?", (name,)
    ).fetchone() is not None


def _normalize_dataset(name: str) -> str:
    """Strip WGS/pilot/reduced/full/reference suffixes used in VCF dataset names."""
    import re

    return re.sub(r"_(wgs|full|reduced|pilot|reference)$", "", name)


def _load_species_csv(repo: Path | None) -> dict[str, str]:
    """Load the spa_1000_sp species.csv mapping dataset -> canonical binomial.

    Also builds reverse mappings from snp_vcf / haps_vcf directory names so
    VCF datasets like `bono2018_reference` can be linked to their species, and
    uses `parent_dataset` to catch sibling datasets (e.g. `bono2018_stacks`).
    """
    mapping: dict[str, str] = {}
    if not repo:
        return mapping
    path = repo / "config" / "species.csv"
    if not path.exists():
        return mapping
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            canonical = row.get("canonical_name", "").strip()
            dataset = row.get("dataset", "").strip()
            if not canonical or not dataset:
                continue
            # The canonical_name may contain multiple species separated by ;
            binomial = canonical.split(";")[0].strip()
            mapping[_normalize_dataset(dataset)] = binomial
            # First dir component of VCF paths often is the VCF dataset name
            for vcf_col in ("snp_vcf", "haps_vcf", "popmap"):
                val = row.get(vcf_col, "") or ""
                for part in val.split(";"):
                    part = part.strip()
                    if not part or part.startswith("hpc/"):
                        continue
                    first = part.split("/")[0].strip()
                    if first:
                        mapping[first] = binomial
            # Sibling datasets under the same parent project (e.g. bono2018_*)
            parent = row.get("parent_dataset", "").strip()
            if parent:
                mapping[parent] = binomial
                norm_parent = _normalize_dataset(parent)
                if norm_parent != parent:
                    mapping[norm_parent] = binomial
    return mapping


def _binomial_for_dataset(
    dataset: str,
    species_map: dict[str, str],
    db_mapping: dict[str, str],
) -> str:
    """Resolve a VCF dataset name to a canonical binomial."""
    if not dataset:
        return ""
    norm = _normalize_dataset(dataset)
    for key in (dataset, norm):
        if key in species_map:
            return species_map[key]
        if key in db_mapping:
            return db_mapping[key]
    # Fall back to prefix matches so sibling datasets (e.g. bono2018_stacks)
    # are mapped to the same species as their parent project.
    for d in (dataset, norm):
        for src in (species_map, db_mapping):
            for key in sorted(src, key=len, reverse=True):
                if d == key or d.startswith(key + "_"):
                    return src[key]
    return ""


def _dataset_to_binomial(conn: sqlite3.Connection) -> dict[str, str]:
    """Build a dataset -> binomial map from genoseek tables."""
    mapping: dict[str, str] = {}
    old = conn.row_factory
    conn.row_factory = sqlite3.Row
    try:
        for table in ("genetic_ingest_report", "spa_ready_candidates"):
            try:
                for row in conn.execute(
                    f"SELECT dataset, binomial FROM {table} WHERE dataset IS NOT NULL AND binomial IS NOT NULL"
                ):
                    mapping[str(row["dataset"])] = str(row["binomial"])
            except sqlite3.OperationalError:
                pass
    finally:
        conn.row_factory = old
    return mapping


def _find_file(candidates: list[Path], prefer: list[str] = None) -> Path | None:
    """Pick the best candidate from a list, preferring names in order."""
    if not candidates:
        return None
    prefer = prefer or []
    for name in prefer:
        for c in candidates:
            if c.name == name:
                return c
    for sub in prefer:
        for c in candidates:
            if sub in c.name:
                return c
    return candidates[0]


def _find_popmap(
    dataset: str,
    binomial: str,
    src_vcf: Path,
    repo: Path | None,
    genoseek_data: Path | None,
) -> Path | None:
    """Locate a popmap TSV for a VCF dataset."""
    norm = _normalize_dataset(dataset)
    places: list[Path] = []
    # Next to the source VCF (often the dataset directory)
    if src_vcf:
        places.extend([src_vcf.parent, src_vcf.parent / "metadata"])
    # spa_1000_sp results locations
    for base_name in (norm, dataset, f"{dataset}_wgs", f"{norm}_wgs"):
        if repo:
            places.append(repo / "results" / base_name)
            places.append(repo / "data" / "external" / base_name)
    # genoseek data/genetic/<binomial>
    if binomial and genoseek_data:
        key = binomial.replace(" ", "_").lower()
        places.append(genoseek_data / "genetic" / key)
    candidates: list[Path] = []
    for place in places:
        if not place or not place.exists():
            continue
        for pat in ("wgs_popmap.tsv", "popmap.tsv", f"{norm}_popmap.tsv", f"{dataset}_popmap.tsv", "*_popmap.tsv"):
            if pat.endswith(".tsv") and "*" not in pat:
                f = place / pat
                if f.exists():
                    candidates.append(f)
            else:
                candidates.extend(sorted(place.glob(pat)))
    # de-duplicate while preserving order; avoid accidentally picking files from the publish tree
    seen: set[Path] = set()
    unique: list[Path] = []
    for c in candidates:
        if c not in seen and "vcf_database" not in str(c):
            seen.add(c)
            unique.append(c)
    return _find_file(unique, prefer=["wgs_popmap.tsv", "popmap.tsv"])


def _find_pairwise_csv(
    dataset: str,
    src_vcf: Path,
    repo: Path | None,
    kind: str = "gst",
) -> Path | None:
    """Locate a pairwise GST/FST CSV for a VCF dataset."""
    norm = _normalize_dataset(dataset)
    places: list[Path] = []
    if src_vcf:
        places.extend([src_vcf.parent, src_vcf.parent / "metadata"])
    if repo:
        for base_name in (norm, dataset, f"{dataset}_wgs", f"{norm}_wgs"):
            places.append(repo / "results" / base_name)
    candidates: list[Path] = []
    for place in places:
        if not place or not place.exists():
            continue
        # Search a few common patterns
        for pat in (
            f"pairwise_{kind}.csv",
            f"*{kind}_pairwise.csv",
            f"*_{kind}_pairwise.csv",
            f"*_{kind}.csv",
            f"*{kind}*.csv",
        ):
            candidates.extend(sorted(place.rglob(pat)))
    seen: set[Path] = set()
    unique: list[Path] = []
    for c in candidates:
        if c not in seen and "vcf_database" not in str(c):
            seen.add(c)
            unique.append(c)
    return _find_file(unique, prefer=[f"pairwise_{kind}.csv"])


def _copy_sidecars(
    conn: sqlite3.Connection,
    publish_root: Path,
    repo: Path | None,
    genoseek_data: Path | None,
):
    """Copy popmap + GST/FST sidecars next to published VCFs and record their relative paths."""
    # Ensure columns exist
    for col in ("popmap_relative_path", "gst_relative_path", "fst_relative_path", "metrics_relative_path"):
        try:
            conn.execute(f"ALTER TABLE vcf_outputs ADD COLUMN {col} TEXT")
        except sqlite3.OperationalError:
            pass

    mapping = _dataset_to_binomial(conn)
    old = conn.row_factory
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            "SELECT rowid, dataset, binomial, source_path, relative_path, popmap_path, pairwise_gst_path, metrics_path FROM vcf_outputs"
        ).fetchall()
    finally:
        conn.row_factory = old

    for row in rows:
        rowid = row["rowid"]
        dataset = row["dataset"] or ""
        binomial = (row["binomial"] or "") or mapping.get(_normalize_dataset(dataset), "")
        src_vcf = Path(row["source_path"]) if row["source_path"] else None
        rel = row["relative_path"] or ""
        if not rel:
            continue
        dst_dir = (publish_root / rel).parent
        dst_dir.mkdir(parents=True, exist_ok=True)

        sidecars: dict[str, Path | None] = {}

        # popmap
        popmap_src = Path(row["popmap_path"]) if row["popmap_path"] else None
        if not popmap_src or not popmap_src.exists():
            popmap_src = _find_popmap(dataset, binomial, src_vcf, repo, genoseek_data)
        if popmap_src and popmap_src.exists():
            popmap_dst = dst_dir / "popmap.tsv"
            shutil.copy2(popmap_src, popmap_dst)
            sidecars["popmap_relative_path"] = str(popmap_dst.relative_to(publish_root))
            print(f"[sidecar] popmap for {dataset}: {popmap_src} -> {popmap_dst}")

        # GST
        gst_src = Path(row["pairwise_gst_path"]) if row["pairwise_gst_path"] else None
        if not gst_src or not gst_src.exists():
            gst_src = _find_pairwise_csv(dataset, src_vcf, repo, kind="gst")
        if gst_src and gst_src.exists():
            gst_dst = dst_dir / "pairwise_gst.csv"
            shutil.copy2(gst_src, gst_dst)
            sidecars["gst_relative_path"] = str(gst_dst.relative_to(publish_root))
            print(f"[sidecar] GST for {dataset}: {gst_src} -> {gst_dst}")

        # FST
        fst_src = _find_pairwise_csv(dataset, src_vcf, repo, kind="fst")
        if fst_src and fst_src.exists():
            fst_dst = dst_dir / "pairwise_fst.csv"
            shutil.copy2(fst_src, fst_dst)
            sidecars["fst_relative_path"] = str(fst_dst.relative_to(publish_root))
            print(f"[sidecar] FST for {dataset}: {fst_src} -> {fst_dst}")

        # Metrics bundle (diversity, PCA, coverage, mapping)
        metrics_src = Path(row["metrics_path"]) if row["metrics_path"] else None
        if metrics_src and metrics_src.exists():
            metrics_tar = dst_dir / "metrics.tar.gz"
            if metrics_src.is_file() and metrics_src.suffixes == [".tar", ".gz"]:
                shutil.copy2(metrics_src, metrics_tar)
            else:
                # Pack directory into a single downloadable tarball
                import tarfile
                if metrics_tar.exists():
                    metrics_tar.unlink()
                with tarfile.open(metrics_tar, "w:gz") as tf:
                    for p in metrics_src.rglob("*"):
                        if p.is_file():
                            tf.add(p, arcname=p.relative_to(metrics_src))
                # Also copy directory for easy browsing
                metrics_dst = dst_dir / "metrics"
                if metrics_dst.exists():
                    shutil.rmtree(metrics_dst, ignore_errors=True)
                shutil.copytree(metrics_src, metrics_dst)
            sidecars["metrics_relative_path"] = str(metrics_tar.relative_to(publish_root))
            print(f"[sidecar] metrics for {dataset}: {metrics_src} -> {metrics_tar}")

        if sidecars:
            sets = ", ".join(f"{k} = ?" for k in sidecars)
            conn.execute(
                f"UPDATE vcf_outputs SET {sets} WHERE rowid = ?",
                list(sidecars.values()) + [rowid],
            )


def _tsv_to_table(conn: sqlite3.Connection, table: str, path: Path, pk: str | None = None):
    """Create a table from a TSV and insert all rows."""
    if not path.exists():
        print(f"[warn] {path} not found, skipping table {table}")
        return 0
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            return 0
        cols = [_safe_col(h) for h in fieldnames]
        col_defs = [f"{c} TEXT" for c in cols]
        if pk and _safe_col(pk) in cols:
            col_defs[cols.index(_safe_col(pk))] += " PRIMARY KEY"
        conn.execute(f"DROP TABLE IF EXISTS {table}")
        conn.execute(f"CREATE TABLE {table} ({', '.join(col_defs)})")
        placeholders = ", ".join(["?"] * len(cols))
        conn.executemany(
            f"INSERT INTO {table} ({', '.join(cols)}) VALUES ({placeholders})",
            [[row.get(h, "") for h in fieldnames] for row in reader],
        )
    print(f"[table] {table}: {conn.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0]} rows from {path}")
    return len(cols)


def _vcf_outputs_table(
    conn: sqlite3.Connection,
    manifest: Path,
    species_map: dict[str, str] | None = None,
    db_mapping: dict[str, str] | None = None,
):
    """Load VCF manifest into vcf_outputs table, adding relative_path + binomial."""
    species_map = species_map or {}
    db_mapping = db_mapping or {}
    conn.execute("DROP TABLE IF EXISTS vcf_outputs")
    if not manifest.exists():
        print(f"[warn] VCF manifest not found: {manifest}")
        return
    with open(manifest, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            return
        # Make publish_path portable: store path relative to the publish root
        # (parent of vcf_database dir) so the same DB works on any sync/archive host.
        publish_root = manifest.parent.parent
        fieldnames.append("relative_path")
        fieldnames.append("binomial")
        cols = [_safe_col(h) for h in fieldnames]
        conn.execute(
            "CREATE TABLE vcf_outputs ("
            + ", ".join(f"{c} TEXT" for c in cols)
            + ")"
        )
        placeholders = ", ".join(["?"] * len(cols))
        rows = []
        for row in reader:
            rel = ""
            if row.get("publish_path"):
                publish_path = Path(row["publish_path"])
                try:
                    rel = str(publish_path.relative_to(publish_root))
                except ValueError:
                    # Manifest may contain absolute paths from another host.
                    # Align on the 'vcf_database' directory suffix so the
                    # relative path still works on the local publish tree.
                    parts = list(publish_path.parts)
                    for i, part in enumerate(parts):
                        if part == "vcf_database":
                            rel = str(Path(*parts[i:]))
                            break
                    if not rel:
                        rel = row["publish_path"]
            row["relative_path"] = rel
            row["binomial"] = _binomial_for_dataset(
                row.get("dataset", ""), species_map, db_mapping
            )
            rows.append([row.get(h, "") for h in fieldnames])
        conn.executemany(
            f"INSERT INTO vcf_outputs ({', '.join(cols)}) VALUES ({placeholders})",
            rows,
        )
    n = conn.execute("SELECT COUNT(*) FROM vcf_outputs").fetchone()[0]
    print(f"[table] vcf_outputs: {n} rows from {manifest}")


def _ensure_candidate_rows(conn: sqlite3.Connection) -> None:
    """Back-fill candidate rows for any VCF binomial not already in candidates."""
    if not _table_exists(conn, "candidates") or not _table_exists(conn, "vcf_outputs"):
        return
    existing = {r[0] for r in conn.execute("SELECT binomial FROM candidates")}
    new_rows = []
    for row in conn.execute("SELECT DISTINCT binomial FROM vcf_outputs WHERE binomial IS NOT NULL AND binomial != ''"):
        binomial = row[0]
        if binomial not in existing:
            genus = binomial.split()[0] if " " in binomial else binomial
            new_rows.append((binomial, genus))
    if not new_rows:
        return
    # Insert with default values for the other catalog columns.
    cols = [
        "binomial", "genus", "score", "n_papers", "markers", "geographies",
        "places", "gbif_taxon_key", "gbif_status", "gbif_match_type",
        "gbif_occurrence_count", "gbif_n_countries", "gbif_extent_km2",
        "ncbi_sra_count", "ncbi_biosample_count", "ncbi_bioproject_count",
    ]
    placeholders = ", ".join(["?"] * len(cols))
    defaults = ("0", "0", "WGS", "", "", "", "", "", "", "", "", "", "0", "0")
    conn.executemany(
        f"INSERT INTO candidates ({', '.join(cols)}) VALUES ({placeholders})",
        [(b, g) + defaults for b, g in new_rows],
    )
    print(f"[candidates] back-filled {len(new_rows)} species from VCF outputs")


def _create_catalog(conn: sqlite3.Connection, publish_path: Path):
    """Export a flat species-level catalog joining candidates, spa_ready, refs and VCFs."""
    publish_path.parent.mkdir(parents=True, exist_ok=True)
    conn.execute("DROP VIEW IF EXISTS catalog")
    conn.execute("""
        CREATE VIEW catalog AS
        SELECT
            c.binomial,
            c.genus,
            CAST(c.score AS REAL) AS score,
            c.n_papers,
            c.markers,
            c.geographies,
            c.places,
            s.dataset,
            s.n_sra_runs,
            s.n_with_coords,
            s.n_populations,
            s.wgs_with_coords,
            s.wgs_total,
            s.radseq_or_gbs_with_coords,
            s.targeted_capture_with_coords,
            s.ref_genome_accession,
            r.ref_total_sequence_length,
            r.ref_ftp_genomic_fasta,
            s.status_flags,
            s.ready_for_spa,
            s.ready_for_wgs_gst,
            v.vcf_type,
            v.n_samples,
            v.n_variants,
            v.publish_path AS vcf_publish_path,
            v.relative_path AS vcf_relative_path,
            v.popmap_relative_path,
            v.gst_relative_path,
            v.fst_relative_path,
            v.metrics_relative_path,
            v.pairwise_gst_path,
            v.md5 AS vcf_md5
        FROM candidates c
        LEFT JOIN spa_ready_candidates s USING (binomial)
        LEFT JOIN reference_genome_report r USING (binomial, dataset)
        LEFT JOIN vcf_outputs v ON v.binomial = c.binomial
        ORDER BY CAST(c.score AS REAL) DESC, c.binomial, v.vcf_type
    """)
    # Export the view to TSV
    cursor = conn.execute("SELECT * FROM catalog")
    headers = [d[0] for d in cursor.description]
    with open(publish_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(cursor)
    print(f"[publish] {publish_path}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", default="data", help="Directory containing genoseek TSVs")
    p.add_argument("--vcf-manifest", default=None, help="Path to spa_1000_sp vcf_manifest.tsv")
    default_repo = os.environ.get(
        "SPA_REPO",
        os.environ.get("GENOSEEK_WORK", str(Path.home() / "work")) + "/spa_1000_sp",
    )
    p.add_argument("--repo", default=default_repo, help="Path to spa_1000_sp repo for results/scratch")
    p.add_argument("--genoseek-data", default="data", help="Path to genoseek data dir (contains genetic/)")
    p.add_argument("--out-db", default="data/genoseek.db", help="Output SQLite path")
    p.add_argument("--publish-tsv", default="publish/genoseek_catalog.tsv", help="Flat catalog TSV")
    args = p.parse_args()

    data_dir = Path(args.data_dir)
    out_db = Path(args.out_db)
    out_db.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(out_db, isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("BEGIN")

    # Core tables generated by the literature-mining pipeline.
    _tsv_to_table(conn, "candidates", data_dir / "candidates.tsv", pk="binomial")
    _tsv_to_table(conn, "spa_ready_candidates", data_dir / "spa_ready_candidates.tsv", pk="binomial")
    _tsv_to_table(conn, "reference_genome_report", data_dir / "reference_genome_report.tsv")
    _tsv_to_table(conn, "genetic_ingest_report", data_dir / "genetic_ingest_report.tsv")

    # Build combined dataset -> binomial map for VCF rows.
    db_mapping = _dataset_to_binomial(conn)
    species_map = _load_species_csv(Path(args.repo) if args.repo else None)

    # VCF outputs produced by the spa_1000_sp WGS pipeline.
    if args.vcf_manifest:
        manifest = Path(args.vcf_manifest)
        _vcf_outputs_table(conn, manifest, species_map=species_map, db_mapping=db_mapping)

        # Remove exact duplicate VCF rows (e.g., a scratch VCF and its copy under
        # repo/results) while keeping the first occurrence.
        conn.execute("""
            DELETE FROM vcf_outputs
            WHERE rowid NOT IN (
                SELECT MIN(rowid)
                FROM vcf_outputs
                GROUP BY COALESCE(binomial, ''), md5, vcf_type
            )
        """)
        # Remove rows that resolve to the SAME published VCF path but carry
        # different metadata (e.g. the file re-exported with a recomputed md5,
        # or listed twice in the manifest under the same publish_path). The
        # relative_path identifies the actual archive file, so keep one row
        # per published path per (binomial, vcf_type).
        conn.execute("""
            DELETE FROM vcf_outputs
            WHERE relative_path != ''
              AND rowid NOT IN (
                  SELECT MIN(rowid)
                  FROM vcf_outputs
                  WHERE relative_path != ''
                  GROUP BY COALESCE(binomial, ''), vcf_type, relative_path
              )
        """)

        # publish_root is the parent of the vcf_database directory.
        publish_root = manifest.parent.parent
        _copy_sidecars(
            conn,
            publish_root,
            Path(args.repo) if args.repo else None,
            Path(args.genoseek_data) if args.genoseek_data else None,
        )
    else:
        print("[info] No --vcf-manifest provided, skipping vcf_outputs table")

    indexes = [
        ("idx_spa_ready_binomial", "spa_ready_candidates", "binomial"),
        ("idx_ref_binomial", "reference_genome_report", "binomial"),
        ("idx_vcf_dataset", "vcf_outputs", "dataset"),
        ("idx_vcf_binomial", "vcf_outputs", "binomial"),
    ]
    for idx_name, table, col in indexes:
        if _table_exists(conn, table):
            conn.execute(f"CREATE INDEX IF NOT EXISTS {idx_name} ON {table}({col})")
        else:
            print(f"[warn] table {table} missing, skipping index {idx_name}")

    _ensure_candidate_rows(conn)

    try:
        _create_catalog(conn, Path(args.publish_tsv))
    except sqlite3.OperationalError as exc:
        print(f"[warn] catalog view could not be created: {exc}")

    conn.commit()
    # Disable WAL so the .db file is self-contained and can be served from
    # a read-only volume (e.g., an external archive drive).
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    conn.close()
    print(f"[done] {out_db}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
