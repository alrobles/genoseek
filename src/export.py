"""Export mined candidates to SQLite and TSV."""
from __future__ import annotations
import csv
import sqlite3
from pathlib import Path


def _unique_key(row: dict) -> tuple:
    return (row.get("source", ""), row.get("external_id", ""), row.get("doi", ""))


def save(db_path: str, tsv_path: str, records: list[dict], enrichment: dict | None = None):
    """records: list of dicts with paper-level fields and binomials/places/score.
       enrichment: optional binomial -> dict of extra candidate-level fields.
    """
    enrichment = enrichment or {}
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    Path(tsv_path).parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path, isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("BEGIN")
    conn.executescript("""
        DROP TABLE IF EXISTS paper_candidate;
        DROP TABLE IF EXISTS candidates;
        DROP TABLE IF EXISTS papers;
        CREATE TABLE papers (
            paper_id INTEGER PRIMARY KEY AUTOINCREMENT,
            source TEXT,
            external_id TEXT,
            doi TEXT,
            title TEXT,
            abstract TEXT,
            authors TEXT,
            pub_year INTEGER,
            journal TEXT,
            keywords TEXT,
            UNIQUE(source, external_id, doi)
        );
        CREATE TABLE candidates (
            candidate_id INTEGER PRIMARY KEY AUTOINCREMENT,
            binomial TEXT UNIQUE,
            genus TEXT,
            score REAL DEFAULT 0,
            n_papers INTEGER DEFAULT 0,
            markers TEXT,
            geographies TEXT,
            places TEXT,
            gbif_taxon_key INTEGER,
            gbif_status TEXT,
            gbif_occurrence_count INTEGER,
            gbif_n_countries INTEGER,
            gbif_extent_km2 REAL,
            ncbi_sra_count INTEGER,
            ncbi_biosample_count INTEGER,
            ncbi_bioproject_count INTEGER
        );
        CREATE TABLE paper_candidate (
            paper_id INTEGER,
            candidate_id INTEGER,
            score REAL,
            markers TEXT,
            places TEXT,
            PRIMARY KEY (paper_id, candidate_id)
        );
    """)

    # Deduplicate papers by key to avoid ambiguity with INSERT OR IGNORE
    seen_keys: set[tuple] = set()
    paper_rows = []
    for r in records:
        key = _unique_key(r)
        if key in seen_keys:
            continue
        seen_keys.add(key)
        paper_rows.append((
            key[0], key[1], key[2], r.get("title", ""), r.get("abstract", ""),
            r.get("authors", ""), r.get("pub_year") or None,
            r.get("journal", ""), r.get("keywords", ""),
        ))

    conn.executemany(
        """INSERT OR IGNORE INTO papers
           (source, external_id, doi, title, abstract, authors, pub_year, journal, keywords)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        paper_rows,
    )

    # Map paper key -> id
    key_to_id: dict[tuple, int] = {}
    cur = conn.execute(
        "SELECT source, external_id, doi, paper_id FROM papers"
    )
    for source, external_id, doi, paper_id in cur:
        key_to_id[(source, external_id, doi)] = paper_id

    # Build candidates
    candidates: dict[str, dict] = {}
    for r in records:
        key = _unique_key(r)
        paper_id = key_to_id.get(key)
        if paper_id is None:
            continue
        for binomial in r.get("binomials", []):
            genus = binomial.split()[0]
            cand = candidates.setdefault(binomial, {
                "binomial": binomial,
                "genus": genus,
                "score": 0.0,
                "n_papers": 0,
                "markers": set(),
                "geographies": set(),
                "places": set(),
                "paper_ids": set(),
            })
            cand["score"] += r.get("score", 0.0)
            if paper_id not in cand["paper_ids"]:
                cand["n_papers"] += 1
                cand["paper_ids"].add(paper_id)
            for m in r.get("markers", []):
                cand["markers"].add(m)
            if r.get("geography"):
                cand["geographies"].add(r["geography"])
            if r.get("places"):
                for p in r["places"].split(","):
                    if p:
                        cand["places"].add(p)

    candidate_rows = []
    for cand in candidates.values():
        extra = enrichment.get(cand["binomial"], {})
        candidate_rows.append((
            cand["binomial"], cand["genus"], cand["score"], cand["n_papers"],
            ",".join(sorted(cand["markers"])),
            ",".join(sorted(cand["geographies"])),
            ",".join(sorted(cand["places"])),
            extra.get("gbif_taxon_key"), extra.get("gbif_status"),
            extra.get("gbif_occurrence_count"), extra.get("gbif_n_countries"),
            extra.get("gbif_extent_km2"),
            extra.get("ncbi_sra_count"), extra.get("ncbi_biosample_count"),
            extra.get("ncbi_bioproject_count"),
        ))

    conn.executemany(
        """INSERT OR REPLACE INTO candidates
           (binomial, genus, score, n_papers, markers, geographies, places,
            gbif_taxon_key, gbif_status, gbif_occurrence_count,
            gbif_n_countries, gbif_extent_km2,
            ncbi_sra_count, ncbi_biosample_count, ncbi_bioproject_count)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        candidate_rows,
    )

    # Map binomial -> id
    binomial_to_id: dict[str, int] = {}
    cur = conn.execute("SELECT binomial, candidate_id FROM candidates")
    for binomial, candidate_id in cur:
        binomial_to_id[binomial] = candidate_id

    # Insert paper-candidate links
    link_rows = []
    for r in records:
        key = _unique_key(r)
        paper_id = key_to_id.get(key)
        if paper_id is None:
            continue
        for binomial in r.get("binomials", []):
            candidate_id = binomial_to_id.get(binomial)
            if candidate_id is None:
                continue
            link_rows.append((
                paper_id, candidate_id, r.get("score", 0.0),
                ",".join(r.get("markers", [])), r.get("places", ""),
            ))

    conn.executemany(
        """INSERT OR REPLACE INTO paper_candidate
           (paper_id, candidate_id, score, markers, places)
           VALUES (?, ?, ?, ?, ?)""",
        link_rows,
    )
    conn.commit()
    conn.close()

    fieldnames = [
        "binomial", "genus", "score", "n_papers", "markers", "geographies", "places",
        "gbif_taxon_key", "gbif_status", "gbif_match_type",
        "gbif_occurrence_count", "gbif_n_countries", "gbif_extent_km2",
        "ncbi_sra_count", "ncbi_biosample_count", "ncbi_bioproject_count"
    ]
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for cand in sorted(candidates.values(), key=lambda x: x["score"], reverse=True):
            extra = enrichment.get(cand["binomial"], {})
            writer.writerow({
                "binomial": cand["binomial"],
                "genus": cand["genus"],
                "score": round(cand["score"], 3),
                "n_papers": cand["n_papers"],
                "markers": ",".join(sorted(cand["markers"])),
                "geographies": ",".join(sorted(cand["geographies"])),
                "places": ",".join(sorted(cand["places"])),
                "gbif_taxon_key": extra.get("gbif_taxon_key", ""),
                "gbif_status": extra.get("gbif_status", ""),
                "gbif_match_type": extra.get("gbif_match_type", ""),
                "gbif_occurrence_count": extra.get("gbif_occurrence_count", ""),
                "gbif_n_countries": extra.get("gbif_n_countries", ""),
                "gbif_extent_km2": extra.get("gbif_extent_km2", ""),
                "ncbi_sra_count": extra.get("ncbi_sra_count", ""),
                "ncbi_biosample_count": extra.get("ncbi_biosample_count", ""),
                "ncbi_bioproject_count": extra.get("ncbi_bioproject_count", ""),
            })
