"""Search PubMed and GBIF Literature FTS5 indexes."""
from __future__ import annotations
import os
import sqlite3
from pathlib import Path
from typing import Optional


# Allow one env var to set the whole work tree, or set the specific index dirs.
_GENOSEEK_WORK = os.environ.get("GENOSEEK_WORK", str(Path.home() / "work"))
PUBMED_DB = os.path.join(
    os.environ.get("PUBMED_INDEX_DIR", os.path.join(_GENOSEEK_WORK, "pubmed", "index")),
    "pubmed_fts.db",
)
GBIF_DB = os.path.join(
    os.environ.get("GBIF_LITERATURE_DIR", os.path.join(_GENOSEEK_WORK, "gbif_literature")),
    "gbif_literature_fts.db",
)


def search_pubmed(
    query: str,
    limit: int = 200,
    year_min: Optional[int] = None,
    year_max: Optional[int] = None,
) -> list[dict]:
    if not os.path.exists(PUBMED_DB):
        raise FileNotFoundError(
            f"PubMed FTS index not found at {PUBMED_DB}. "
            "Set PUBMED_INDEX_DIR to the directory containing pubmed_fts.db."
        )
    conn = sqlite3.connect(PUBMED_DB, timeout=30)
    conn.row_factory = sqlite3.Row
    sql = """
        SELECT a.pmid, a.doi, a.title, a.abstract, a.journal,
               a.pub_year, a.authors, a.mesh_terms, a.keywords, rank
        FROM articles_fts f
        JOIN articles a ON a.pmid = f.rowid
        WHERE articles_fts MATCH ?
    """
    params: list = [query]
    if year_min:
        sql += " AND a.pub_year >= ?"
        params.append(year_min)
    if year_max:
        sql += " AND a.pub_year <= ?"
        params.append(year_max)
    sql += " ORDER BY rank LIMIT ?"
    params.append(limit)
    try:
        rows = conn.execute(sql, params).fetchall()
    except sqlite3.OperationalError:
        conn.close()
        return []
    conn.close()
    return [dict(r) for r in rows]


def search_gbif_literature(
    query: str,
    limit: int = 200,
    year_min: Optional[int] = None,
    year_max: Optional[int] = None,
) -> list[dict]:
    if not os.path.exists(GBIF_DB):
        raise FileNotFoundError(
            f"GBIF Literature FTS index not found at {GBIF_DB}. "
            "Set GBIF_LITERATURE_DIR to the directory containing gbif_literature_fts.db."
        )
    conn = sqlite3.connect(GBIF_DB, timeout=30)
    conn.row_factory = sqlite3.Row
    sql = """
        SELECT gbif_id, title, abstract, authors, keywords, topics,
               doi, year, source, publisher, language, websites,
               gbif_taxon_key, rank
        FROM literature
        WHERE literature MATCH ?
    """
    params: list = [query]
    conditions = []
    if year_min is not None:
        conditions.append("CAST(year AS INTEGER) >= ?")
        params.append(year_min)
    if year_max is not None:
        conditions.append("CAST(year AS INTEGER) <= ?")
        params.append(year_max)
    if conditions:
        sql += " AND " + " AND ".join(conditions)
    sql += " ORDER BY rank LIMIT ?"
    params.append(limit)
    try:
        rows = conn.execute(sql, params).fetchall()
    except sqlite3.OperationalError:
        conn.close()
        return []
    conn.close()
    return [dict(r) for r in rows]
