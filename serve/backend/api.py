"""EcoSeek Genominer — population-genetic candidate + VCF database server.

Mirrors the EcoSeek Bioclim file-server pattern:
  https://genominer.ecoseek.org  (public via Cloudflare Tunnel)
  http://127.0.0.1:8651          (local Docker)

The server is read-only: it serves the SQLite catalog built by
scripts/build_genominer_database.py and the VCF files stored next to it.
"""
from __future__ import annotations

import csv
from collections import Counter
import html
import io
import mimetypes
import os
import re
import sqlite3
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Path as PathParam, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles

DATA_DIR = Path(os.environ.get("GENOMINER_DATA_DIR", "/data"))
PORT = int(os.environ.get("GENOMINER_PORT", "8651"))
BASE_URL = os.environ.get("GENOMINER_BASE_URL", f"http://127.0.0.1:{PORT}")

DB_PATH = DATA_DIR / "genominer.db"
STATIC_DIR = Path(__file__).parent.parent / "frontend" / "public"

app = FastAPI(title="EcoSeek Genominer", version="0.2.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


def _db() -> sqlite3.Connection:
    if not DB_PATH.exists():
        raise HTTPException(status_code=503, detail="genominer.db not found")
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


def _query(sql: str, params: tuple = ()) -> list[sqlite3.Row]:
    conn = _db()
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _file_size_human(path: Path) -> str:
    size = path.stat().st_size
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def _normalize_dataset(name: str) -> str:
    """Strip WGS/pilot/reduced/full suffixes used in VCF dataset names."""
    return re.sub(r"_(wgs|full|reduced|pilot)$", "", name)


_CONTINENT_KEYWORDS: dict[str, str] = {
    "north america": "North America",
    "central america": "North America",
    "caribbean": "North America",
    "nearctic": "North America",
    "south america": "South America",
    "neotropical": "South America",
    "europe": "Europe",
    "palearctic": "Europe",
    "asia": "Asia",
    "indomalayan": "Asia",
    "oriental": "Asia",
    "africa": "Africa",
    "afrotropical": "Africa",
    "oceania": "Oceania",
    "australia": "Oceania",
    "australasian": "Oceania",
    "antarctica": "Antarctica",
}


def _normalize_geography(value: str | None) -> str:
    """Collapse geography strings to continent-level labels."""
    if not value:
        return ""
    seen: set[str] = set()
    for raw in value.split(","):
        token = raw.strip().lower()
        for key, label in _CONTINENT_KEYWORDS.items():
            if key in token:
                seen.add(label)
    return ", ".join(sorted(seen))


def _vcf_chart_svg() -> str:
    """Return a horizontal bar-chart SVG of VCF counts per species."""
    rows = _query(
        """SELECT binomial, COUNT(*) AS n FROM vcf_outputs
           WHERE binomial IS NOT NULL AND binomial != ''
           GROUP BY binomial ORDER BY n DESC, binomial"""
    )
    if not rows:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="400" height="60" style="background:#fff"><text x="20" y="30" fill="#64748b">No VCF data available</text></svg>'

    data = [(str(r["binomial"]), int(r["n"])) for r in rows]
    max_n = max(n for _, n in data)
    bar_height = 22
    left_margin = 210
    right_margin = 70
    top_margin = 40
    bottom_margin = 24
    width = 720
    height = top_margin + len(data) * bar_height + bottom_margin
    plot_width = width - left_margin - right_margin

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" style="background:#fff">',
        '<style>text{font-family:ui-sans-serif,system-ui,sans-serif;font-size:12px;fill:#334155}.label{fill:#334155;text-anchor:end}.count{fill:#0f172a;font-weight:600}.title{font-size:14px;font-weight:600;fill:#0f172a}</style>',
        f'<text class="title" x="{left_margin}" y="24">VCF files per species ({sum(n for _, n in data)} VCFs, {len(data)} species)</text>',
    ]

    for i, (name, count) in enumerate(data):
        y = top_margin + i * bar_height
        bar_w = (count / max_n) * plot_width if max_n else 0
        label = html.escape(name)
        svg.append(f'<rect x="{left_margin}" y="{y + 4}" width="{bar_w:.1f}" height="{bar_height - 8}" rx="4" fill="var(--primary,#3b82f6)" opacity="0.85"/>')
        svg.append(f'<text class="label" x="{left_margin - 10}" y="{y + bar_height / 2 + 4}">{label}</text>')
        svg.append(f'<text class="count" x="{left_margin + bar_w + 8}" y="{y + bar_height / 2 + 4}">{count}</text>')

    # x-axis ticks
    svg.append(f'<line x1="{left_margin}" y1="{height - bottom_margin}" x2="{left_margin + plot_width}" y2="{height - bottom_margin}" stroke="#cbd5e1"/>')
    for tick in range(0, max_n + 1):
        if tick % max(1, max_n // 5) == 0 or tick == max_n:
            x = left_margin + (tick / max_n) * plot_width
            svg.append(f'<text x="{x}" y="{height - 6}" text-anchor="middle" fill="#64748b" font-size="10">{tick}</text>')

    svg.append("</svg>")
    return "\n".join(svg)


def _dataset_to_binomial() -> dict[str, str]:
    """Build a dataset -> binomial map from genominer tables."""
    mapping: dict[str, str] = {}
    for table in ("genetic_ingest_report", "spa_ready_candidates"):
        try:
            for row in _query(f"SELECT dataset, binomial FROM {table} WHERE dataset IS NOT NULL AND binomial IS NOT NULL"):
                mapping[str(row["dataset"])] = str(row["binomial"])
        except sqlite3.OperationalError:
            pass
    return mapping


def _safe_relative_path(rel: str) -> Path:
    """Resolve a relative path under DATA_DIR and prevent directory traversal."""
    if ".." in rel or rel.startswith("/") or re.search(r"[;<>&|`$]", rel):
        raise HTTPException(status_code=400, detail="Invalid relative path")
    target = (DATA_DIR / rel).resolve()
    # Must be inside DATA_DIR and under vcf_database/files
    try:
        target.relative_to(DATA_DIR)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="Path outside data dir") from exc
    if not str(target).startswith(str(DATA_DIR / "vcf_database" / "files")):
        raise HTTPException(status_code=400, detail="Only vcf_database files may be served")
    return target


# ── v0.2 UI ─────────────────────────────────────────────────────────────────
CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
:root {
  --primary: #166534; --primary-light: #dcfce7; --primary-dark: #052e16;
  --bg: #fafafa; --bg-card: #ffffff; --text: #171717;
  --text-muted: #737373; --border: #e5e5e5; --radius: 12px;
  --shadow: 0 1px 3px rgba(0,0,0,0.1), 0 1px 2px rgba(0,0,0,0.06);
  --font: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
  --mono: 'JetBrains Mono', 'SF Mono', Monaco, monospace;
}
[data-theme="dark"] {
  --bg: #0a0a0a; --bg-card: #171717; --text: #fafafa;
  --text-muted: #a3a3a3; --border: #262626;
  --primary: #22c55e; --primary-light: #14532d; --primary-dark: #052e16;
  --shadow: 0 1px 3px rgba(0,0,0,0.4);
}
* { margin: 0; padding: 0; box-sizing: border-box; }
html { scroll-behavior: smooth; }
body { font-family: var(--font); background: var(--bg); color: var(--text); line-height: 1.6; }
a { color: var(--primary); text-decoration: none; }
a:hover { opacity: .8; }
.container { max-width: 1200px; margin: 0 auto; padding: 2rem 1rem; }
.hero { text-align: center; padding: 3rem 1rem 2rem; }
.hero img { width: 64px; height: 64px; }
.hero h1 { font-size: 2.2rem; margin: .5rem 0; }
.hero p { color: var(--text-muted); max-width: 640px; margin: 0 auto; }
.card { background: var(--bg-card); border: 1px solid var(--border); border-radius: var(--radius); box-shadow: var(--shadow); overflow: hidden; margin: 1.5rem 0; }
.card h2 { font-size: 1.1rem; padding: 1rem 1.25rem; border-bottom: 1px solid var(--border); }
.table-wrap { overflow-x: auto; }
table { width: 100%; border-collapse: collapse; font-size: .9rem; }
th, td { padding: .75rem 1rem; border-bottom: 1px solid var(--border); text-align: left; }
th { font-weight: 600; color: var(--text-muted); font-size: .8rem; text-transform: uppercase; letter-spacing: .03em; }
tr:hover { background: rgba(0,0,0,.02); }
.mono { font-family: var(--mono); font-size: .85rem; }
.badge { display: inline-block; padding: .15rem .5rem; border-radius: 999px; background: var(--primary-light); color: var(--primary-dark); font-size: .75rem; font-weight: 600; }
.dl-btn { display: inline-flex; align-items: center; gap: .25rem; font-size: .75rem; padding: .2rem .5rem; border-radius: 6px; background: var(--primary); color: #fff; margin-right: .25rem; margin-bottom: .25rem; }
.dl-btn.gray { background: #737373; }
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 1rem; margin: 2rem 0; }
.stat { background: var(--bg-card); border: 1px solid var(--border); border-radius: var(--radius); padding: 1rem; text-align: center; box-shadow: var(--shadow); }
.stat .num { font-size: 1.8rem; font-weight: 700; color: var(--primary); }
.stat .label { font-size: .8rem; color: var(--text-muted); }
.footer { text-align: center; padding: 2rem; color: var(--text-muted); font-size: .85rem; }
.theme-toggle { position: fixed; top: 1rem; right: 1rem; border: 1px solid var(--border); background: var(--bg-card); border-radius: 999px; padding: .5rem 1rem; cursor: pointer; font-size: .85rem; }
.facet { display: inline-block; padding: .15rem .5rem; border-radius: 999px; background: var(--primary-light); color: var(--primary-dark); font-size: .75rem; font-weight: 500; margin: .1rem; }
.search-form { display: flex; gap: .5rem; flex-wrap: wrap; align-items: center; padding: 1rem 1.25rem; border-bottom: 1px solid var(--border); }
.search-form input, .search-form select { padding: .4rem .6rem; border: 1px solid var(--border); border-radius: 6px; font-family: var(--font); }
.search-form button { padding: .4rem 1rem; background: var(--primary); color: #fff; border: none; border-radius: 6px; cursor: pointer; }
</style>
"""

THEME_JS = """
<script>
(function(){
  const key='genominer-theme';
  const root=document.documentElement;
  const saved=localStorage.getItem(key);
  if(saved){ root.setAttribute('data-theme', saved); }
  else if(window.matchMedia('(prefers-color-scheme: dark)').matches){ root.setAttribute('data-theme','dark'); }
  window.toggleTheme=function(){
    const next = root.getAttribute('data-theme')==='dark' ? 'light' : 'dark';
    root.setAttribute('data-theme', next); localStorage.setItem(key, next);
  };
})();
</script>
"""


def _html_page(title: str, body: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<link rel="icon" type="image/svg+xml" href="/static/favicon.svg">
{CSS}
</head>
<body>
<button class="theme-toggle" onclick="toggleTheme()">🌓 Theme</button>
{body}
{THEME_JS}
</body>
</html>"""


# ── API endpoints ──────────────────────────────────────────────────────────
@app.get("/api/health")
def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "db_exists": DB_PATH.exists(),
        "db_path": str(DB_PATH),
        "data_dir": str(DATA_DIR),
    }


@app.get("/api/catalog")
def api_catalog(
    q: str = Query(None, description="Search binomial/genus/markers/geographies"),
    genus: str = Query(None),
    has_vcf: bool = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    sort: str = Query("score DESC", description="Column direction, e.g. score DESC, binomial ASC"),
) -> dict[str, Any]:
    where = []
    params: list[Any] = []
    if q:
        where.append("""
            (binomial LIKE ? OR genus LIKE ? OR markers LIKE ?
             OR geographies LIKE ? OR places LIKE ? OR status_flags LIKE ?)
        """)
        like = f"%{q}%"
        params.extend([like, like, like, like, like, like])
    if genus:
        where.append("genus = ?")
        params.append(genus)
    if has_vcf is not None:
        where.append("vcf_relative_path IS NOT NULL" if has_vcf else "vcf_relative_path IS NULL")
    # Validate sort to avoid injection; allow alnum, underscore, space, ASC/DESC
    if not re.match(r"^[A-Za-z0-9_ ]+(?: (ASC|DESC))?$", sort):
        sort = "score DESC"
    count_sql = "SELECT COUNT(*) FROM catalog"
    if where:
        count_sql += " WHERE " + " AND ".join(where)
    total = _query(count_sql, tuple(params))
    total = total[0][0] if total else 0

    sql = "SELECT * FROM catalog"
    if where:
        sql += " WHERE " + " AND ".join(where)
    sql += f" ORDER BY {sort} LIMIT ? OFFSET ?"
    params.extend([limit, offset])
    rows = [dict(r) for r in _query(sql, tuple(params))]
    for r in rows:
        r["geographies"] = _normalize_geography(r.get("geographies"))
    return {"count": total, "limit": limit, "offset": offset, "total": total, "returned": len(rows), "base_url": BASE_URL, "items": rows}


@app.get("/api/candidates")
def api_candidates(
    q: str = Query(None, description="Search binomial/genus/markers/geographies"),
    genus: str = Query(None),
    markers: str = Query(None),
    has_vcf: bool = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    sort: str = Query("score DESC"),
) -> dict[str, Any]:
    return api_catalog(q=q, genus=genus, has_vcf=has_vcf, limit=limit, offset=offset, sort=sort)


@app.get("/api/candidate/{binomial}")
def api_candidate(binomial: str) -> dict[str, Any]:
    rows = _query("SELECT * FROM catalog WHERE binomial = ?", (binomial,))
    if not rows:
        raise HTTPException(status_code=404, detail="Candidate not found")
    items = [dict(r) for r in rows]
    for r in items:
        r["geographies"] = _normalize_geography(r.get("geographies"))
    return {"binomial": binomial, "items": items}


@app.get("/api/facets")
def api_facets() -> dict[str, Any]:
    """Datasette-style facet counts."""
    genus = _query(
        "SELECT genus, COUNT(*) AS n FROM catalog WHERE genus IS NOT NULL GROUP BY genus ORDER BY n DESC LIMIT 50"
    )
    markers_counter: Counter[str] = Counter()
    for row in _query("SELECT markers FROM catalog WHERE markers IS NOT NULL AND markers != ''"):
        for marker in str(row["markers"]).split(","):
            marker = marker.strip()
            if marker:
                markers_counter[marker] += 1
    markers = [{"marker": k, "n": v} for k, v in markers_counter.most_common(50)]
    vcf_status = _query(
        """SELECT CASE WHEN vcf_relative_path IS NOT NULL THEN 'has_vcf' ELSE 'no_vcf' END AS status,
            COUNT(*) AS n FROM catalog GROUP BY status"""
    )
    geo_rows = _query("SELECT geographies FROM catalog WHERE geographies IS NOT NULL AND geographies != ''")
    geo_counter: Counter[str] = Counter()
    for row in geo_rows:
        for label in _normalize_geography(row["geographies"]).split(", "):
            if label:
                geo_counter[label] += 1
    top_geographies = [{"geography": k, "n": v} for k, v in geo_counter.most_common(50)]
    return {
        "genus": [dict(r) for r in genus],
        "markers": [dict(r) for r in markers],
        "vcf_status": [dict(r) for r in vcf_status],
        "geographies": [dict(r) for r in top_geographies],
    }


@app.get("/api/vcfs")
def api_vcfs(
    dataset: str = Query(None),
    binomial: str = Query(None),
    vcf_type: str = Query(None),
    has_sidecar: str = Query(None, description="popmap|gst|fst|metrics"),
) -> dict[str, Any]:
    where = []
    params: list[Any] = []
    if dataset:
        where.append("dataset = ?")
        params.append(dataset)
    if binomial:
        where.append("binomial = ?")
        params.append(binomial)
    if vcf_type:
        where.append("vcf_type = ?")
        params.append(vcf_type)
    if has_sidecar == "popmap":
        where.append("popmap_relative_path IS NOT NULL")
    elif has_sidecar == "gst":
        where.append("gst_relative_path IS NOT NULL")
    elif has_sidecar == "fst":
        where.append("fst_relative_path IS NOT NULL")
    elif has_sidecar == "metrics":
        where.append("metrics_relative_path IS NOT NULL")

    sql = "SELECT * FROM vcf_outputs"
    if where:
        sql += " WHERE " + " AND ".join(where)
    sql += " ORDER BY dataset, vcf_type"
    rows = []
    for r in _query(sql, tuple(params)):
        d = dict(r)
        d["filename"] = Path(d["relative_path"]).name if d.get("relative_path") else Path(d["publish_path"]).name
        rows.append(d)
    return {"count": len(rows), "items": rows}


@app.get("/api/vcf-chart")
def api_vcf_chart() -> Response:
    """SVG bar chart of VCF counts per species."""
    return Response(content=_vcf_chart_svg(), media_type="image/svg+xml")


@app.get("/api/vcf/{dataset}/{vcf_type}/{filename}")
def download_vcf(
    dataset: str,
    vcf_type: str,
    filename: str,
) -> FileResponse:
    if ".." in dataset or ".." in vcf_type or ".." in filename:
        raise HTTPException(status_code=400, detail="Invalid path")
    if re.search(r"[;<>&|`$]", filename):
        raise HTTPException(status_code=400, detail="Invalid filename")

    rows = _query(
        "SELECT relative_path FROM vcf_outputs WHERE dataset = ? AND vcf_type = ?",
        (dataset, vcf_type),
    )
    if not rows:
        raise HTTPException(status_code=404, detail="VCF not found")

    rel_path = ""
    for r in rows:
        candidate = DATA_DIR / r["relative_path"]
        if candidate.name == filename:
            rel_path = r["relative_path"]
            break
    if not rel_path:
        rel_path = rows[0]["relative_path"]

    return _serve_file(rel_path)


@app.get("/api/download/{relative_path:path}")
def download_relative(relative_path: str) -> FileResponse:
    """Generic file download by relative path (must live under vcf_database/files)."""
    return _serve_file(relative_path)


def _serve_file(rel_path: str) -> FileResponse:
    file_path = _safe_relative_path(rel_path)
    if not file_path.exists() or not file_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    media = mimetypes.guess_type(str(file_path))[0] or "application/octet-stream"
    return FileResponse(file_path, media_type=media, filename=file_path.name)


@app.get("/api/vcf_manifest.tsv")
def vcf_manifest_tsv() -> PlainTextResponse:
    rows = [dict(r) for r in _query("SELECT * FROM vcf_outputs")]
    if not rows:
        return PlainTextResponse("")
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return PlainTextResponse(out.getvalue(), media_type="text/tab-separated-values")


# ── HTML landing page ───────────────────────────────────────────────────────
@app.get("/", response_class=HTMLResponse)
def landing(
    q: str = Query(None),
    genus: str = Query(None),
    has_vcf: bool = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
) -> str:
    stats = _query(
        """
        SELECT
            (SELECT COUNT(*) FROM candidates) AS n_candidates,
            (SELECT COUNT(*) FROM spa_ready_candidates WHERE ready_for_spa = '1') AS n_spa_ready,
            (SELECT COUNT(*) FROM vcf_outputs) AS n_vcfs,
            (SELECT COUNT(*) FROM vcf_outputs WHERE vcf_type = 'clean_biallelic_snps') AS n_clean_vcfs,
            (SELECT COUNT(DISTINCT binomial) FROM catalog WHERE vcf_relative_path IS NOT NULL) AS n_species_with_vcf
        """
    )[0]

    catalog_data = api_catalog(q=q, genus=genus, has_vcf=has_vcf, limit=limit, offset=offset, sort="score DESC")
    catalog = catalog_data["items"]
    facets = api_facets()

    facet_html = ""
    facet_html += "<div style='padding:0.5rem 1.25rem;'><strong>Genus</strong> "
    for f in facets["genus"][:15]:
        facet_html += f'<a class="facet" href="/?genus={f["genus"]}">{f["genus"]} ({f["n"]})</a>'
    facet_html += "</div>"
    facet_html += "<div style='padding:0.5rem 1.25rem;'><strong>Markers</strong> "
    for f in facets["markers"][:15]:
        facet_html += f'<span class="facet">{f["marker"]} ({f["n"]})</span>'
    facet_html += "</div>"
    facet_html += "<div style='padding:0.5rem 1.25rem;'><strong>Geography</strong> "
    for f in facets["geographies"][:15]:
        facet_html += f'<span class="facet">{f["geography"]} ({f["n"]})</span>'
    facet_html += "</div>"

    catalog_rows = ""
    for r in catalog:
        links = ""
        if r.get("vcf_relative_path"):
            vcf_name = Path(str(r["vcf_relative_path"])).name
            url = f"/api/vcf/{r['dataset']}/{r['vcf_type']}/{vcf_name}"
            links += f'<a class="dl-btn" href="{url}">VCF</a>'
        if r.get("popmap_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{r["popmap_relative_path"]}">popmap</a>'
        if r.get("gst_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{r["gst_relative_path"]}">G\'ST</a>'
        if r.get("fst_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{r["fst_relative_path"]}">FST</a>'
        if r.get("metrics_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{r["metrics_relative_path"]}">metrics</a>'
        status = " ".join(f'<span class="badge">{s}</span>' for s in str(r.get("status_flags") or "").split(",") if s)
        catalog_rows += f"""
        <tr>
          <td><a href="/api/candidate/{r['binomial']}">{r['binomial']}</a></td>
          <td>{r.get('genus','')}</td>
          <td class="mono">{r.get('score','')}</td>
          <td>{r.get('markers','')}</td>
          <td>{r.get('geographies','')}</td>
          <td>{status}</td>
          <td>{links}</td>
        </tr>"""

    vcf_rows = api_vcfs(dataset=None, binomial=None, vcf_type=None)["items"]
    vcf_table_rows = ""
    for d in vcf_rows:
        url = f"/api/vcf/{d['dataset']}/{d['vcf_type']}/{d['filename']}"
        size = ""
        fpath = DATA_DIR / d["relative_path"] if d.get("relative_path") else None
        if fpath and fpath.exists():
            size = _file_size_human(fpath)
        links = f'<a class="dl-btn" href="{url}">VCF</a>'
        if d.get("popmap_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{d["popmap_relative_path"]}">popmap</a>'
        if d.get("gst_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{d["gst_relative_path"]}">G\'ST</a>'
        if d.get("fst_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{d["fst_relative_path"]}">FST</a>'
        if d.get("metrics_relative_path"):
            links += f'<a class="dl-btn gray" href="/api/download/{d["metrics_relative_path"]}">metrics</a>'
        vcf_table_rows += f"""
        <tr>
          <td class="mono">{d['dataset']}</td>
          <td>{d.get('binomial','')}</td>
          <td>{d['vcf_type']}</td>
          <td>{d.get('n_samples','')}</td>
          <td>{d.get('n_variants','')}</td>
          <td class="mono">{size}</td>
          <td>{links}</td>
        </tr>"""

    search_value = q or ""
    genus_value = genus or ""
    has_vcf_checked = "checked" if has_vcf else ""
    body = f"""
    <div class="hero">
      <img src="/static/ecoseek-logo.svg" alt="EcoSeek">
      <h1>EcoSeek Genominer</h1>
      <p>Population-genetic candidate species, reference genomes, and clean VCFs for landscape-genetic / refugia analyses.</p>
    </div>
    <div class="container">
      <div class="stats">
        <div class="stat"><div class="num">{stats['n_candidates']}</div><div class="label">Candidates</div></div>
        <div class="stat"><div class="num">{stats['n_spa_ready']}</div><div class="label">SPA-ready</div></div>
        <div class="stat"><div class="num">{stats['n_vcfs']}</div><div class="label">VCF files</div></div>
        <div class="stat"><div class="num">{stats['n_species_with_vcf']}</div><div class="label">Species with VCF</div></div>
      </div>

      <div class="card">
        <h2>VCF performance</h2>
        <div style="padding:1rem 1.25rem; overflow:auto;">
          <img src="/api/vcf-chart" alt="VCF files per species" style="max-width:100%;height:auto;">
        </div>
      </div>

      <div class="card">
        <h2>Facets</h2>
        {facet_html}
      </div>

      <div class="card">
        <h2>Species catalog</h2>
        <form class="search-form" method="get" action="/">
          <input type="text" name="q" value="{search_value}" placeholder="Search species, genus, markers, geography...">
          <input type="text" name="genus" value="{genus_value}" placeholder="Genus">
          <label><input type="checkbox" name="has_vcf" value="true" {has_vcf_checked}> Has VCF</label>
          <input type="number" name="limit" value="{limit}" min="1" max="1000" style="width:70px">
          <button type="submit">Filter</button>
          <a href="/" style="font-size:.85rem">clear</a>
        </form>
        <div class="table-wrap">
          <table>
            <thead>
              <tr><th>Species</th><th>Genus</th><th>Score</th><th>Markers</th><th>Continents</th><th>Flags</th><th>Downloads</th></tr>
            </thead>
            <tbody>
              {catalog_rows}
            </tbody>
          </table>
        </div>
      </div>

      <div class="card">
        <h2>VCF files</h2>
        <div class="table-wrap">
          <table>
            <thead>
              <tr><th>Dataset</th><th>Binomial</th><th>Type</th><th>Samples</th><th>Variants</th><th>Size</th><th>Downloads</th></tr>
            </thead>
            <tbody>
              {vcf_table_rows}
            </tbody>
          </table>
        </div>
      </div>

      <div class="card">
        <h2>API</h2>
        <div style="padding:1rem 1.25rem; font-family: var(--mono); font-size:.85rem;">
          <p>GET <a href="/api/catalog">/api/catalog?q=...&amp;genus=...&amp;has_vcf=true&amp;limit=100</a></p>
          <p>GET <a href="/api/candidates">/api/candidates</a></p>
          <p>GET <a href="/api/facets">/api/facets</a></p>
          <p>GET <a href="/api/vcfs">/api/vcfs</a></p>
          <p>GET <a href="/api/vcf_manifest.tsv">/api/vcf_manifest.tsv</a></p>
          <p>GET /api/download/&lt;relative_path&gt;</p>
        </div>
      </div>
    </div>
    <div class="footer">
      <a href="https://ecoseek.org">EcoSeek</a> ·
      <a href="https://github.com/alrobles/genominer">Source</a>
    </div>
    """
    return _html_page("EcoSeek Genominer", body)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=PORT)
