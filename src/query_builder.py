"""Build FTS5 query strings from seeds and target taxonomic groups."""
from __future__ import annotations
import csv
import json
from pathlib import Path


def load_seeds(path: str):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_target_groups(path: str):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _group_for_taxon(taxon: str, groups: dict):
    for gk, g in groups.items():
        if any(taxon.startswith(gen) for gen in g.get("genera", [])):
            return gk
    return None


def build_queries(seed, groups, max_per_seed: int = 100):
    """Yield a list of (FTS5 query string, geography hint) for one seed."""
    taxon = seed.get("taxon", "").strip()
    marker = seed.get("marker", "").strip()
    seed_geo = seed.get("geography", "").strip()

    group_key = _group_for_taxon(taxon, groups)
    group_geo = ",".join(groups[group_key].get("geographies", [])) if group_key else ""

    marker_terms = [m.strip() for m in marker.split("|") if m.strip()]

    queries = []
    # Taxon-only broad
    if taxon:
        queries.append((f'"{taxon}"', seed_geo))
        queries.append((f'"{taxon}" population genetics', seed_geo))
        queries.append((f'"{taxon}" phylogeography', seed_geo))
    # Marker-specific for the seed taxon
    for m in marker_terms[:3]:
        if taxon:
            queries.append((f'"{taxon}" "{m}"', seed_geo))
    # Geography-boosted for the seed taxon
    if taxon and seed_geo:
        # Use only the first region to keep queries focused
        first_region = seed_geo.split(",")[0].strip()
        queries.append((f'"{taxon}" "{first_region}"', seed_geo))

    # Group-level similar taxa: use group geographies, not the seed's specific one
    if group_key:
        for gen in groups[group_key].get("genera", []):
            if gen != taxon.split()[0]:
                queries.append((f'"{gen}" population genetics', group_geo))
                for m in marker_terms[:2]:
                    queries.append((f'"{gen}" "{m}"', group_geo))

    # Deduplicate preserving order, cap
    seen = set()
    out = []
    for q, geo in queries:
        if q not in seen:
            seen.add(q)
            out.append((q, geo))
        if len(out) >= max_per_seed:
            break
    return out
