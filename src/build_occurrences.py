"""Build an occurrence/coordinates CSV from GBIF for a species."""
from __future__ import annotations
import csv
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

GBIF_API = "https://api.gbif.org/v1"


def _get(url: str, retries: int = 3, delay: float = 0.2):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code in (429, 503, 504) and attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
                continue
            return {}
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
                continue
            return {}
    return {}


def _basis_ok(record: dict) -> bool:
    ok = {"HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "MATERIAL_SAMPLE", "OBSERVATION"}
    return record.get("basisOfRecord", "").upper() in {b.upper() for b in ok}


_INTRODUCED = frozenset(t.lower() for t in [
    "introduced", "invasive", "cultivated", "managed", "escaped",
    "transplanted", "released", "captivity", "captiv", "domesticated",
    "planted", "sown", "non-native", "nonnative", "alien", "exotic",
])


def _establishment_ok(record: dict) -> bool:
    val = record.get("establishmentMeans")
    if val is None:
        return True
    v = str(val).strip().lower()
    return not any(term in v for term in _INTRODUCED)


def _occurrence_status_ok(record: dict) -> bool:
    val = record.get("occurrenceStatus")
    if val is None:
        return True
    return str(val).strip().lower() == "present"


def _individual_count_ok(record: dict) -> bool:
    val = record.get("individualCount")
    if val is None:
        return True
    try:
        return float(val) > 0
    except (TypeError, ValueError):
        return True


def _uncertainty_ok(record: dict, max_m: float = 10000.0) -> bool:
    unc = record.get("coordinateUncertaintyInMeters")
    if unc is None:
        return True
    try:
        return float(unc) <= max_m
    except (TypeError, ValueError):
        return True


def _year_ok(record: dict, year_min: int = 1950) -> bool:
    year = record.get("year")
    if year is None:
        return False
    try:
        return int(year) >= year_min
    except (TypeError, ValueError):
        return False


def _thin_records(records: list[dict], res_deg: float, max_points: int):
    """Spatially thin records to one per cell and then cap to max_points."""
    cells: dict[tuple, dict] = {}
    for r in records:
        lon = r.get("decimalLongitude")
        lat = r.get("decimalLatitude")
        if lon is None or lat is None:
            continue
        try:
            lon = float(lon)
            lat = float(lat)
        except (TypeError, ValueError):
            continue
        cell = (round(lon / res_deg), round(lat / res_deg))
        if cell not in cells:
            r["_lon"] = lon
            r["_lat"] = lat
            cells[cell] = r
    kept = list(cells.values())
    if len(kept) > max_points:
        # Deterministic sample by sorting on cell key for reproducibility
        kept.sort(key=lambda x: (x["_lat"], x["_lon"]))
        kept = kept[:max_points]
    return kept


def build_for_species(
    canonical_name: str,
    taxon_key: int,
    out_dir: str,
    max_raw: int = 2000,
    res_deg: float = 0.5,
    max_points: int = 150,
    year_min: int = 1950,
    sleep: float = 0.1,
) -> Optional[Path]:
    """Query GBIF and write an occurrence CSV. Returns path or None on failure."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"{canonical_name.replace(' ', '_').lower()}_occurrences.csv"
    coords_csv = out_dir / f"{canonical_name.replace(' ', '_').lower()}_all_markers_coordinates.csv"

    dataset_key = out_csv.stem
    records: list[dict] = []
    offset = 0
    limit = 300
    while len(records) < max_raw:
        url = (
            f"{GBIF_API}/occurrence/search?taxonKey={taxon_key}"
            f"&limit={limit}&offset={offset}&hasCoordinate=true"
            f"&hasGeospatialIssue=false&year={year_min},2025"
        )
        data = _get(url, delay=sleep)
        time.sleep(sleep)
        if not data:
            break
        results = data.get("results", [])
        if not results:
            break
        for r in results:
            if (_basis_ok(r) and _uncertainty_ok(r) and _year_ok(r, year_min)
                    and _establishment_ok(r) and _occurrence_status_ok(r)
                    and _individual_count_ok(r)):
                records.append(r)
        offset += len(results)
        if data.get("endOfRecords"):
            break
        if offset >= max_raw:
            break

    if len(records) < 5:
        print(f"[{canonical_name}] not enough GBIF records: {len(records)}")
        return None

    # Spatial thin to target points
    thinned = _thin_records(records, res_deg, max_points)

    fieldnames = [
        "population", "longitude", "latitude", "species", "taxonKey",
        "year", "basisOfRecord", "establishmentMeans", "occurrenceStatus",
        "individualCount", "source", "n_samples"
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for i, r in enumerate(thinned, start=1):
            w.writerow({
                "population": f"GBIF_{i:04d}",
                "longitude": r["_lon"],
                "latitude": r["_lat"],
                "species": canonical_name,
                "taxonKey": taxon_key,
                "year": r.get("year", ""),
                "basisOfRecord": r.get("basisOfRecord", ""),
                "establishmentMeans": r.get("establishmentMeans", ""),
                "occurrenceStatus": r.get("occurrenceStatus", ""),
                "individualCount": r.get("individualCount", ""),
                "source": "GBIF",
                "n_samples": 1,
            })

    # Coordinates file is identical for nicher-only runs
    coords_csv.write_text(out_csv.read_text())
    return coords_csv
