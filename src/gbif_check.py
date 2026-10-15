"""Query GBIF for occurrence count and spatial extent of a binomial."""
from __future__ import annotations
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request

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


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    R = 6371.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def check(binomial: str, sample_limit: int = 300, sleep: float = 0.1) -> dict:
    """Return dict with taxon_key, status, occurrence_count, n_countries, extent_km2."""
    result = {
        "gbif_taxon_key": None,
        "gbif_status": None,
        "gbif_occurrence_count": None,
        "gbif_n_countries": None,
        "gbif_extent_km2": None,
        "gbif_match_type": None,
    }
    match_url = f"{GBIF_API}/species/match?name={urllib.parse.quote(binomial)}"
    match = _get(match_url)
    time.sleep(sleep)
    if not match or match.get("matchType") == "NONE":
        return result

    result["gbif_taxon_key"] = match.get("usageKey")
    result["gbif_status"] = match.get("status")
    result["gbif_match_type"] = match.get("matchType")

    taxon_key = result["gbif_taxon_key"]
    count_url = f"{GBIF_API}/occurrence/search?taxonKey={taxon_key}&limit=0"
    count_data = _get(count_url)
    time.sleep(sleep)
    result["gbif_occurrence_count"] = count_data.get("count")

    if result["gbif_occurrence_count"]:
        sample_url = (
            f"{GBIF_API}/occurrence/search?taxonKey={taxon_key}"
            f"&limit={sample_limit}&fields=decimalLatitude,decimalLongitude,countryCode"
        )
        sample_data = _get(sample_url)
        time.sleep(sleep)
        records = sample_data.get("results", [])
        if records:
            lats = [r.get("decimalLatitude") for r in records if r.get("decimalLatitude") is not None]
            lons = [r.get("decimalLongitude") for r in records if r.get("decimalLongitude") is not None]
            if lats and lons:
                min_lat, max_lat = min(lats), max(lats)
                min_lon, max_lon = min(lons), max(lons)
                mid_lat = (min_lat + max_lat) / 2
                width = _haversine_km(mid_lat, min_lon, mid_lat, max_lon)
                height = _haversine_km(min_lat, 0, max_lat, 0)
                result["gbif_extent_km2"] = round(width * height, 2)
            countries = {r.get("countryCode") for r in records if r.get("countryCode")}
            result["gbif_n_countries"] = len(countries)

    return result
