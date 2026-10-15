"""Download SRA run metadata and BioSample attributes, then build a popmap
with coordinates for species that have population-genetic sampling data.

Uses only the Python standard library (no pandas).
"""
from __future__ import annotations
import csv
import io
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Optional

from src import ncbi_check

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url: str, retries: int = 3, delay: float = 0.2) -> str:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (429, 503, 504) and attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
                continue
            return ""
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
                continue
            return ""
    return ""


def _api_key() -> str:
    return ncbi_check._api_key()


def parse_lat_lon(text: str) -> tuple[Optional[float], Optional[float]]:
    """Best-effort parsing of lat/lon strings.

    Handles formats such as:
        '45.67 N 85.53 W'
        '40.5, -120.2'
        '40.5;-120.2'
        '-35.5 150.2'
    Returns (None, None) when parsing fails.
    """
    if not text:
        return None, None
    text = text.strip()
    # Remove common non-numeric qualifiers
    text = re.sub(r"[°'\"]", " ", text)
    text = re.sub(r"(?i)\s*(degrees|deg)\s*", " ", text)

    # Decimal pair separated by comma or semicolon
    m = re.match(r"^\s*(-?\d+\.?\d*)\s*[,;]\s*(-?\d+\.?\d*)\s*$", text)
    if m:
        try:
            lat, lon = float(m.group(1)), float(m.group(2))
            if -90 <= lat <= 90 and -180 <= lon <= 180:
                return lat, lon
        except ValueError:
            pass

    # "lat N/S lon E/W"
    m = re.search(
        r"(-?\d+\.?\d*)\s*([NnSs])[^0-9.-]*(-?\d+\.?\d*)\s*([EeWw])",
        text,
    )
    if m:
        try:
            lat = float(m.group(1))
            lon = float(m.group(3))
            if m.group(2).upper() == "S":
                lat = -lat
            if m.group(4).upper() == "W":
                lon = -lon
            if -90 <= lat <= 90 and -180 <= lon <= 180:
                return lat, lon
        except ValueError:
            pass

    # "45.5 N 120.2 W" (no punctuation)
    tokens = re.findall(r"(-?\d+\.?\d*)\s*([NnSsEeWw]?)", text)
    nums = []
    dirs = []
    for num, direc in tokens:
        if num:
            nums.append(float(num))
            dirs.append(direc.upper())
    if len(nums) >= 2:
        lat, lon = nums[0], nums[1]
        if dirs[0] == "S":
            lat = -lat
        if dirs[1] == "W":
            lon = -lon
        if -90 <= lat <= 90 and -180 <= lon <= 180:
            return lat, lon
    return None, None


def _biosample_url(ids: list[str], api_key: str) -> str:
    id_list = ",".join(urllib.parse.quote(x) for x in ids)
    url = f"{EUTILS}/efetch.fcgi?db=biosample&id={id_list}"
    if api_key:
        url += f"&api_key={api_key}"
    return url


def fetch_biosample_attributes(
    biosample_ids: list[str],
    api_key: str = "",
    sleep: float = 0.34,
    chunk: int = 100,
) -> dict[str, dict]:
    """Fetch and parse BioSample attributes keyed by BioSample accession."""
    if not api_key:
        api_key = _api_key()
    result: dict[str, dict] = {}
    for i in range(0, len(biosample_ids), chunk):
        chunk_ids = biosample_ids[i : i + chunk]
        xml = _get(_biosample_url(chunk_ids, api_key))
        if not xml:
            continue
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            continue
        for bs in root.findall("BioSample"):
            acc = bs.get("accession")
            if not acc:
                continue
            attrs: dict[str, str] = {}
            attrs["accession"] = acc
            title = bs.findtext("Description/Title")
            if title:
                attrs["title"] = title.strip()
            taxon = bs.find("Description/Organism")
            if taxon is not None:
                attrs["taxonomy_name"] = taxon.get("taxonomy_name", "")
            attrs_el = bs.find("Attributes")
            if attrs_el is not None:
                for attr in attrs_el.findall("Attribute"):
                    name = attr.get("attribute_name") or attr.get("harmonized_name")
                    if name and attr.text:
                        attrs[name] = attr.text.strip()
            result[acc] = attrs
        time.sleep(sleep)
    return result


def fetch_sra_runinfo(
    binomial: str,
    api_key: str = "",
    max_runs: int = 2000,
    sleep: float = 0.34,
) -> list[dict]:
    """Fetch SRA runinfo CSV for a binomial.

    Returns a list of dicts keyed by the CSV header.
    """
    if not api_key:
        api_key = _api_key()
    term = f'"{binomial}"[Organism]'
    search_url = (
        f"{EUTILS}/esearch.fcgi?db=sra&term={urllib.parse.quote(term)}"
        f"&retmax={max_runs}&retmode=json"
    )
    if api_key:
        search_url += f"&api_key={api_key}"
    js = _get(search_url)
    time.sleep(sleep)
    if not js:
        return []
    try:
        data = json.loads(js)
    except json.JSONDecodeError:
        return []
    idlist = data.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return []

    rows: list[dict] = []
    # efetch runinfo works best with small batches
    for i in range(0, len(idlist), 100):
        ids = ",".join(idlist[i : i + 100])
        fetch_url = (
            f"{EUTILS}/efetch.fcgi?db=sra&id={ids}&rettype=runinfo&retmode=text"
        )
        if api_key:
            fetch_url += f"&api_key={api_key}"
        csv_text = _get(fetch_url)
        time.sleep(sleep)
        if not csv_text:
            continue
        reader = csv.DictReader(io.StringIO(csv_text))
        for row in reader:
            if row:
                rows.append(row)
    return rows


def _cell(lat: float, lon: float, res_deg: float) -> tuple[float, float]:
    return (round(lat / res_deg) * res_deg, round(lon / res_deg) * res_deg)


def build_popmap(
    binomial: str,
    out_dir: str,
    max_runs: int = 2000,
    coord_res_deg: float = 0.5,
    min_samples_per_pop: int = 1,
    sleep: float = 0.34,
) -> dict:
    """Build popmap for a species from SRA + BioSample metadata.

    Returns a summary dict and writes files under `out_dir`.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    dataset = binomial.replace(" ", "_").lower()

    runinfo = fetch_sra_runinfo(binomial, max_runs=max_runs, sleep=sleep)
    # Write runinfo TSV
    runinfo_path = out_path / f"{dataset}_sra_runinfo.csv"
    if runinfo:
        with open(runinfo_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=runinfo[0].keys())
            w.writeheader()
            w.writerows(runinfo)

    # Accept NCBI (SAMN), EBI (SAME, SAMEA, ...), DDBJ (SAMD) and
    # other BioSample prefixes such as SAMS.
    biosample_re = re.compile(r"^SAM[A-Z]{1,2}\d+", re.IGNORECASE)
    biosample_ids = sorted({
        row["BioSample"] for row in runinfo
        if row.get("BioSample") and biosample_re.match(row["BioSample"])
    })
    bio_attrs = fetch_biosample_attributes(biosample_ids, sleep=sleep)

    attr_path = out_path / f"{dataset}_biosample_attributes.csv"
    if bio_attrs:
        # collect all attribute keys
        keys = set()
        for a in bio_attrs.values():
            keys.update(a.keys())
        fieldnames = ["accession", "taxonomy_name", "title"] + sorted(keys - {"accession", "taxonomy_name", "title"})
        with open(attr_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for acc in sorted(bio_attrs):
                w.writerow({k: bio_attrs[acc].get(k, "") for k in fieldnames})

    # Combine runinfo + biosample attributes to get coordinates per BioSample
    sample_records: list[dict] = []
    seen_biosamples: set[str] = set()
    for row in runinfo:
        bs = row.get("BioSample", "")
        if not bs or bs in seen_biosamples:
            continue
        seen_biosamples.add(bs)
        attrs = bio_attrs.get(bs, {})
        lat, lon = None, None
        for key in ("lat_lon", "geographic location (latitude and longitude)", "latitude and longitude"):
            if key in attrs:
                lat, lon = parse_lat_lon(attrs[key])
                if lat is not None:
                    break
        # Fallback: separate lat/lon keys
        if lat is None and "latitude" in attrs and "longitude" in attrs:
            try:
                lat = float(attrs["latitude"])
                lon = float(attrs["longitude"])
            except (ValueError, TypeError):
                pass
        # Fallback: geo_loc_name might be a place name; skip for coords
        sample_name = row.get("SampleName") or row.get("Sample") or bs
        sample_records.append({
            "sample": sample_name,
            "biosample": bs,
            "srr": row.get("Run", ""),
            "bioproject": row.get("BioProject", ""),
            "library_strategy": row.get("LibraryStrategy", ""),
            "library_source": row.get("LibrarySource", ""),
            "scientific_name": row.get("ScientificName", ""),
            "lat": lat,
            "lon": lon,
            "lat_lon_raw": attrs.get("lat_lon", ""),
            "geo_loc_name": attrs.get("geo_loc_name", ""),
            "collection_date": attrs.get("collection_date", ""),
        })

    coords_path = out_path / f"{dataset}_sample_coordinates.csv"
    with open(coords_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "sample", "biosample", "srr", "bioproject", "library_strategy",
            "library_source", "scientific_name", "lat", "lon", "lat_lon_raw",
            "geo_loc_name", "collection_date"
        ])
        w.writeheader()
        w.writerows(sample_records)

    # Build popmap by grouping rounded coordinates
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for rec in sample_records:
        if rec["lat"] is None or rec["lon"] is None:
            continue
        cell = _cell(rec["lat"], rec["lon"], coord_res_deg)
        grouped[cell].append(rec)

    populations: list[dict] = []
    popmap_rows: list[dict] = []
    for idx, (cell, recs) in enumerate(sorted(grouped.items()), start=1):
        if len(recs) < min_samples_per_pop:
            continue
        pop_name = f"POP_{idx:03d}"
        lats = [r["lat"] for r in recs]
        lons = [r["lon"] for r in recs]
        populations.append({
            "population": pop_name,
            "longitude": round(sum(lons) / len(lons), 6),
            "latitude": round(sum(lats) / len(lats), 6),
            "n_samples": len(recs),
        })
        for r in recs:
            popmap_rows.append({
                "sample": r["sample"],
                "biosample": r["biosample"],
                "population": pop_name,
                "longitude": r["lon"],
                "latitude": r["lat"],
            })

    popmap_path = out_path / f"{dataset}_popmap.tsv"
    with open(popmap_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "biosample", "population", "longitude", "latitude"], delimiter="\t")
        w.writeheader()
        w.writerows(popmap_rows)

    pop_coords_path = out_path / f"{dataset}_population_coordinates.csv"
    with open(pop_coords_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["population", "longitude", "latitude", "n_samples"])
        w.writeheader()
        w.writerows(populations)

    summary = {
        "binomial": binomial,
        "dataset": dataset,
        "n_sra_runs": len(runinfo),
        "n_biosamples": len(sample_records),
        "n_with_coords": sum(1 for r in sample_records if r["lat"] is not None),
        "n_populations": len(populations),
        "out_dir": str(out_path),
        "popmap": str(popmap_path),
        "sample_coordinates": str(coords_path),
        "population_coordinates": str(pop_coords_path),
    }
    return summary
