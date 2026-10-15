"""Query NCBI Entrez for sequence/genomic evidence of a binomial."""
from __future__ import annotations
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get_xml(url: str, retries: int = 3, delay: float = 0.2):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
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
    key = os.environ.get("ENTREZ_API_KEY", "")
    if key:
        return key
    path = os.environ.get("ENTREZ_KEY_FILE", os.path.expanduser("~/env/nih-entrez-key"))
    if os.path.exists(path):
        raw = open(path).read().strip()
        if raw.startswith("ENTREZ_API_KEY="):
            return raw.split("=", 1)[1]
        return raw
    return ""


def _count(db: str, term: str, api_key: str) -> int:
    url = f"{EUTILS}/esearch.fcgi?db={db}&term={urllib.parse.quote(term)}&retmax=0"
    if api_key:
        url += f"&api_key={api_key}"
    xml = _get_xml(url)
    if not xml:
        return 0
    try:
        root = ET.fromstring(xml)
        count = root.find("Count")
        return int(count.text) if count is not None and count.text else 0
    except ET.ParseError:
        return 0


_key_warning_shown = False


def check(binomial: str, sleep: float | None = None) -> dict:
    """Return SRA, BioSample and BioProject counts for a binomial."""
    global _key_warning_shown
    api_key = _api_key()
    if sleep is None:
        # 3 eutils requests per call; be polite: 0.34 s with key, 1.0 s without.
        sleep = 0.34 if api_key else 1.0
    elif not api_key and sleep < 0.34:
        sleep = 0.34
    if not api_key and not _key_warning_shown:
        print("[warn] ENTREZ_API_KEY not set; using slower public NCBI rate.")
        _key_warning_shown = True
    term = f'"{binomial}"[Organism]'
    result = {}
    for db in ("sra", "biosample", "bioproject"):
        result[f"ncbi_{db}_count"] = _count(db, term, api_key)
        time.sleep(sleep / 3)
    return result
