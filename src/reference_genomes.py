"""Find and rank reference genomes for a binomial using the NCBI Datasets API.
"""
from __future__ import annotations
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Optional

DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon"


def _safe_int(value) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _ftp_url(accession: str, assembly_name: str) -> str:
    """Return the canonical NCBI FTP URL for the genomic FASTA.

    Path convention:
      intermediate dirs use accession number without version, grouped by 3 digits.
      final directory and filename use full accession including version.
    """
    prefix = accession.split("_")[0]
    acc_no_version = accession.split(".")[0]
    digits = re.sub(r"[^0-9]", "", acc_no_version)
    parts = "/".join(digits[i : i + 3] for i in range(0, len(digits), 3))
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{parts}/"
        f"{accession}_{assembly_name}/{accession}_{assembly_name}_genomic.fna.gz"
    )


def _head_status(url: str, timeout: int = 20) -> int:
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status
    except Exception:
        return 0


def _parse_report(report: dict) -> dict:
    org = report.get("organism", {})
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    acc = report.get("accession") or report.get("current_accession", "")
    asm_name = info.get("assembly_name", "")
    ftp = _ftp_url(acc, asm_name) if acc and asm_name else ""
    return {
        "accession": acc,
        "assembly_name": asm_name,
        "organism_name": org.get("organism_name", ""),
        "tax_id": org.get("tax_id"),
        "common_name": org.get("common_name", ""),
        "source_database": report.get("source_database", ""),
        "assembly_level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category", ""),
        "release_date": info.get("release_date", ""),
        "bioproject_accession": info.get("bioproject_accession", ""),
        "ftp_genomic_fasta": ftp,
        "total_sequence_length": _safe_int(stats.get("total_sequence_length")),
        "total_ungapped_length": _safe_int(stats.get("total_ungapped_length")),
        "number_of_scaffolds": _safe_int(stats.get("number_of_scaffolds")),
        "scaffold_n50": _safe_int(stats.get("scaffold_n50")),
        "number_of_chromosomes": _safe_int(stats.get("total_number_of_chromosomes")),
        "gc_percent": stats.get("gc_percent"),
        "genome_coverage": stats.get("genome_coverage", ""),
    }


def _score(report: dict, target_binomial: str) -> int:
    score = 0
    target_genus = target_binomial.split()[0]
    org = report.get("organism_name", "")
    org_parts = org.split()
    if org == target_binomial:
        score += 100
    elif org_parts and org_parts[0] == target_genus:
        score += 50
    level = report.get("assembly_level", "").lower()
    if level == "chromosome":
        score += 40
    elif level == "scaffold":
        score += 25
    elif level == "contig":
        score += 10
    if report.get("refseq_category", "").lower() == "reference genome":
        score += 30
    return score


def fetch_assemblies(
    taxon: str,
    page_size: int = 20,
    retries: int = 3,
    sleep: float = 0.5,
) -> list[dict]:
    """Query NCBI Datasets for genome assembly reports for a taxon name."""
    name = urllib.parse.quote(taxon)
    url = f"{DATASETS_BASE}/{name}/dataset_report?filters.assembly_source=all&page_size={page_size}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode("utf-8"))
            return [_parse_report(rep) for rep in data.get("reports", [])]
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            if attempt < retries - 1:
                time.sleep(sleep * (attempt + 1))
                continue
            return []
        except Exception:
            if attempt < retries - 1:
                time.sleep(sleep * (attempt + 1))
                continue
            return []
    return []


def select_reference(binomial: str, sleep: float = 0.5) -> Optional[dict]:
    """Return the best reference genome report for a binomial.

    Searches the species name first, then the genus.
    """
    import json  # local import to avoid top-level if not needed
    parts = binomial.split()
    if len(parts) < 2:
        return None
    species_name = binomial
    genus_name = parts[0]

    reports = fetch_assemblies(species_name, page_size=20, sleep=sleep)
    time.sleep(sleep)
    if not reports:
        reports = fetch_assemblies(genus_name, page_size=20, sleep=sleep)
        time.sleep(sleep)
    if not reports:
        return None

    for rep in reports:
        rep["score"] = _score(rep, binomial)
    reports.sort(key=lambda x: x["score"], reverse=True)
    best = reports[0]
    # Verify FASTA URL is reachable when it is set.
    if best.get("ftp_genomic_fasta"):
        if _head_status(best["ftp_genomic_fasta"], timeout=10) != 200:
            best["ftp_genomic_fasta"] = ""
    best["n_alternative_assemblies"] = len(reports) - 1
    return best
