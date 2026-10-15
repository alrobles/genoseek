#!/usr/bin/env python3
"""Combine genominer candidate, genetic-ingest and reference-genome reports
into an actionable SPA-ready candidate table.
"""
from __future__ import annotations
import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.spa_scheduler import candidate_priority


def _vcf_type_to_strategy(vcf_type: str) -> str:
    """Map a VCF type directory name to a library strategy string."""
    if vcf_type.lower() == "stacks":
        return "GBS"
    return "WGS"


def _identifier(row: dict) -> str:
    """Best stable sample identifier for matching SRA and VCF sources."""
    return (row.get("srr") or row.get("biosample") or row.get("sample") or "").strip()


def _load_sra_records(sample_csv: Path) -> list[dict]:
    if not sample_csv.is_file():
        return []
    with open(sample_csv, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _load_vcf_popmap_records(vcf_db: Path, dataset: str) -> list[dict]:
    """Read popmap.tsv files from the VCF database.

    Returns rows with 'sample', 'srr', 'biosample', 'lat', 'lon',
    'library_strategy' and 'population' keys so they can be merged with SRA.
    """
    vcf_dir = vcf_db / "files" / dataset
    if not vcf_dir.is_dir():
        return []
    records: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for vcf_type_dir in sorted(vcf_dir.iterdir()):
        if not vcf_type_dir.is_dir():
            continue
        popmap = vcf_type_dir / "popmap.tsv"
        if not popmap.is_file():
            continue
        strategy = _vcf_type_to_strategy(vcf_type_dir.name)
        with open(popmap, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                raw_sample = (row.get("sample") or "").strip()
                if not raw_sample:
                    continue
                biosample = ""
                srr = ""
                if raw_sample.upper().startswith(("SAMN", "SAME", "SAMD", "SAMS")):
                    biosample = raw_sample
                elif raw_sample.upper().startswith(("SRR", "ERR", "DRR")):
                    srr = raw_sample
                else:
                    srr = raw_sample
                if (raw_sample, strategy) in seen:
                    continue
                seen.add((raw_sample, strategy))
                records.append({
                    "sample": raw_sample,
                    "srr": srr,
                    "biosample": biosample,
                    "population": (row.get("population") or "").strip(),
                    "library_strategy": strategy,
                    "lat": (row.get("latitude") or row.get("lat") or "").strip(),
                    "lon": (row.get("longitude") or row.get("lon") or "").strip(),
                })
    return records


def _merge_records(sra_records: list[dict], vcf_records: list[dict]) -> list[dict]:
    """Merge SRA metadata with VCF popmap coordinates.

    VCF popmaps often have more complete coordinates for WGS samples than the
    SRA CSV. If the same sample appears in both, we keep the SRA library
    strategy and the VCF lat/lon (when SRA lat/lon is missing) without double
    counting.
    """
    if not vcf_records:
        return sra_records

    # Indexes for matching VCF samples back to SRA rows.
    by_srr: dict[str, dict] = {}
    by_biosample: dict[str, dict] = {}
    by_sample: dict[str, dict] = {}
    for r in sra_records:
        srr = (r.get("srr") or "").strip()
        biosample = (r.get("biosample") or "").strip()
        sample = (r.get("sample") or "").strip()
        if srr:
            by_srr[srr] = r
        if biosample:
            by_biosample[biosample] = r
        if sample:
            by_sample[sample] = r

    merged: dict[tuple[str, str], dict] = {}
    for r in sra_records:
        rid = _identifier(r)
        strat = (r.get("library_strategy") or "").strip()
        if not rid or not strat:
            continue
        merged[(rid, strat)] = dict(r)

    for r in vcf_records:
        raw = r["sample"]
        sra = by_srr.get(raw) or by_biosample.get(raw) or by_sample.get(raw)
        if sra:
            rid = _identifier(sra)
            # Prefer SRA library_strategy, fall back to VCF inference.
            strat = (sra.get("library_strategy") or r["library_strategy"]).strip()
        else:
            rid = raw
            strat = r["library_strategy"]
        if not rid or not strat:
            continue
        key = (rid, strat)
        if key in merged:
            existing = merged[key]
            # Prefer VCF coordinates when the SRA record lacks them.
            if not (existing.get("lat") and existing.get("lon")):
                if r.get("lat") and r.get("lon"):
                    existing["lat"] = r["lat"]
                    existing["lon"] = r["lon"]
            if not existing.get("population") and r.get("population"):
                existing["population"] = r["population"]
        else:
            merged[key] = dict(r)
            merged[key]["srr"] = sra.get("srr") if sra else r.get("srr", "")
            merged[key]["biosample"] = sra.get("biosample") if sra else r.get("biosample", "")
            merged[key]["sample"] = rid
    return list(merged.values())


def _load_sample_records(genetic_dir: str, dataset: str, vcf_db: Path | None) -> list[dict]:
    """Load SRA sample rows and merge with VCF popmap rows when available."""
    sample_csv = Path(genetic_dir) / dataset / f"{dataset}_sample_coordinates.csv"
    sra_records = _load_sra_records(sample_csv)
    vcf_records = _load_vcf_popmap_records(vcf_db, dataset) if vcf_db else []
    if not sra_records:
        return vcf_records
    if not vcf_records:
        return sra_records
    return _merge_records(sra_records, vcf_records)


def _load_n_populations(genetic_dir: str, dataset: str, vcf_db: Path | None) -> int:
    """Return population count from population_coordinates.csv or VCF popmaps."""
    pop_csv = Path(genetic_dir) / dataset / f"{dataset}_population_coordinates.csv"
    if pop_csv.is_file():
        with open(pop_csv, newline="", encoding="utf-8") as f:
            return max(0, sum(1 for _ in f) - 1)
    if vcf_db:
        vcf_dir = vcf_db / "files" / dataset
        if vcf_dir.is_dir():
            populations: set[str] = set()
            for vcf_type_dir in vcf_dir.iterdir():
                if not vcf_type_dir.is_dir():
                    continue
                popmap = vcf_type_dir / "popmap.tsv"
                if not popmap.is_file():
                    continue
                with open(popmap, newline="", encoding="utf-8") as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        pop = (row.get("population") or "").strip()
                        if pop:
                            populations.add(pop)
            return len(populations)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", default="data/candidates.tsv")
    parser.add_argument("--genetic", default="data/genetic_ingest_report.tsv")
    parser.add_argument("--references", default="data/reference_genome_report.tsv")
    parser.add_argument("--genetic-dir", default="data/genetic")
    parser.add_argument("--vcf-database", default="data/vcf_database",
                        help="Root of the VCF database (used as fallback for coordinates)")
    parser.add_argument("--out", default="data/spa_ready_candidates.tsv")
    args = parser.parse_args()

    vcf_db = Path(args.vcf_database) if args.vcf_database else None

    candidates: dict[str, dict] = {}
    with open(args.candidates, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            candidates[row["binomial"]] = row

    genetic: dict[str, dict] = {}
    with open(args.genetic, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genetic[row["binomial"]] = row

    references: dict[str, dict] = {}
    with open(args.references, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            references[row["binomial"]] = row

    out_rows = []
    for binomial in sorted(candidates):
        cand = candidates[binomial]
        gen = genetic.get(binomial, {})
        ref = references.get(binomial, {})
        dataset = binomial.replace(" ", "_").lower()

        # Count samples with coordinates by library strategy.
        # Prefer data/genetic sample_coordinates.csv; fall back to VCF popmaps
        # so species processed by the WGS pipeline are also counted.
        sample_records = _load_sample_records(args.genetic_dir, dataset, vcf_db)
        strategy_counts: dict[str, dict[str, int]] = {}
        for row in sample_records:
            strat = (row.get("library_strategy") or "").strip()
            has = bool(row.get("lat") and row.get("lon"))
            if strat not in strategy_counts:
                strategy_counts[strat] = {"with_coords": 0, "total": 0}
            strategy_counts[strat]["total"] += 1
            if has:
                strategy_counts[strat]["with_coords"] += 1

        wgs_coords = sum(v["with_coords"] for k, v in strategy_counts.items() if "WGS" in k.upper())
        wgs_total = sum(v["total"] for k, v in strategy_counts.items() if "WGS" in k.upper())
        radseq_coords = sum(v["with_coords"] for k, v in strategy_counts.items() if k.upper() in {"RAD-SEQ", "GBS"})
        capture_coords = sum(v["with_coords"] for k, v in strategy_counts.items() if "TARGETED" in k.upper() or "CAPTURE" in k.upper())

        # Use the CSV/VCF popmaps as a fallback/override in case the
        # genetic_ingest_report is stale, partial, or was produced before the
        # WGS pipeline populated the coordinate files.
        csv_with_coords = sum(v["with_coords"] for v in strategy_counts.values())
        csv_n_pops = _load_n_populations(args.genetic_dir, dataset, vcf_db)

        ref_acc = ref.get("ref_genome_accession", "")
        ref_level = ref.get("ref_assembly_level", "")
        ref_size = ref.get("ref_total_sequence_length", "")
        has_ref = bool(ref_acc and ref_level)

        n_coords = max(
            int(gen.get("n_with_coords", "0") or 0),
            csv_with_coords,
        )
        n_pops = max(
            int(gen.get("n_populations", "0") or 0),
            csv_n_pops,
        )
        gbif_occ = int(cand.get("gbif_occurrence_count", "0") or 0)

        status = []
        if n_coords >= 3:
            status.append("coords")
        if gbif_occ >= 100:
            status.append("gbif")
        if has_ref:
            status.append("ref_genome")
        if wgs_coords >= 3:
            status.append("wgs_coords")
        if radseq_coords + capture_coords >= 3:
            status.append("reduced_rep_coords")

        ready_spa = ("coords" in status and "gbif" in status and ("wgs_coords" in status or "reduced_rep_coords" in status))
        ready_wgs_gst = ("wgs_coords" in status and "ref_genome" in status)

        priority = candidate_priority(cand)

        out_rows.append({
            "binomial": binomial,
            "dataset": dataset,
            "priority": round(priority, 2),
            "gbif_occurrence_count": gbif_occ,
            "n_sra_runs": gen.get("n_sra_runs", "0"),
            "n_with_coords": n_coords,
            "n_populations": n_pops,
            "wgs_with_coords": wgs_coords,
            "wgs_total": wgs_total,
            "radseq_or_gbs_with_coords": radseq_coords,
            "targeted_capture_with_coords": capture_coords,
            "ref_genome_accession": ref_acc,
            "ref_organism_name": ref.get("ref_organism_name", ""),
            "ref_assembly_level": ref_level,
            "ref_refseq_category": ref.get("ref_refseq_category", ""),
            "ref_total_sequence_length": ref_size,
            "ref_ftp_genomic_fasta": ref.get("ref_ftp_genomic_fasta", ""),
            "ready_for_spa": int(ready_spa),
            "ready_for_wgs_gst": int(ready_wgs_gst),
            "status_flags": ",".join(status) if status else "none",
        })

    out_rows.sort(key=lambda x: (x["ready_for_wgs_gst"], x["ready_for_spa"], x["priority"]), reverse=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "binomial", "dataset", "priority", "gbif_occurrence_count", "n_sra_runs",
            "n_with_coords", "n_populations", "wgs_with_coords", "wgs_total",
            "radseq_or_gbs_with_coords", "targeted_capture_with_coords",
            "ref_genome_accession", "ref_organism_name", "ref_assembly_level",
            "ref_refseq_category", "ref_total_sequence_length", "ref_ftp_genomic_fasta",
            "ready_for_spa", "ready_for_wgs_gst", "status_flags"
        ], delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {args.out} with {len(out_rows)} candidates")


if __name__ == "__main__":
    main()
