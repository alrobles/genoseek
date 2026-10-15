"""Extract candidate species and genetic-marker evidence from text."""
from __future__ import annotations
import json
import re
from pathlib import Path


GENETIC_MARKERS = {
    "microsatellite": re.compile(r"\bmicrosatellite", re.IGNORECASE),
    "SNP": re.compile(r"\bSNP\b|single.?nucleotide polymorphism", re.IGNORECASE),
    "RADseq": re.compile(r"\bRAD[- ]?seq|restriction[- ]site", re.IGNORECASE),
    "mtDNA": re.compile(r"\bmtDNA\b|mitochondrial DNA|cytochrome[- ]?b|cox1|COI\b|ND\d+", re.IGNORECASE),
    "WGS": re.compile(r"\bwhole.?genome|low.?coverage|genome.?resequencing", re.IGNORECASE),
    "STR": re.compile(r"\bshort tandem repeat|STR\b", re.IGNORECASE),
    "allozyme": re.compile(r"\ballozyme", re.IGNORECASE),
}

_INVALID_SPECIES = {
    "sp", "spp", "spp.", "sp.", "species", "with", "and", "or", "from", "for",
    "in", "the", "a", "an", "of", "to", "by", "on", "at", "is", "are", "was",
    "were", "be", "been", "have", "has", "had", "do", "does", "did", "will",
    "would", "could", "should", "may", "might", "can", "must", "shall", "this",
    "that", "these", "those", "such", "some", "any", "all", "each", "every",
    "other", "another", "many", "much", "most", "more", "less", "few", "several",
    "both", "either", "neither", "one", "two", "three", "gen", "nov", "sensu",
    "lato", "stricto", "complex", "group", "clade", "lineage", "population",
    "populations", "individuals", "specimens", "samples", "sample", "about",
    "around", "among", "between", "within", "without", "across", "along",
    "against", "into", "onto", "upon", "over", "under", "above", "below",
    "near", "during", "before", "after", "while", "since", "until", "because",
    "although", "though", "if", "unless", "whether", "once", "whenever",
    "wherever", "however", "therefore", "thus", "hence", "meanwhile",
    "furthermore", "moreover", "nevertheless", "nonetheless", "otherwise",
    "instead", "accordingly", "consequently", "besides", "indeed", "likewise",
    "similarly", "et", "al", "van", "von", "de", "la", "da", "et", "not",
    "salamanders", "as", "mitochondria", "hybrid", "chromosomes", "genome", "dna",
    "rna", "protein", "gene", "genes", "sequence", "sequences", "clone", "clones",
    "subsection", "subgenus", "section", "subsp", "ssp", "ssp.", "var", "var.",
    "f", "f.", "forma", "phylogeny", "phylogenetic", "morphology", "taxonomy",
    "systematics", "biogeography", "evolution", "genetics", "genomics", "biology",
    "ecology", "ecological", "molecular", "genetic", "genomic", "geographic",
    "geographical", "distribution", "range", "habitat", "diversity", "conservation",
    "history", "patterns", "processes", "analysis", "analyses", "data", "study",
    "studies", "using", "based", "structure", "differentiation", "flow",
    "divergence", "migration", "selection", "adaptation", "demography", "expansion",
    "contraction", "refugium", "refugia", "bottleneck", "fragmentation", "isolation",
    "connectivity", "dispersal", "colonization", "recolonization", "introgression",
    "hybridization", "speciation", "diversification", "lineages", "subspecies",
    "morphological", "anatomical", "physiological", "behavioral", "ecological",
    "geological", "climatic", "environmental", "seasonal", "annual", "perennial",
    "terrestrial", "aquatic", "marine", "freshwater", "coastal", "montane",
    "alpine", "boreal", "tropical", "subtropical", "temperate", "arid",
    "continental", "insular", "endemic", "native", "introduced", "invasive",
    "naturalized", "cultivated", "domesticated", "wild", "captive", "feral",
    "herbivorous", "carnivorous", "omnivorous", "insectivorous", "granivorous",
    "frugivorous", "folivorous", "nectarivorous", "piscivorous",
}


POPULATION_TERMS = re.compile(
    r"\b(?:population genetics|genetic structure|gene flow|phylogeography|"
    r"F[_\-]?ST|G[_\-]?ST|differentiation|refug|bottleneck|dispersal)\b",
    re.IGNORECASE,
)


def load_genera(target_groups_path: str):
    with open(target_groups_path, encoding="utf-8") as f:
        groups = json.load(f)
    genera = set()
    for g in groups.values():
        genera.update(g.get("genera", []))
    return genera


_BINOMIAL_RE = re.compile(r"\b([A-Z][a-z]+)\s+([a-z]{2,})\b")

def extract_binomials(text: str, allowed_genera: set[str] | None = None):
    """Return a set of binomial names found in text."""
    if not text:
        return set()
    found = set()
    for genus, species in _BINOMIAL_RE.findall(text):
        if allowed_genera and genus not in allowed_genera:
            continue
        # Exclude common false positives
        if genus in {"Fig", "Table", "Supplementary", "Electronic", "This", "Our", "Download"}:
            continue
        if species.lower() in _INVALID_SPECIES or species.endswith((".", ",", ";")):
            continue
        if re.match(r"^\d", species):
            continue
        found.add(f"{genus} {species}")
    return found


def detect_markers(text: str):
    """Return list of marker categories detected in text."""
    return sorted([m for m, pat in GENETIC_MARKERS.items() if pat.search(text)])


def evidence_score(text: str, markers: list[str]):
    """Simple relevance score based on population terms and markers."""
    if not text:
        return 0.0
    pop_hits = len(POPULATION_TERMS.findall(text))
    marker_hits = len(markers)
    return pop_hits + 2.0 * marker_hits
