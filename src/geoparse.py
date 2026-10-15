"""Extract geographic place names from text using a curated gazetteer."""
from __future__ import annotations
import re
from pathlib import Path

# Curated list of places relevant to temperate refugia / landscape genetics.
# Can be replaced by a full gazetteer in production.
_GAZETTEER = {
    # Continents / large regions
    "North America", "South America", "Europe", "Asia", "Africa",
    "North America", "South America", "Eurasia", "Holarctic", "Neartic",
    "Mediterranean", "boreal", "temperate", "arctic", "alpine",
    "Nearctic", "Palearctic", "Afrotropical", "Neotropical",

    # Countries (selection)
    "Canada", "United States", "USA", "Mexico", "Guatemala", "Costa Rica",
    "Panama", "Colombia", "Ecuador", "Peru", "Bolivia", "Chile", "Argentina",
    "United Kingdom", "Ireland", "France", "Spain", "Portugal", "Italy",
    "Greece", "Turkey", "Germany", "Poland", "Sweden", "Norway", "Finland",
    "Denmark", "Netherlands", "Belgium", "Switzerland", "Austria", "Czech Republic",
    "Slovakia", "Hungary", "Romania", "Bulgaria", "Serbia", "Croatia", "Bosnia",
    "Ukraine", "Belarus", "Lithuania", "Latvia", "Estonia", "Russia",
    "China", "Japan", "Korea", "India", "Iran", "Iraq", "Syria", "Israel",
    "Morocco", "Algeria", "Tunisia", "Libya", "Egypt", "Ethiopia", "Kenya",
    "Tanzania", "South Africa", "Madagascar", "Australia", "New Zealand",

    # US states / regions
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
    "Connecticut", "Delaware", "Florida", "Georgia", "Idaho", "Illinois",
    "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland",
    "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana",
    "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York",
    "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania",
    "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah",
    "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming",
    "New England", "Midwest", "Great Plains", "Pacific Northwest", "Southeast",
    "Southwest", "Appalachia",

    # Canadian provinces
    "Alberta", "British Columbia", "Manitoba", "New Brunswick", "Newfoundland",
    "Nova Scotia", "Ontario", "Quebec", "Saskatchewan", "Yukon", "Northwest Territories",

    # Mountain ranges / ecoregions
    "Appalachian Mountains", "Appalachians", "Blue Ridge", "Great Smoky Mountains",
    "Rocky Mountains", "Rockies", "Sierra Nevada", "Cascade Range", "Cascade Mountains",
    "Coast Range", "Coast Mountains", "Siskiyou", "Klamath Mountains",
    "Sierra Madre", "Cordillera", "Pyrenees", "Alps", "Apennines",
    "Carpathians", "Balkans", "Caucasus", "Himalayas", "Hindu Kush", "Altai",
    "Andes", "Ural", "Scandes", "Scandinavia", "Baltic", "Mediterranean Basin",
    "Iberian Peninsula", "Italian Peninsula", "Balkan Peninsula", "Anatolia",
    "Beringia", "Laurentide", "Wisconsin", "Illinoian",
}

# Build a regex that matches whole phrases from the gazetteer, longest first.
_PLACES = sorted(_GAZETTEER, key=len, reverse=True)
_PLACE_RE = re.compile(
    r"\b(?:" + "|".join(re.escape(p) for p in _PLACES) + r")\b",
    re.IGNORECASE,
)


def extract_places(text: str) -> set[str]:
    """Return the set of gazetteer place names found in text."""
    if not text:
        return set()
    return {m.group(0) for m in _PLACE_RE.finditer(text)}
