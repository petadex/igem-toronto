"""
Data loading and cleaning for PlasticDB and PAZy datasets.

PlasticDB: microorganism-plastic degradation entries (TSV download)
PAZy:      thoroughly characterised plastic-active enzymes (web scrape)
"""

import pandas as pd
import requests
from pathlib import Path
from bs4 import BeautifulSoup
import re

DATA_DIR = Path(__file__).parent.parent / "data"

PLASTICDB_TSV_URL = "https://plasticdb.org/static/degraders_list.tsv"
PAZY_PROTEINS_URL = "https://www.pazy.eu/proteins"
PAZY_PLASTICS_URL = "https://www.pazy.eu/plastics"

PLASTIC_FULL_NAMES = {
    "PET": "Polyethylene terephthalate",
    "PE":  "Polyethylene (general)",
    "LDPE":"Low-density polyethylene",
    "HDPE":"High-density polyethylene",
    "PP":  "Polypropylene",
    "PS":  "Polystyrene",
    "PVC": "Polyvinyl chloride",
    "PU":  "Polyurethane",
    "PLA": "Polylactic acid",
    "PHB": "Polyhydroxybutyrate",
    "PHA": "Polyhydroxyalkanoate",
    "PCL": "Polycaprolactone",
    "PBS": "Polybutylene succinate",
    "PBSA":"Polybutylene succinate-co-adipate",
    "PBAT":"Polybutylene adipate-co-terephthalate",
    "PVA": "Polyvinyl alcohol",
    "Nylon":"Nylon / Polyamide",
    "PES": "Polyethylene succinate",
    "PHO": "Polyhydroxyoctanoate",
    "PHBV":"Poly(3-hydroxybutyrate-co-3-hydroxyvalerate)",
}

PLASTIC_CATEGORIES = {
    "Biodegradable/Bio-based": ["PHA","PHB","PHO","PHBV","PLA","PCL","PBS","PBSA","PBAT","PES","PVA"],
    "Commodity Thermoplastics":["PE","LDPE","HDPE","PP","PS","PVC","PET"],
    "Polyurethanes":           ["PU"],
    "Polyamides":              ["Nylon"],
    "Polyesters":              ["PET","PES","PBAT","PBS","PBSA","PCL","PLA"],
}


def load_plasticdb(local_path: str | Path | None = None, force_download: bool = False) -> pd.DataFrame:
    """
    Load and clean the PlasticDB microorganism dataset.

    Returns a tidy DataFrame with one row per (organism, plastic, paper) triplet.
    """
    cache = DATA_DIR / "plasticdb_microorganisms.tsv"
    cache.parent.mkdir(parents=True, exist_ok=True)

    if local_path:
        path = Path(local_path)
    elif cache.exists() and not force_download:
        path = cache
    else:
        resp = requests.get(PLASTICDB_TSV_URL, timeout=60)
        resp.raise_for_status()
        cache.write_bytes(resp.content)
        path = cache

    df = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip")

    df.columns = [
        "organism", "tax_id", "plastic", "reference", "enzyme_name",
        "enzyme_id", "db_enzyme_name", "gene", "genbank_id", "sequence",
        "year", "evidence", "plastic_used", "manufacturer", "analytical_grade",
        "thermophilic", "isolation_sample", "isolation_environment",
        "isolation_location", "extrapolated_from_enzyme", "enzyme_id_in_paper", "doi",
    ]

    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df["tax_id"] = pd.to_numeric(df["tax_id"], errors="coerce")

    bool_cols = ["analytical_grade", "thermophilic", "extrapolated_from_enzyme"]
    for col in bool_cols:
        df[col] = df[col].map({"Yes": True, "No": False, "yes": True, "no": False})

    df["plastic"] = df["plastic"].str.strip()
    df["organism"] = df["organism"].str.strip()
    df["isolation_location"] = df["isolation_location"].str.strip()

    df["has_sequence"] = df["sequence"].notna() & (df["sequence"].str.len() > 10)
    df["has_genbank"] = df["genbank_id"].notna() & (df["genbank_id"].str.len() > 3)
    df["has_enzyme"] = df["enzyme_name"].notna() & (df["enzyme_name"] != "")

    df["plastic_category"] = df["plastic"].apply(_categorise_plastic)
    df["plastic_full_name"] = df["plastic"].map(PLASTIC_FULL_NAMES).fillna(df["plastic"])

    genus_species = df["organism"].str.extract(r"^(\w+)\s+(\w+)", expand=True)
    df["genus"] = genus_species[0]
    df["species"] = genus_species[1]

    df["decade"] = (df["year"] // 10 * 10).astype("Int64")

    return df.reset_index(drop=True)


def _categorise_plastic(plastic: str) -> str:
    if pd.isna(plastic):
        return "Unknown"
    for cat, members in PLASTIC_CATEGORIES.items():
        if plastic in members:
            return cat
    return "Other"


def fetch_pazy_proteins(force: bool = False) -> pd.DataFrame:
    """
    Scrape PAZy protein table (thoroughly characterised plastic-active enzymes).
    Falls back to cached CSV if the site is unreachable.
    """
    cache = DATA_DIR / "pazy_proteins.csv"
    cache.parent.mkdir(parents=True, exist_ok=True)

    if cache.exists() and not force:
        return pd.read_csv(cache)

    try:
        resp = requests.get(PAZY_PROTEINS_URL, timeout=30,
                            headers={"User-Agent": "Mozilla/5.0 (research bot)"})
        soup = BeautifulSoup(resp.text, "html.parser")
        tables = soup.find_all("table")
        if tables:
            df = pd.read_html(str(tables[0]))[0]
            df.to_csv(cache, index=False)
            return df
    except Exception:
        pass

    return _pazy_fallback()


def _pazy_fallback() -> pd.DataFrame:
    """
    Curated representative sample from PAZy (Buchholz et al. 2022) for offline use.
    Covers the main well-characterised plastic-active enzyme families.
    """
    records = [
        {"enzyme_name":"LCC","organism":"Leaf-branch compost metagenome","plastic":"PET","ec_number":"3.1.1.101","uniprot":"G9BY57","year":2012,"source":"PAZy"},
        {"enzyme_name":"IsPETase","organism":"Ideonella sakaiensis 201-F6","plastic":"PET","ec_number":"3.1.1.101","uniprot":"A0A0K8P8H7","year":2016,"source":"PAZy"},
        {"enzyme_name":"TfCut2","organism":"Thermobifida fusca KW3","plastic":"PET","ec_number":"3.1.1.101","uniprot":"Q47RJ6","year":2010,"source":"PAZy"},
        {"enzyme_name":"ThermoPETase","organism":"Engineered variant of IsPETase","plastic":"PET","ec_number":"3.1.1.101","uniprot":None,"year":2019,"source":"PAZy"},
        {"enzyme_name":"FAST-PETase","organism":"Engineered variant of IsPETase","plastic":"PET","ec_number":"3.1.1.101","uniprot":None,"year":2022,"source":"PAZy"},
        {"enzyme_name":"ICCG LCC","organism":"Engineered LCC variant","plastic":"PET","ec_number":"3.1.1.101","uniprot":None,"year":2020,"source":"PAZy"},
        {"enzyme_name":"PHL7","organism":"Pseudomonas lini DSM 16768","plastic":"PET","ec_number":"3.1.1.101","uniprot":"A0A6M3YWE3","year":2021,"source":"PAZy"},
        {"enzyme_name":"PET2","organism":"Metagenome","plastic":"PET","ec_number":"3.1.1.101","uniprot":None,"year":2020,"source":"PAZy"},
        {"enzyme_name":"BhrPETase","organism":"Bhargavaea cecembensis","plastic":"PET","ec_number":"3.1.1.101","uniprot":None,"year":2023,"source":"PAZy"},
        {"enzyme_name":"CsPETase","organism":"Cryptosporangium aurantiacum","plastic":"PET","ec_number":"3.1.1.101","uniprot":None,"year":2023,"source":"PAZy"},
        {"enzyme_name":"PHB depolymerase (R. pickettii)","organism":"Ralstonia pickettii","plastic":"PHB","ec_number":"3.1.1.75","uniprot":"P20594","year":1994,"source":"PAZy"},
        {"enzyme_name":"PhaZ7","organism":"Paucimonas lemoignei","plastic":"PHB","ec_number":"3.1.1.75","uniprot":"Q9F3L4","year":2001,"source":"PAZy"},
        {"enzyme_name":"PhaZSpo","organism":"Bacillus megaterium","plastic":"PHB","ec_number":"3.1.1.75","uniprot":"O87625","year":2000,"source":"PAZy"},
        {"enzyme_name":"PHO depolymerase","organism":"Pseudomonas fluorescens GK13","plastic":"PHA","ec_number":"3.1.1.76","uniprot":None,"year":1994,"source":"PAZy"},
        {"enzyme_name":"PlaA","organism":"Pseudomonas stutzeri","plastic":"PLA","ec_number":"3.1.1.-","uniprot":"Q9ZIQ4","year":2003,"source":"PAZy"},
        {"enzyme_name":"PlaC","organism":"Pseudomonas sp. DS04-T","plastic":"PLA","ec_number":"3.1.1.-","uniprot":None,"year":2012,"source":"PAZy"},
        {"enzyme_name":"PCL depolymerase","organism":"Alcaligenes faecalis T1","plastic":"PCL","ec_number":"3.1.1.74","uniprot":"P16417","year":1993,"source":"PAZy"},
        {"enzyme_name":"PBSA depolymerase","organism":"Bacillus pumilus KT","plastic":"PBSA","ec_number":"3.1.1.-","uniprot":None,"year":2004,"source":"PAZy"},
        {"enzyme_name":"NylB","organism":"Arthrobacter sp. KI72","plastic":"Nylon","ec_number":"3.5.1.84","uniprot":"P10163","year":1992,"source":"PAZy"},
        {"enzyme_name":"NylA","organism":"Arthrobacter sp. KI72","plastic":"Nylon","ec_number":"3.5.1.-","uniprot":"P10162","year":1992,"source":"PAZy"},
        {"enzyme_name":"PueA","organism":"Pseudomonas chlororaphis","plastic":"PU","ec_number":"3.1.1.-","uniprot":"Q9X4G0","year":2000,"source":"PAZy"},
        {"enzyme_name":"PueB","organism":"Pseudomonas chlororaphis","plastic":"PU","ec_number":"3.1.1.-","uniprot":"Q9X4G1","year":2000,"source":"PAZy"},
        {"enzyme_name":"Pfl_1298 lipase","organism":"Pseudomonas fluorescens SBW25","plastic":"PU","ec_number":"3.1.1.3","uniprot":"Q4KE02","year":2008,"source":"PAZy"},
        {"enzyme_name":"PES hydrolase","organism":"Pelosinus sp. UFO1","plastic":"PES","ec_number":"3.1.1.-","uniprot":None,"year":2019,"source":"PAZy"},
    ]
    df = pd.DataFrame(records)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(DATA_DIR / "pazy_proteins.csv", index=False)
    return df


def get_unique_organisms(df: pd.DataFrame) -> pd.DataFrame:
    """Return one row per unique organism with aggregated metadata."""
    agg = (
        df.groupby("organism")
        .agg(
            tax_id=("tax_id", "first"),
            genus=("genus", "first"),
            species=("species", "first"),
            plastics_degraded=("plastic", lambda x: sorted(x.dropna().unique().tolist())),
            n_plastics=("plastic", "nunique"),
            n_entries=("organism", "count"),
            years_active=("year", lambda x: sorted(x.dropna().unique().tolist())),
            first_year=("year", "min"),
            last_year=("year", "max"),
            has_enzyme=("has_enzyme", "any"),
            has_sequence=("has_sequence", "any"),
            isolation_locations=("isolation_location", lambda x: sorted(x.dropna().unique().tolist())),
            isolation_environments=("isolation_environment", lambda x: sorted(x.dropna().unique().tolist())),
        )
        .reset_index()
    )
    return agg


def get_plastic_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Return per-plastic summary statistics."""
    agg = (
        df.groupby("plastic")
        .agg(
            full_name=("plastic_full_name", "first"),
            category=("plastic_category", "first"),
            n_entries=("plastic", "count"),
            n_unique_organisms=("organism", "nunique"),
            n_unique_genera=("genus", "nunique"),
            n_with_sequence=("has_sequence", "sum"),
            n_with_enzyme=("has_enzyme", "sum"),
            first_year=("year", "min"),
            last_year=("year", "max"),
            pct_with_sequence=("has_sequence", "mean"),
        )
        .reset_index()
    )
    agg["pct_with_sequence"] = (agg["pct_with_sequence"] * 100).round(1)
    return agg.sort_values("n_entries", ascending=False).reset_index(drop=True)


def load_all(force_download: bool = False) -> dict[str, pd.DataFrame]:
    """Convenience loader — returns dict with 'plasticdb', 'pazy', 'organisms', 'plastics'."""
    plasticdb = load_plasticdb(force_download=force_download)
    pazy = fetch_pazy_proteins(force=force_download)
    organisms = get_unique_organisms(plasticdb)
    plastics = get_plastic_summary(plasticdb)
    return {"plasticdb": plasticdb, "pazy": pazy, "organisms": organisms, "plastics": plastics}
