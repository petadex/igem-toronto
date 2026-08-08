"""
enhance.py — Expands organism-profiles with additional data and visualizations.

Modes:
  --mode genome        Fetch NCBI Assembly stats (genome size, GC%, assembly level)
  --mode pubmed        Fetch PubMed literature counts (organism + plastic/bioplastic)
  --mode sra-expand    Search SRA for bioplastic-specific experiments; find new organisms
  --mode pages         Regenerate all HTML pages with charts and all cached data
  --mode status        Show current cache coverage

Usage (from organism-profiles/ directory):
  python enhance.py --mode genome  --batch 80
  python enhance.py --mode pubmed  --batch 150
  python enhance.py --mode sra-expand
  python enhance.py --mode pages
"""

import argparse
import json
import re
import time
import unicodedata
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT        = Path(__file__).parent
DATA_DIR    = ROOT / "data"
PAGES_DIR   = ROOT / "pages"
PDB_TSV     = ROOT.parent / "plastic-biodegradation-analysis" / "data" / "plasticdb_microorganisms.tsv"

SRA_CSV     = DATA_DIR / "sra_stats.csv"
BACDIVE_CSV = DATA_DIR / "bacdive_data.csv"
GENOME_CSV  = DATA_DIR / "genome_data.csv"
PUBMED_CSV  = DATA_DIR / "pubmed_counts.csv"
EXTRA_CSV   = DATA_DIR / "extra_organisms.csv"   # organisms added via SRA bioplastic search

DATA_DIR.mkdir(exist_ok=True)
PAGES_DIR.mkdir(exist_ok=True)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
session   = requests.Session()
session.headers["User-Agent"] = "igem-toronto-research/1.0 (academic, bioplastic-profiles)"

# ---------------------------------------------------------------------------
# Bioplastic classification
# ---------------------------------------------------------------------------
BIOPLASTICS = {
    # PHA family
    "PHA", "PHA Blend", "PHB", "PHB-Blend", "PHBV", "PHBVH", "PHO", "PHBH",
    "PHBHHx", "PHC", "PHN", "PHV", "PHPV",
    "P3HP", "P3HV", "P4HB", "P34HB", "PMCL",
    "P(3HB-co-4HB)", "P(3HB-co-3MP)", "P(3HB-co-HV)", "P(3HV-co-4HB)",
    "P(3HB-co-3HV)", "mcl-PHA", "scl-PHA",
    # PLA family
    "PLA", "PLA Blend", "PLLA", "PDLA", "poly(lactic acid)", "polylactic acid",
    # PCL family
    "PCL", "PCL Blend",
    # PBS / PBSA / PBAT family
    "PBS", "PBS Blend", "PBS-Blend", "PBSA", "PBSA-Blend", "PBAT", "PBAT-Blend",
    "PBSeT", "PBST55", "PBSTIL", "PBSeT",
    # Other biodegradables
    "PVA", "PVA Blend", "O-PVA",
    "PEF",
    "PEA",
    "PPL",
    "PTS",
    "Ecovio-FT",
    "P(3HB-co-HV)",
    "PHBHHx",
}

CONVENTIONAL = {
    "PE", "LDPE", "HDPE", "LLDPE", "LDPE Blend", "LLDPE Blend", "PE Blend", "O-PE",
    "PET", "PETG",
    "PP",
    "PS", "PS Blend",
    "PVC", "PVC Blend",
    "PU", "Impranil", "PU Blend",
    "Nylon", "PA",
    "PC",
    "PES",
    "PEG",
    "NR",
    "PSS", "PTC",
}


def classify_plastic(p: str) -> str:
    p = (p or "").strip()
    if p in BIOPLASTICS:
        return "bioplastic"
    if p in CONVENTIONAL:
        return "conventional"
    return "other"


# ---------------------------------------------------------------------------
# Rate-limited requests
# ---------------------------------------------------------------------------
_last_req = 0.0


def ncbi_get(url: str, sleep: float = 0.4) -> requests.Response | None:
    global _last_req
    wait = sleep - (time.time() - _last_req)
    if wait > 0:
        time.sleep(wait)
    _last_req = time.time()
    try:
        r = session.get(url, timeout=20)
        r.raise_for_status()
        return r
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Load PlasticDB
# ---------------------------------------------------------------------------
def load_plasticdb():
    df = pd.read_csv(PDB_TSV, sep="\t", dtype=str, on_bad_lines="skip")
    df.columns = [
        "organism","tax_id","plastic","reference","enzyme_name",
        "enzyme_id","db_enzyme_name","gene","genbank_id","sequence",
        "year","evidence","plastic_used","manufacturer","analytical_grade",
        "thermophilic","isolation_sample","isolation_environment",
        "isolation_location","extrapolated_from_enzyme","enzyme_id_in_paper","doi",
    ]
    df["organism"]     = df["organism"].str.strip()
    df["plastic"]      = df["plastic"].str.strip()
    df["year"]         = pd.to_numeric(df["year"], errors="coerce")
    df["has_sequence"] = df["sequence"].notna() & (df["sequence"].str.len() > 10)
    df["has_enzyme"]   = df["enzyme_name"].notna() & (df["enzyme_name"] != "")
    df["has_genbank"]  = df["genbank_id"].notna() & (df["genbank_id"].str.len() > 3)
    df["plastic_class"] = df["plastic"].apply(classify_plastic)
    return df


def build_org_summary(df: pd.DataFrame) -> pd.DataFrame:
    orgs = (
        df.groupby("organism")
        .agg(
            tax_id         = ("tax_id", "first"),
            n_entries      = ("organism", "count"),
            plastics       = ("plastic", lambda x: sorted(x.dropna().unique().tolist())),
            n_plastics     = ("plastic", "nunique"),
            has_sequence   = ("has_sequence", "any"),
            has_enzyme     = ("has_enzyme", "any"),
            has_genbank    = ("has_genbank", "any"),
            first_year     = ("year", "min"),
            last_year      = ("year", "max"),
            thermophilic   = ("thermophilic", lambda x: x.mode()[0] if len(x.mode()) else ""),
            isolation_envs = ("isolation_environment",
                              lambda x: "; ".join(sorted(x.dropna().unique()))),
            isolation_locs = ("isolation_location",
                              lambda x: "; ".join(sorted(x.dropna().unique()))),
            n_bioplastic   = ("plastic_class", lambda x: (x == "bioplastic").sum()),
            n_conventional = ("plastic_class", lambda x: (x == "conventional").sum()),
        )
        .reset_index()
        .sort_values("n_entries", ascending=False)
    )
    orgs["bioplastic_relevant"] = orgs["n_bioplastic"] > 0
    return orgs


# ---------------------------------------------------------------------------
# NCBI Assembly genome data
# ---------------------------------------------------------------------------
def _parse_meta_stat(meta: str, category: str) -> str:
    """Extract a value from NCBI Assembly esummary meta XML blob."""
    m = re.search(
        rf'category="{re.escape(category)}"[^>]*>([^<]*)</Stat>',
        meta)
    return m.group(1).strip() if m else ""


def _fetch_genome_from_uid(uid: str) -> dict:
    """Fetch genome stats from a known NCBI Assembly UID."""
    url = f"{NCBI_BASE}/esummary.fcgi?db=assembly&id={uid}&retmode=json"
    r = ncbi_get(url)
    if r is None:
        return {}
    try:
        result = r.json().get("result", {})
    except Exception:
        return {}
    for k, entry in result.items():
        if k != "uids":
            return _parse_assembly_entry(entry)
    return {}


def _parse_assembly_entry(entry: dict) -> dict:
    """Extract standardised genome stats from one NCBI Assembly esummary entry."""
    meta = entry.get("meta", "")
    return {
        "genome_size_bp":     _parse_meta_stat(meta, "total_length"),
        "contig_count":       _parse_meta_stat(meta, "contig_count"),
        "contig_n50":         _parse_meta_stat(meta, "contig_n50"),
        "scaffold_n50":       _parse_meta_stat(meta, "scaffold_n50"),
        "assembly_level":     entry.get("assemblystatus", ""),
        "assembly_accession": entry.get("assemblyaccession", ""),
        "assembly_name":      entry.get("assemblyname", ""),
        "coverage":           entry.get("coverage", ""),
        "taxid":              entry.get("taxid", ""),
    }


def _fetch_genome_one(organism: str) -> dict:
    row = {"organism": organism}
    # Search assembly
    url = (f"{NCBI_BASE}/esearch.fcgi?db=assembly"
           f"&term={requests.utils.quote(organism)}[Organism]"
           f"&retmax=5&retmode=json")
    r = ncbi_get(url)
    if r is None:
        return row
    try:
        ids = r.json()["esearchresult"].get("idlist", [])
    except Exception:
        return row
    if not ids:
        return row

    url2 = f"{NCBI_BASE}/esummary.fcgi?db=assembly&id={','.join(ids[:5])}&retmode=json"
    r2 = ncbi_get(url2)
    if r2 is None:
        return row
    try:
        result = r2.json().get("result", {})
    except Exception:
        return row

    best = None
    level_order = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    for uid, entry in result.items():
        if uid == "uids":
            continue
        level = entry.get("assemblystatus", "")
        if best is None or level_order.get(level, 9) < level_order.get(best.get("assemblystatus",""), 9):
            best = entry

    if best:
        parsed = _parse_assembly_entry(best)
        row.update(parsed)

    return row


def run_genome_batch(orgs: pd.DataFrame, batch: int):
    cached = pd.read_csv(GENOME_CSV, dtype=str) if GENOME_CSV.exists() else pd.DataFrame()
    done   = set(cached["organism"].tolist()) if not cached.empty else set()
    todo   = orgs[~orgs["organism"].isin(done)].head(batch)

    if todo.empty:
        print(f"Genome: all {len(cached)} organisms already cached.")
        return

    print(f"Genome: fetching {len(todo)} (cache has {len(done)})...")
    rows = []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        rows.append(_fetch_genome_one(row["organism"]))
        if i % 20 == 0:
            # Checkpoint: save every 20 so timeout doesn't lose progress
            ckpt = pd.concat([cached, pd.DataFrame(rows)], ignore_index=True) if not cached.empty else pd.DataFrame(rows)
            ckpt.to_csv(GENOME_CSV, index=False)
            print(f"  {i}/{len(todo)} (checkpoint saved)")

    new_df   = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(GENOME_CSV, index=False)
    found = int(combined["assembly_accession"].notna().sum()) if "assembly_accession" in combined.columns else 0
    print(f"Genome: saved {len(combined)} rows. With assembly: {found}.")


# ---------------------------------------------------------------------------
# PubMed literature count
# ---------------------------------------------------------------------------
def _fetch_pubmed_one(organism: str) -> dict:
    row = {"organism": organism}

    # Plastic-specific count
    term = (f'"{organism}"[Organism] AND '
            f'(plastic* OR bioplastic* OR polymer* OR polyethylene OR polypropylene '
            f'OR polylactic OR polyhydroxy OR polycaprolactone OR degradation)')
    url = (f"{NCBI_BASE}/esearch.fcgi?db=pubmed"
           f"&term={requests.utils.quote(term)}&retmax=0&retmode=json")
    r = ncbi_get(url, sleep=0.35)
    if r:
        try:
            row["pubmed_plastic_count"] = int(r.json()["esearchresult"].get("count", 0))
        except Exception:
            pass

    # Total publication count
    term2 = f'"{organism}"[Organism]'
    url2 = (f"{NCBI_BASE}/esearch.fcgi?db=pubmed"
            f"&term={requests.utils.quote(term2)}&retmax=0&retmode=json")
    r2 = ncbi_get(url2, sleep=0.35)
    if r2:
        try:
            row["pubmed_total_count"] = int(r2.json()["esearchresult"].get("count", 0))
        except Exception:
            pass

    return row


def run_pubmed_batch(orgs: pd.DataFrame, batch: int):
    cached = pd.read_csv(PUBMED_CSV, dtype=str) if PUBMED_CSV.exists() else pd.DataFrame()
    done   = set(cached["organism"].tolist()) if not cached.empty else set()
    todo   = orgs[~orgs["organism"].isin(done)].head(batch)

    if todo.empty:
        print(f"PubMed: all {len(cached)} organisms already cached.")
        return

    print(f"PubMed: fetching {len(todo)} (cache has {len(done)})...")
    rows = []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        rows.append(_fetch_pubmed_one(row["organism"]))
        if i % 25 == 0:
            # Checkpoint every 25
            ckpt = pd.concat([cached, pd.DataFrame(rows)], ignore_index=True) if not cached.empty else pd.DataFrame(rows)
            ckpt.to_csv(PUBMED_CSV, index=False)
            print(f"  {i}/{len(todo)} (checkpoint saved)")

    new_df   = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(PUBMED_CSV, index=False)
    print(f"PubMed: saved {len(combined)} rows.")


# ---------------------------------------------------------------------------
# SRA bioplastic expansion — find new organisms
# ---------------------------------------------------------------------------
# SRA bioplastic search terms — used with efetch runinfo to get organism names directly
BIOPLASTIC_SRA_SEARCHES = [
    ("PLA_degradation",     "polylactic acid[Title] AND degradation[Title]"),
    ("PHB_degradation",     "polyhydroxybutyrate[Title] AND degradation[Title]"),
    ("PHA_degradation",     "polyhydroxyalkanoate[Title] AND degradation[Title]"),
    ("PCL_degradation",     "polycaprolactone[Title] AND degradation[Title]"),
    ("PBSA_degradation",    "polybutylene succinate[Title] AND degradation[Title]"),
    ("bioplastic_deg",      "bioplastic[Title] AND degradation[Title]"),
    ("PETase",              "PETase[Title]"),
    ("MHETase",             "MHETase[Title]"),
    ("cutinase_plastic",    "cutinase[Title] AND (plastic[Title] OR polyester[Title])"),
    ("lipase_bioplastic",   "lipase[Title] AND (polylactic[Title] OR polycaprolactone[Title])"),
    ("plastic_enzyme",      "plastic degrading enzyme[Title]"),
    ("depolymerase_pha",    "depolymerase[Title] AND polyhydroxy[Title]"),
]


def _efetch_runinfo_organisms(ids: list[str]) -> list[str]:
    """Use efetch runinfo to get scientific names from a list of SRA run/experiment IDs."""
    if not ids:
        return []
    uid_str = ",".join(ids[:500])
    url = (f"{NCBI_BASE}/efetch.fcgi?db=sra&id={uid_str}"
           f"&rettype=runinfo&retmode=text")
    r = ncbi_get(url, sleep=0.5)
    if r is None:
        return []
    orgs = []
    for line in r.text.splitlines()[1:]:  # skip header
        parts = line.split(",")
        # RunInfo CSV columns: Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,
        #   size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,
        #   LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,
        #   Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,
        #   TaxID,ScientificName,...
        if len(parts) > 28:
            name = parts[28].strip().strip('"')
            if name and len(name) > 4 and name not in ("ScientificName", ""):
                orgs.append(name)
    return orgs


def run_sra_expand(existing_organisms: set):
    if EXTRA_CSV.exists():
        existing_extra = pd.read_csv(EXTRA_CSV, dtype=str)
        done_keys = set(existing_extra.get("search_key", pd.Series()).tolist())
    else:
        existing_extra = pd.DataFrame()
        done_keys = set()

    new_orgs: dict[str, dict] = {}

    for key, term in BIOPLASTIC_SRA_SEARCHES:
        if key in done_keys:
            print(f"  (already searched): {key}")
            continue

        print(f"  Searching SRA: {term}...")
        url = (f"{NCBI_BASE}/esearch.fcgi?db=sra"
               f"&term={requests.utils.quote(term)}&retmax=500&retmode=json&usehistory=y")
        r = ncbi_get(url, sleep=0.4)
        if r is None:
            done_keys.add(key)
            continue
        try:
            data  = r.json()["esearchresult"]
            count = int(data.get("count", 0))
            ids   = data.get("idlist", [])
        except Exception:
            done_keys.add(key)
            continue

        print(f"    {count} records found, fetching organism names for {len(ids)}...")

        if ids:
            found_orgs = _efetch_runinfo_organisms(ids)
            for org in found_orgs:
                if (org not in existing_organisms
                        and org not in new_orgs
                        and len(org) > 5
                        and not org[0].isdigit()
                        and org.replace(" ","").replace(".","").isalpha()):
                    new_orgs[org] = {
                        "organism":   org,
                        "search_key": key,
                        "search_term": term,
                        "sra_result_count": count,
                        "source": "SRA_bioplastic_efetch",
                    }

        done_keys.add(key)
        time.sleep(0.5)

    if not new_orgs:
        print("No new organisms found via SRA bioplastic search.")
        # Still mark searches as done so they don't repeat
        if not EXTRA_CSV.exists():
            placeholder = pd.DataFrame([{"organism":"","search_key":k,"search_term":"","sra_result_count":0,"source":"done"}
                                        for k in done_keys])
            placeholder.to_csv(EXTRA_CSV, index=False)
        return

    new_df   = pd.DataFrame(list(new_orgs.values()))
    combined = pd.concat([existing_extra, new_df], ignore_index=True) if not existing_extra.empty else new_df
    combined   = combined[combined["organism"].str.len() > 0]
    combined.to_csv(EXTRA_CSV, index=False)
    print(f"SRA expand: {len(new_orgs)} new organisms. Sample:")
    print(new_df["organism"].head(20).to_string())


# ---------------------------------------------------------------------------
# Enzyme family classification
# ---------------------------------------------------------------------------
ENZYME_FAMILIES = {
    "cutinase":    "Cutinase",
    "lipase":      "Lipase",
    "esterase":    "Esterase",
    "laccase":     "Laccase",
    "peroxidase":  "Peroxidase",
    "protease":    "Protease",
    "depolymerase":"Depolymerase",
    "hydrolase":   "Hydrolase",
    "petase":      "PETase",
    "mhetase":     "MHETase",
    "oxidase":     "Oxidase",
    "alkb":        "AlkB (alkane hydroxylase)",
    "cytochrome":  "Cytochrome P450",
    "monooxygenase":"Monooxygenase",
}


def classify_enzyme(name: str) -> str:
    if not name or pd.isna(name):
        return ""
    n = name.lower()
    for key, label in ENZYME_FAMILIES.items():
        if key in n:
            return label
    return "Other enzyme"


# ---------------------------------------------------------------------------
# HTML page generation with Chart.js
# ---------------------------------------------------------------------------
CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       background: #f2f5fb; color: #1a1a2e; line-height: 1.6; }
.header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
          color: #fff; padding: 28px 40px; }
.header h1 { font-size: 1.85rem; font-weight: 700; letter-spacing: -0.4px; }
.header .meta { color: #9fb3cc; font-size: 0.88rem; margin-top: 8px; }
.tag-row { margin-top: 10px; display: flex; flex-wrap: wrap; gap: 6px; }
.tag { display: inline-block; border-radius: 4px; padding: 3px 10px; font-size: 0.78rem;
       font-weight: 600; }
.tag-bio  { background: #1a7a4a; color: #d4f5e6; }
.tag-conv { background: #3a5cc0; color: #d4e3ff; }
.tag-both { background: #7a4a1a; color: #f5e6d4; }
.breadcrumb { padding: 10px 40px; background: #e6eaf3; font-size: 0.84rem; }
.breadcrumb a { color: #3a7bd5; text-decoration: none; }
.container { max-width: 1140px; margin: 0 auto; padding: 28px 40px; }
.section { background: #fff; border-radius: 10px; border: 1px solid #d8dff0;
           margin-bottom: 22px; overflow: hidden; box-shadow: 0 1px 4px rgba(0,0,0,.04); }
.section-header { background: #f0f4fb; border-bottom: 1px solid #d8dff0;
                  padding: 13px 22px; font-weight: 600; font-size: 0.98rem;
                  color: #1a1a2e; display: flex; align-items: center; gap: 8px; }
.section-body { padding: 20px 22px; }
.grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
.grid3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; }
.kv-row { display: flex; gap: 12px; margin-bottom: 9px; font-size: 0.9rem; }
.kv-label { color: #6b7a99; min-width: 175px; flex-shrink: 0; }
.kv-value { font-weight: 500; }
.pill { display: inline-block; border-radius: 4px; padding: 2px 8px; font-size: 0.78rem;
        margin: 2px; }
.pill-bio  { background: #d4f5e6; color: #155724; }
.pill-conv { background: #e8f0fe; color: #3a5cc0; }
.pill-other { background: #f0f0f0; color: #555; }
.badge { display: inline-block; border-radius: 4px; padding: 2px 8px;
         font-size: 0.77rem; font-weight: 600; }
.badge-yes { background: #d4edda; color: #155724; }
.badge-no  { background: #f8d7da; color: #721c24; }
.badge-na  { background: #e2e3e5; color: #383d41; }
table { width: 100%; border-collapse: collapse; font-size: 0.86rem; }
th { background: #f0f4fb; padding: 9px 13px; text-align: left;
     font-weight: 600; border-bottom: 2px solid #d8dff0; }
td { padding: 8px 13px; border-bottom: 1px solid #edf0f5; vertical-align: top; }
tr:last-child td { border-bottom: none; }
.stat-box { background: #f0f4fb; border-radius: 8px; padding: 16px; text-align: center; }
.stat-val { font-size: 1.75rem; font-weight: 700; color: #3a5cc0; }
.stat-lbl { font-size: 0.78rem; color: #6b7a99; margin-top: 4px; }
.no-data  { color: #9aa3b5; font-style: italic; font-size: 0.88rem; }
.note     { font-size: 0.78rem; color: #9aa3b5; margin-top: 10px; }
.chart-wrap { position: relative; height: 220px; margin-bottom: 6px; }
.chart-wrap-sm { position: relative; height: 160px; }
a { color: #3a7bd5; }
#search-input { width: 100%; padding: 10px 14px; font-size: 0.95rem;
  border: 1px solid #c8d0e4; border-radius: 6px; margin-bottom: 14px; outline: none; }
.filter-bar { display: flex; gap: 8px; margin-bottom: 12px; flex-wrap: wrap; }
.filter-btn { padding: 5px 14px; border: 1px solid #c8d0e4; border-radius: 20px;
              cursor: pointer; font-size: 0.83rem; background: #fff; }
.filter-btn.active { background: #3a5cc0; color: #fff; border-color: #3a5cc0; }
@media (max-width: 680px) {
  .grid2, .grid3 { grid-template-columns: 1fr; }
  .container { padding: 14px; }
  .header { padding: 18px 16px; }
}
"""

CHART_JS = "https://cdn.jsdelivr.net/npm/chart.js@4/dist/chart.umd.min.js"


def _slug(name: str) -> str:
    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def _badge(val) -> str:
    s = str(val).lower()
    if s in ("true","yes","1","true"): return '<span class="badge badge-yes">Yes</span>'
    if s in ("false","no","0","false"): return '<span class="badge badge-no">No</span>'
    return '<span class="badge badge-na">N/A</span>'


def _kv(label: str, value, pill: bool = False) -> str:
    if pill and value:
        parts = [v.strip() for v in str(value).split(";") if v.strip()]
        pills = "".join(f'<span class="pill pill-other">{p}</span>' for p in parts)
        return (f'<div class="kv-row"><span class="kv-label">{label}</span>'
                f'<span class="kv-value">{pills}</span></div>')
    return (f'<div class="kv-row"><span class="kv-label">{label}</span>'
            f'<span class="kv-value">{value if value else "N/A"}</span></div>')


def _fmt_bases(b) -> str:
    try:
        b = float(b)
        if b >= 1e9:  return f"{b/1e9:.2f} Gbp"
        if b >= 1e6:  return f"{b/1e6:.2f} Mbp"
        if b >= 1e3:  return f"{b/1e3:.0f} kbp"
        return f"{b:.0f} bp"
    except Exception: return "N/A"


def _fmt_genome(b) -> str:
    try:
        b = float(b)
        if b >= 1e9:  return f"{b/1e9:.2f} Gbp"
        if b >= 1e6:  return f"{b/1e6:.2f} Mbp"
        return f"{b/1e3:.0f} kbp"
    except Exception: return "N/A"


def _plastic_pills(plastics_list, df_org=None) -> str:
    if not plastics_list:
        return "<span class='no-data'>none recorded</span>"
    parts = []
    for p in plastics_list:
        cls = classify_plastic(p)
        css = {"bioplastic": "pill-bio", "conventional": "pill-conv"}.get(cls, "pill-other")
        parts.append(f'<span class="pill {css}" title="{cls}">{p}</span>')
    return "".join(parts)


def _make_charts(ent: pd.DataFrame, org_name: str) -> tuple[str, str]:
    """Return (chart_html, chart_script) for the organism."""
    # 1. Plastic breakdown (horizontal bar)
    plastic_counts = ent.groupby("plastic")["organism"].count().sort_values(ascending=True)
    p_labels  = [p for p in plastic_counts.index.tolist()]
    p_values  = plastic_counts.values.tolist()
    p_colors  = [
        "#27ae60" if classify_plastic(p) == "bioplastic"
        else "#3a5cc0" if classify_plastic(p) == "conventional"
        else "#95a5a6"
        for p in p_labels
    ]

    # 2. Research timeline (bar by year)
    year_counts = ent.dropna(subset=["year"])
    if not year_counts.empty:
        year_counts = year_counts.groupby("year")["organism"].count()
        year_counts.index = year_counts.index.astype(int)
        year_counts = year_counts.sort_index()
        y_labels = year_counts.index.tolist()
        y_values = year_counts.values.tolist()
    else:
        y_labels, y_values = [], []

    # 3. Evidence methods (top methods, simplified)
    ev_series = ent["evidence"].dropna()
    methods   = {}
    for ev in ev_series:
        for m in re.split(r"[;,]", ev):
            m = m.strip()
            if m and len(m) > 1:
                methods[m] = methods.get(m, 0) + 1
    top_methods = sorted(methods.items(), key=lambda x: -x[1])[:8]
    ev_labels = [x[0] for x in top_methods]
    ev_values = [x[1] for x in top_methods]
    ev_colors = [
        "#3a5cc0","#27ae60","#e67e22","#e74c3c","#9b59b6",
        "#1abc9c","#f39c12","#95a5a6",
    ][:len(ev_labels)]

    # Enzyme families
    enz_series = ent["enzyme_name"].dropna()
    fam_counts: dict = {}
    for en in enz_series:
        fam = classify_enzyme(en)
        if fam:
            fam_counts[fam] = fam_counts.get(fam, 0) + 1
    fam_labels = list(fam_counts.keys())
    fam_values = list(fam_counts.values())

    chart_html = """
<div class="grid2" style="margin-bottom:20px">
  <div>
    <div style="font-size:0.82rem;font-weight:600;color:#6b7a99;margin-bottom:8px">
      Plastic types degraded
      <span style="color:#27ae60">&#9632;</span> bioplastic
      <span style="color:#3a5cc0">&#9632;</span> conventional
    </div>
    <div class="chart-wrap"><canvas id="plastic-chart"></canvas></div>
  </div>
  <div>
    <div style="font-size:0.82rem;font-weight:600;color:#6b7a99;margin-bottom:8px">
      Research publications by year
    </div>
    <div class="chart-wrap"><canvas id="year-chart"></canvas></div>
  </div>
</div>
<div class="grid2">
  <div>
    <div style="font-size:0.82rem;font-weight:600;color:#6b7a99;margin-bottom:8px">
      Evidence methods used
    </div>
    <div class="chart-wrap-sm"><canvas id="ev-chart"></canvas></div>
  </div>
  <div>
    <div style="font-size:0.82rem;font-weight:600;color:#6b7a99;margin-bottom:8px">
      Enzyme families identified
    </div>
    <div class="chart-wrap-sm"><canvas id="enz-chart"></canvas></div>
  </div>
</div>"""

    p_labels_js  = json.dumps(p_labels)
    p_values_js  = json.dumps(p_values)
    p_colors_js  = json.dumps(p_colors)
    y_labels_js  = json.dumps(y_labels)
    y_values_js  = json.dumps(y_values)
    ev_labels_js = json.dumps(ev_labels)
    ev_values_js = json.dumps(ev_values)
    ev_colors_js = json.dumps(ev_colors)
    fam_labels_js = json.dumps(fam_labels)
    fam_values_js = json.dumps(fam_values)

    chart_script = f"""
<script src="{CHART_JS}"></script>
<script>
Chart.defaults.font.family = "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif";
Chart.defaults.font.size = 11;

// Plastic types chart
new Chart(document.getElementById('plastic-chart'), {{
  type: 'bar',
  data: {{
    labels: {p_labels_js},
    datasets: [{{
      data: {p_values_js},
      backgroundColor: {p_colors_js},
      borderWidth: 0
    }}]
  }},
  options: {{
    indexAxis: 'y',
    responsive: true,
    maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ grid: {{ color: '#edf0f5' }}, title: {{ display: true, text: 'Entries', font: {{size:10}} }} }},
      y: {{ grid: {{ display: false }}, ticks: {{ font: {{ size: 10 }} }} }}
    }}
  }}
}});

// Research timeline chart
new Chart(document.getElementById('year-chart'), {{
  type: 'bar',
  data: {{
    labels: {y_labels_js},
    datasets: [{{
      data: {y_values_js},
      backgroundColor: '#3a5cc0',
      borderWidth: 0,
      borderRadius: 2
    }}]
  }},
  options: {{
    responsive: true,
    maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ grid: {{ display: false }}, ticks: {{ maxRotation: 45, font: {{ size: 9 }} }} }},
      y: {{ grid: {{ color: '#edf0f5' }}, ticks: {{ stepSize: 1, font: {{ size: 10 }} }},
            title: {{ display: true, text: 'Entries', font: {{ size: 10 }} }} }}
    }}
  }}
}});

// Evidence methods chart
new Chart(document.getElementById('ev-chart'), {{
  type: 'bar',
  data: {{
    labels: {ev_labels_js},
    datasets: [{{
      data: {ev_values_js},
      backgroundColor: {ev_colors_js},
      borderWidth: 0,
      borderRadius: 2
    }}]
  }},
  options: {{
    indexAxis: 'y',
    responsive: true,
    maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ grid: {{ color: '#edf0f5' }}, ticks: {{ stepSize: 1, font: {{ size: 9 }} }} }},
      y: {{ grid: {{ display: false }}, ticks: {{ font: {{ size: 9 }} }} }}
    }}
  }}
}});

// Enzyme families chart
new Chart(document.getElementById('enz-chart'), {{
  type: 'doughnut',
  data: {{
    labels: {fam_labels_js},
    datasets: [{{
      data: {fam_values_js},
      backgroundColor: ['#3a5cc0','#27ae60','#e67e22','#e74c3c','#9b59b6','#1abc9c'],
      borderWidth: 2,
      borderColor: '#fff'
    }}]
  }},
  options: {{
    responsive: true,
    maintainAspectRatio: false,
    plugins: {{
      legend: {{ position: 'right', labels: {{ boxWidth: 10, font: {{ size: 10 }} }} }}
    }}
  }}
}});
</script>"""

    return chart_html, chart_script


def _make_page(
    org_row: pd.Series,
    full_df: pd.DataFrame,
    sra_row,
    bd_row,
    genome_row,
    pubmed_row,
    is_extra: bool = False,
) -> str:
    import ast

    name = org_row["organism"]
    ent  = full_df[full_df["organism"] == name].copy()

    # Bioplastic relevance
    n_bio  = int(org_row.get("n_bioplastic", 0) or 0)
    n_conv = int(org_row.get("n_conventional", 0) or 0)
    if n_bio > 0 and n_conv > 0:
        bio_tag = '<span class="tag tag-both">Bioplastic + Conventional</span>'
    elif n_bio > 0:
        bio_tag = '<span class="tag tag-bio">Bioplastic research</span>'
    elif n_conv > 0:
        bio_tag = '<span class="tag tag-conv">Conventional plastic</span>'
    else:
        bio_tag = ""

    extra_tag = '<span class="tag" style="background:#5a1a7a;color:#f0d4f5">SRA-expanded</span>' if is_extra else ""

    # Plastics list
    plastics = org_row.get("plastics", [])
    if isinstance(plastics, str):
        try:    plastics = ast.literal_eval(plastics)
        except: plastics = [p.strip() for p in plastics.split(",")]

    # Entries table
    rows_html = ""
    for _, e in ent.iterrows():
        doi = str(e.get("doi","") or "")
        doi_link = (f'<a href="https://doi.org/{doi}" target="_blank">{doi[:35]}</a>'
                    if len(doi) > 5 else doi)
        yr  = int(e["year"]) if pd.notna(e.get("year")) else "N/A"
        enz = str(e.get("enzyme_name","") or "")
        enz_fam = classify_enzyme(enz)
        cls = classify_plastic(str(e.get("plastic","") or ""))
        plastic_pill_css = {"bioplastic":"pill-bio","conventional":"pill-conv"}.get(cls,"pill-other")
        rows_html += (
            f"<tr>"
            f"<td><span class='pill {plastic_pill_css}'>{e.get('plastic','')}</span></td>"
            f"<td>{yr}</td>"
            f"<td>{enz[:45] or 'N/A'}</td>"
            f"<td style='font-size:0.78rem'>{enz_fam}</td>"
            f"<td>{_badge(e.get('has_sequence',False))}</td>"
            f"<td>{_badge(e.get('has_genbank',False))}</td>"
            f"<td style='font-size:0.78rem'>{doi_link}</td>"
            f"</tr>"
        )
    if rows_html:
        entries_html = (
            "<table><thead><tr>"
            "<th>Plastic</th><th>Year</th><th>Enzyme</th><th>Enzyme family</th>"
            "<th>Sequence</th><th>GenBank</th><th>DOI</th>"
            "</tr></thead><tbody>" + rows_html + "</tbody></table>"
        )
    else:
        entries_html = '<p class="no-data">No PlasticDB entries (organism added via SRA expansion).</p>'

    # SRA section
    if sra_row is not None:
        try: rc = int(float(sra_row.get("sra_run_count","")))
        except: rc = "N/A"
        sra_html = (
            f'<div class="grid3" style="margin-bottom:18px">'
            f'<div class="stat-box"><div class="stat-val">{rc}</div>'
            f'<div class="stat-lbl">SRA runs deposited</div></div>'
            f'<div class="stat-box"><div class="stat-val">{_fmt_bases(sra_row.get("sra_total_bases"))}</div>'
            f'<div class="stat-lbl">Bases (top 5 runs sampled)</div></div>'
            f'<div class="stat-box"><div class="stat-val">{sra_row.get("sra_date_range","") or "N/A"}</div>'
            f'<div class="stat-lbl">Deposit year range</div></div>'
            f'</div>'
            + _kv("Sequencing platforms", sra_row.get("sra_platforms","") or "None in sample", pill=True)
            + _kv("Library strategies",  sra_row.get("sra_strategies","") or "None in sample", pill=True)
            + '<p class="note">Source: NCBI SRA E-utilities esearch + esummary. Run count is exact from esearch. Platform/strategy sampled from top 5 runs.</p>'
        )
    else:
        sra_html = '<p class="no-data">SRA data not yet fetched for this organism.</p>'

    # BacDive section
    if bd_row is not None and str(bd_row.get("bacdive_found","")).lower() == "yes":
        bd_url  = bd_row.get("bacdive_url","")
        bd_link = (f'<a href="{bd_url}" target="_blank">BacDive strain {bd_row.get("bacdive_strain_id","")}</a>'
                   if bd_url else "N/A")
        bd_html = (
            _kv("BacDive record",      bd_link)
            + _kv("Culture temp (C)", bd_row.get("bacdive_temp_c","") or "Not recorded")
            + _kv("pH",               bd_row.get("bacdive_ph","")     or "Not recorded")
            + _kv("Oxygen tolerance", bd_row.get("bacdive_oxygen","") or "Not recorded")
            + _kv("Morphology",       bd_row.get("bacdive_morphology","") or "Not recorded")
            + _kv("Isolation source", bd_row.get("bacdive_isolation","")  or "Not recorded")
            + '<p class="note">Source: BacDive public strain page (bacdive.dsmz.de).</p>'
        )
    elif bd_row is not None:
        bd_html = '<p class="no-data">Organism not found in BacDive public database.</p>'
    else:
        bd_html = '<p class="no-data">BacDive data not yet fetched.</p>'

    # Genome section
    if genome_row is not None and genome_row.get("assembly_accession",""):
        acc        = genome_row.get("assembly_accession","")
        acc_link   = f'<a href="https://www.ncbi.nlm.nih.gov/assembly/{acc}" target="_blank">{acc}</a>'
        gsize      = _fmt_genome(genome_row.get("genome_size_bp",""))
        alevel     = genome_row.get("assembly_level","") or "N/A"
        n50_val    = genome_row.get("contig_n50","") or genome_row.get("scaffold_n50","")
        n50_str    = _fmt_genome(n50_val) if n50_val else "N/A"
        ctg_count  = genome_row.get("contig_count","") or "N/A"
        cov_raw    = genome_row.get("coverage","")
        coverage   = str(cov_raw) if (cov_raw and str(cov_raw) not in ("","nan","N/A")) else "N/A"
        genome_html = (
            f'<div class="grid3" style="margin-bottom:18px">'
            f'<div class="stat-box"><div class="stat-val">{gsize}</div>'
            f'<div class="stat-lbl">Genome size</div></div>'
            f'<div class="stat-box"><div class="stat-val" style="font-size:1.1rem">{alevel}</div>'
            f'<div class="stat-lbl">Assembly level</div></div>'
            f'<div class="stat-box"><div class="stat-val">{n50_str}</div>'
            f'<div class="stat-lbl">Contig N50</div></div>'
            f'</div>'
            + _kv("Assembly accession", acc_link)
            + _kv("Assembly name",      genome_row.get("assembly_name","") or "N/A")
            + _kv("Contig count",       ctg_count)
            + _kv("Sequencing coverage", coverage + "x" if coverage != "N/A" else "N/A")
            + '<p class="note">Source: NCBI Assembly. Best assembly selected by completeness (Complete > Chromosome > Scaffold > Contig).</p>'
        )
    elif genome_row is not None:
        genome_html = '<p class="no-data">No genome assembly found in NCBI Assembly for this organism.</p>'
    else:
        genome_html = '<p class="no-data">Genome data not yet fetched.</p>'

    # PubMed section
    if pubmed_row is not None:
        plastic_papers = pubmed_row.get("pubmed_plastic_count","")
        total_papers   = pubmed_row.get("pubmed_total_count","")
        try: plastic_papers = int(float(plastic_papers))
        except: plastic_papers = "N/A"
        try: total_papers = int(float(total_papers))
        except: total_papers = "N/A"
        pubmed_html = (
            f'<div class="grid2" style="margin-bottom:18px">'
            f'<div class="stat-box"><div class="stat-val">{plastic_papers}</div>'
            f'<div class="stat-lbl">Plastic/bioplastic papers</div></div>'
            f'<div class="stat-box"><div class="stat-val">{total_papers}</div>'
            f'<div class="stat-lbl">Total PubMed papers</div></div>'
            f'</div>'
            + '<p class="note">Source: NCBI PubMed esearch. Plastic/bioplastic count uses query: organism + (plastic* OR bioplastic* OR polymer* OR degradation).</p>'
        )
    else:
        pubmed_html = '<p class="no-data">PubMed data not yet fetched.</p>'

    # Charts
    chart_html, chart_script = _make_charts(ent, name) if not ent.empty else ("", "")

    fy = int(org_row["first_year"]) if pd.notna(org_row.get("first_year")) else "N/A"
    ly = int(org_row["last_year"])  if pd.notna(org_row.get("last_year"))  else "N/A"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{name} | Bioplastic Organism Profile</title>
<style>{CSS}</style>
</head>
<body>
<div class="header">
  <h1><em>{name}</em></h1>
  <div class="meta">
    Tax ID: {org_row.get("tax_id","N/A")} | {org_row.get("n_plastics",0)} plastic type(s) |
    {org_row.get("n_entries",0)} PlasticDB entries | Years: {fy} to {ly}
  </div>
  <div class="tag-row">{bio_tag}{extra_tag}</div>
</div>
<div class="breadcrumb"><a href="index.html">All organisms</a> / {name}</div>
<div class="container">

<!-- Charts -->
{"<div class='section'><div class='section-header'>Research Overview</div><div class='section-body'>" + chart_html + "</div></div>" if chart_html else ""}

<!-- PlasticDB Summary -->
<div class="section"><div class="section-header">PlasticDB Summary</div>
<div class="section-body"><div class="grid2"><div>
{_kv("Plastics degraded", _plastic_pills(plastics))}
{_kv("Bioplastic entries", str(n_bio) + " / " + str(org_row.get("n_entries",0)) if n_bio else "0 / " + str(org_row.get("n_entries",0)))}
{_kv("Thermophilic flag", org_row.get("thermophilic","") or "Not recorded")}
{_kv("Has linked sequence", _badge(org_row.get("has_sequence")))}
{_kv("Has named enzyme",    _badge(org_row.get("has_enzyme")))}
{_kv("Has GenBank ID",      _badge(org_row.get("has_genbank")))}
</div><div>
{_kv("Isolation environments", org_row.get("isolation_envs","") or "Not recorded", pill=True)}
{_kv("Isolation locations",    org_row.get("isolation_locs","") or "Not recorded", pill=True)}
</div></div></div></div>

<!-- Entries table -->
<div class="section"><div class="section-header">PlasticDB Entries ({org_row.get("n_entries",0)} total)</div>
<div class="section-body">{entries_html}</div></div>

<!-- Genome data -->
<div class="section"><div class="section-header">Genome Data (NCBI Assembly)</div>
<div class="section-body">{genome_html}</div></div>

<!-- PubMed literature -->
<div class="section"><div class="section-header">PubMed Literature</div>
<div class="section-body">{pubmed_html}</div></div>

<!-- SRA sequencing -->
<div class="section"><div class="section-header">NCBI SRA Sequencing Data</div>
<div class="section-body">{sra_html}</div></div>

<!-- BacDive physiological -->
<div class="section"><div class="section-header">BacDive Physiological Data</div>
<div class="section-body">{bd_html}</div></div>

</div>
{chart_script}
</body>
</html>"""


def run_pages(full_df: pd.DataFrame, orgs: pd.DataFrame,
              sra_df: pd.DataFrame, bd_df: pd.DataFrame,
              genome_df: pd.DataFrame, pubmed_df: pd.DataFrame,
              extra_df: pd.DataFrame):
    import ast

    sra_idx    = sra_df.set_index("organism")    if not sra_df.empty    else pd.DataFrame()
    bd_idx     = bd_df.set_index("organism")     if not bd_df.empty     else pd.DataFrame()
    genome_idx = genome_df.set_index("organism") if not genome_df.empty else pd.DataFrame()
    pubmed_idx = pubmed_df.set_index("organism") if not pubmed_df.empty else pd.DataFrame()

    # Build combined organism list (PlasticDB + extra SRA organisms)
    all_orgs = orgs.copy()
    extra_orgs = set()
    if not extra_df.empty:
        for _, row in extra_df.iterrows():
            org = row["organism"]
            if org not in set(orgs["organism"]):
                extra_orgs.add(org)
                new_row = pd.Series({
                    "organism": org, "tax_id": "", "n_entries": 0,
                    "plastics": [], "n_plastics": 0, "has_sequence": False,
                    "has_enzyme": False, "has_genbank": False,
                    "first_year": None, "last_year": None,
                    "thermophilic": "", "isolation_envs": "", "isolation_locs": "",
                    "n_bioplastic": 0, "n_conventional": 0, "bioplastic_relevant": False,
                })
                all_orgs = pd.concat([all_orgs, new_row.to_frame().T], ignore_index=True)

    print(f"Generating pages for {len(all_orgs)} organisms ({len(extra_orgs)} SRA-expanded)...")

    for i, (_, row) in enumerate(all_orgs.iterrows(), 1):
        name    = row["organism"]
        s       = _slug(name)
        is_ext  = name in extra_orgs
        sra_row    = sra_idx.loc[name].to_dict()    if (not sra_idx.empty    and name in sra_idx.index)    else None
        bd_row     = bd_idx.loc[name].to_dict()     if (not bd_idx.empty     and name in bd_idx.index)     else None
        genome_row = genome_idx.loc[name].to_dict() if (not genome_idx.empty and name in genome_idx.index) else None
        pubmed_row = pubmed_idx.loc[name].to_dict() if (not pubmed_idx.empty and name in pubmed_idx.index) else None
        html = _make_page(row, full_df, sra_row, bd_row, genome_row, pubmed_row, is_ext)
        (PAGES_DIR / f"{s}.html").write_text(html, encoding="utf-8")
        if i % 100 == 0:
            print(f"  {i}/{len(all_orgs)} pages written")

    # Index page
    sra_map    = dict(zip(sra_df["organism"],    sra_df.get("sra_run_count", pd.Series()))) if not sra_df.empty else {}
    bd_map     = dict(zip(bd_df["organism"],     bd_df.get("bacdive_found",  pd.Series()))) if not bd_df.empty  else {}
    genome_map = dict(zip(genome_df["organism"], genome_df.get("assembly_level", pd.Series()))) if not genome_df.empty else {}
    pubmed_map = dict(zip(pubmed_df["organism"], pubmed_df.get("pubmed_plastic_count", pd.Series()))) if not pubmed_df.empty else {}

    rows_html = ""
    for _, row in all_orgs.sort_values("n_entries", ascending=False).iterrows():
        org = row["organism"]
        s   = _slug(org)
        try: rc = int(float(sra_map.get(org,"")))
        except: rc = "N/A"
        bd_found = str(bd_map.get(org,"")).lower() == "yes"
        has_genome = bool(genome_map.get(org,""))
        try: pm = int(float(pubmed_map.get(org,"")))
        except: pm = "N/A"

        plastics = row.get("plastics",[])
        if isinstance(plastics, str):
            try:    plastics = ast.literal_eval(plastics)
            except: plastics = [p.strip() for p in plastics.split(",")]

        n_bio  = int(row.get("n_bioplastic",0) or 0)
        n_conv = int(row.get("n_conventional",0) or 0)
        if n_bio > 0 and n_conv > 0:
            bio_badge = '<span class="badge" style="background:#f5e6d4;color:#7a4a1a">Both</span>'
        elif n_bio > 0:
            bio_badge = '<span class="badge" style="background:#d4f5e6;color:#1a7a4a">Bioplastic</span>'
        elif n_conv > 0:
            bio_badge = '<span class="badge" style="background:#e8f0fe;color:#3a5cc0">Conv.</span>'
        else:
            bio_badge = ""

        pills = ""
        for p in plastics[:4]:
            cls = classify_plastic(p)
            css = {"bioplastic":"pill-bio","conventional":"pill-conv"}.get(cls,"pill-other")
            pills += f'<span class="pill {css}">{p}</span>'
        if len(plastics) > 4:
            pills += f'<span class="pill pill-other">+{len(plastics)-4}</span>'

        fy = int(row["first_year"]) if pd.notna(row.get("first_year")) else "N/A"
        is_ext = org in extra_orgs
        ext_icon = " <span style='color:#9b59b6;font-size:0.7rem'>SRA</span>" if is_ext else ""

        rows_html += (
            f'<tr data-bio="{n_bio > 0}" data-conv="{n_conv > 0}" data-org="{org.lower()}">'
            f'<td><a href="{s}.html">{org}</a>{ext_icon}</td>'
            f'<td style="text-align:center">{bio_badge}</td>'
            f'<td>{pills}</td>'
            f'<td style="text-align:center">{row.get("n_entries",0) or ""}</td>'
            f'<td style="text-align:center">{rc}</td>'
            f'<td style="text-align:center">{pm}</td>'
            f'<td style="text-align:center">{"Yes" if has_genome else ""}</td>'
            f'<td style="text-align:center">{"Yes" if bd_found else ""}</td>'
            f'<td style="text-align:center">{fy}</td>'
            f'</tr>'
        )

    total = len(all_orgs)
    bio_count  = int((all_orgs["n_bioplastic"].fillna(0).astype(int) > 0).sum())
    sra_count  = sum(1 for v in sra_map.values() if float(v or 0) > 0)
    bd_count   = sum(1 for v in bd_map.values() if str(v).lower() == "yes")
    gen_count  = sum(1 for v in genome_map.values() if v)

    index_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Bioplastic Organism Profiles</title>
<style>{CSS}
.stat-strip {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px; margin-bottom: 22px; }}
</style>
</head>
<body>
<div class="header">
  <h1>Bioplastic Organism Profiles</h1>
  <div class="meta">
    {total} organisms sourced from PlasticDB and NCBI SRA bioplastic-specific searches.
    Cross-referenced with NCBI SRA, NCBI Assembly, PubMed, and BacDive.
    All values from live public APIs. Nothing fabricated or estimated.
  </div>
  <div class="tag-row">
    <span class="tag tag-bio">Bioplastic degraders</span>
    <span class="tag tag-conv">Conventional plastic degraders</span>
    <span class="tag tag-both">Both</span>
  </div>
</div>
<div class="container">

<div class="stat-strip">
  <div class="stat-box"><div class="stat-val">{total}</div><div class="stat-lbl">Total organisms</div></div>
  <div class="stat-box"><div class="stat-val">{bio_count}</div><div class="stat-lbl">With bioplastic entries</div></div>
  <div class="stat-box"><div class="stat-val">{sra_count}</div><div class="stat-lbl">With SRA runs</div></div>
  <div class="stat-box"><div class="stat-val">{gen_count}</div><div class="stat-lbl">With genome assembly</div></div>
</div>

<div class="section">
  <div class="section-header">All Organisms</div>
  <div class="section-body">
    <input id="search-input" placeholder="Filter by organism name..." oninput="filterTable()">
    <div class="filter-bar">
      <button class="filter-btn active" onclick="setFilter('all',this)">All ({total})</button>
      <button class="filter-btn" onclick="setFilter('bio',this)">Bioplastic ({bio_count})</button>
      <button class="filter-btn" onclick="setFilter('conv',this)">Conventional only</button>
      <button class="filter-btn" onclick="setFilter('sra',this)">SRA-expanded</button>
    </div>
    <table id="org-table">
      <thead><tr>
        <th>Organism</th>
        <th style="text-align:center">Relevance</th>
        <th>Plastics</th>
        <th style="text-align:center">PlasticDB</th>
        <th style="text-align:center">SRA runs</th>
        <th style="text-align:center">PubMed</th>
        <th style="text-align:center">Genome</th>
        <th style="text-align:center">BacDive</th>
        <th style="text-align:center">First year</th>
      </tr></thead>
      <tbody>{rows_html}</tbody>
    </table>
  </div>
</div>
</div>

<script>
let currentFilter = 'all';

function setFilter(f, btn) {{
  currentFilter = f;
  document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  filterTable();
}}

function filterTable() {{
  const q = document.getElementById('search-input').value.toLowerCase();
  document.querySelectorAll('#org-table tbody tr').forEach(r => {{
    const org  = r.dataset.org || '';
    const bio  = r.dataset.bio === 'true';
    const conv = r.dataset.conv === 'true';
    const isSra = r.cells[0].querySelector('span') !== null;

    let show = org.includes(q);
    if (show) {{
      if (currentFilter === 'bio')  show = bio;
      if (currentFilter === 'conv') show = conv && !bio;
      if (currentFilter === 'sra')  show = isSra;
    }}
    r.style.display = show ? '' : 'none';
  }});
}}
</script>
</body>
</html>"""

    (PAGES_DIR / "index.html").write_text(index_html, encoding="utf-8")
    total_pages = len(list(PAGES_DIR.glob("*.html")))
    print(f"index.html written. Total pages: {total_pages}.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode",  choices=["genome","fix-genome","pubmed","sra-expand","pages","status"], required=True)
    parser.add_argument("--batch", type=int, default=80)
    args = parser.parse_args()

    full_df = load_plasticdb()
    orgs    = build_org_summary(full_df)
    print(f"Loaded {len(full_df):,} PlasticDB entries, {len(orgs):,} organisms.")
    bio_relevant = int(orgs["bioplastic_relevant"].sum())
    print(f"Bioplastic-relevant organisms: {bio_relevant}/{len(orgs)}")

    def _load(csv, *extra_orgs_dfs):
        if csv.exists():
            return pd.read_csv(csv, dtype=str)
        return pd.DataFrame()

    if args.mode == "status":
        print(f"  SRA:      {len(_load(SRA_CSV))}/{len(orgs)}")
        print(f"  BacDive:  {len(_load(BACDIVE_CSV))}/{len(orgs)}")
        print(f"  Genome:   {len(_load(GENOME_CSV))}/{len(orgs)}")
        print(f"  PubMed:   {len(_load(PUBMED_CSV))}/{len(orgs)}")
        extra = _load(EXTRA_CSV)
        print(f"  SRA-extra:{len(extra)} new organisms")
        print(f"  Pages:    {len(list(PAGES_DIR.glob('*.html')))}")

    elif args.mode == "genome":
        run_genome_batch(orgs, args.batch)

    elif args.mode == "fix-genome":
        # Re-fetch stats for organisms whose accession was found but stats are missing.
        # Uses one esummary call per organism (no esearch needed).
        cached = pd.read_csv(GENOME_CSV, dtype=str) if GENOME_CSV.exists() else pd.DataFrame()
        if cached.empty:
            print("No genome cache to fix.")
        else:
            # Identify rows with accession but missing genome_size_bp
            needs_fix = cached[
                cached["assembly_accession"].notna() &
                (cached["assembly_accession"].str.len() > 3) &
                (cached.get("genome_size_bp", pd.Series(dtype=str)).isna() |
                 (cached.get("genome_size_bp", pd.Series(dtype=str)) == ""))
            ]
            if "genome_size_bp" not in cached.columns:
                needs_fix = cached[
                    cached["assembly_accession"].notna() &
                    (cached["assembly_accession"].str.len() > 3)
                ]
            print(f"fix-genome: {len(needs_fix)} rows need stats refresh...")
            todo = needs_fix.head(args.batch)
            updated = 0
            for i, (idx, row) in enumerate(todo.iterrows(), 1):
                acc = row["assembly_accession"]
                url = (f"{NCBI_BASE}/esummary.fcgi?db=assembly"
                       f"&id={requests.utils.quote(acc)}&retmode=json")
                r = ncbi_get(url)
                if r is None:
                    continue
                try:
                    result = r.json().get("result", {})
                except Exception:
                    continue
                for k, entry in result.items():
                    if k != "uids":
                        parsed = _parse_assembly_entry(entry)
                        for col, val in parsed.items():
                            cached.at[idx, col] = val
                        updated += 1
                        break
                if i % 20 == 0:
                    cached.to_csv(GENOME_CSV, index=False)
                    print(f"  {i}/{len(todo)} (checkpoint)")
            cached.to_csv(GENOME_CSV, index=False)
            print(f"fix-genome: updated {updated} rows.")

    elif args.mode == "pubmed":
        run_pubmed_batch(orgs, args.batch)

    elif args.mode == "sra-expand":
        run_sra_expand(set(orgs["organism"].tolist()))

    elif args.mode == "pages":
        sra_df    = _load(SRA_CSV)
        bd_df     = _load(BACDIVE_CSV)
        genome_df = _load(GENOME_CSV)
        pubmed_df = _load(PUBMED_CSV)
        extra_df  = _load(EXTRA_CSV)
        run_pages(full_df, orgs, sra_df, bd_df, genome_df, pubmed_df, extra_df)
