"""
Build script for organism-profiles.

Fetches real data from:
  - NCBI SRA (E-utilities) -- run counts, platforms, library strategies, base totals
  - BacDive (public HTML pages) -- culture temperature, pH, oxygen tolerance,
    Gram stain, cell morphology, isolation source

Generates:
  - data/sra_stats.csv         (cached; re-run skips already-fetched rows)
  - data/bacdive_data.csv      (cached; re-run skips already-fetched rows)
  - pages/<slug>.html          (one HTML page per organism)
  - pages/index.html           (searchable master index)

Run from the organism-profiles/ directory:
    python build.py

Flags:
    --sra-only        skip BacDive fetching
    --pages-only      skip all fetching, only regenerate HTML from cached CSVs
    --limit N         process only the first N organisms (useful for testing)
"""

import argparse
import re
import sys
import time
import threading
import unicodedata
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT      = Path(__file__).parent
DATA_DIR  = ROOT / "data"
PAGES_DIR = ROOT / "pages"
PDB_TSV   = ROOT.parent / "plastic-biodegradation-analysis" / "data" / "plasticdb_microorganisms.tsv"

DATA_DIR.mkdir(exist_ok=True)
PAGES_DIR.mkdir(exist_ok=True)

SRA_CSV     = DATA_DIR / "sra_stats.csv"
BACDIVE_CSV = DATA_DIR / "bacdive_data.csv"

# ---------------------------------------------------------------------------
# Rate-limited HTTP session (shared across threads)
# ---------------------------------------------------------------------------
class RateLimitedSession:
    """Wraps requests.Session with a global rate limit (requests per second)."""

    def __init__(self, rps: float = 2.5, retries: int = 3):
        self._min_interval = 1.0 / rps
        self._last = 0.0
        self._lock = threading.Lock()
        self._session = requests.Session()
        self._session.headers["User-Agent"] = (
            "organism-profiles-builder/1.0 (igem-toronto.org research)"
        )
        self._retries = retries

    def get(self, url: str, **kwargs) -> requests.Response:
        with self._lock:
            wait = self._min_interval - (time.time() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.time()

        kwargs.setdefault("timeout", 20)
        for attempt in range(self._retries):
            try:
                resp = self._session.get(url, **kwargs)
                return resp
            except requests.RequestException:
                if attempt == self._retries - 1:
                    raise
                time.sleep(2 ** attempt)


NCBI_SESSION    = RateLimitedSession(rps=2.5)
BACDIVE_SESSION = RateLimitedSession(rps=1.5)

# ---------------------------------------------------------------------------
# Load PlasticDB
# ---------------------------------------------------------------------------
def load_plasticdb() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(PDB_TSV, sep="\t", dtype=str, on_bad_lines="skip")
    df.columns = [
        "organism","tax_id","plastic","reference","enzyme_name",
        "enzyme_id","db_enzyme_name","gene","genbank_id","sequence",
        "year","evidence","plastic_used","manufacturer","analytical_grade",
        "thermophilic","isolation_sample","isolation_environment",
        "isolation_location","extrapolated_from_enzyme","enzyme_id_in_paper","doi",
    ]
    df["year"]         = pd.to_numeric(df["year"], errors="coerce")
    df["has_sequence"] = df["sequence"].notna() & (df["sequence"].str.len() > 10)
    df["has_enzyme"]   = df["enzyme_name"].notna() & (df["enzyme_name"] != "")
    df["has_genbank"]  = df["genbank_id"].notna() & (df["genbank_id"].str.len() > 3)
    df["organism"]     = df["organism"].str.strip()
    df["plastic"]      = df["plastic"].str.strip()

    orgs = (
        df.groupby("organism")
        .agg(
            tax_id          = ("tax_id", "first"),
            n_entries       = ("organism", "count"),
            plastics        = ("plastic", lambda x: sorted(x.dropna().unique().tolist())),
            n_plastics      = ("plastic", "nunique"),
            has_sequence    = ("has_sequence", "any"),
            has_enzyme      = ("has_enzyme", "any"),
            has_genbank     = ("has_genbank", "any"),
            first_year      = ("year", "min"),
            last_year       = ("year", "max"),
            thermophilic    = ("thermophilic", lambda x: x.mode()[0] if len(x.mode()) else ""),
            isolation_envs  = ("isolation_environment",
                               lambda x: "; ".join(sorted(x.dropna().unique()))),
            isolation_locs  = ("isolation_location",
                               lambda x: "; ".join(sorted(x.dropna().unique()))),
        )
        .reset_index()
    )
    return df, orgs


# ---------------------------------------------------------------------------
# SRA fetch
# ---------------------------------------------------------------------------
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def _sra_search(organism: str) -> dict:
    """Return SRA run count and platform/strategy breakdown for one organism."""
    term = f'"{organism}"[Organism]'
    url  = (f"{NCBI_BASE}/esearch.fcgi"
            f"?db=sra&term={requests.utils.quote(term)}&retmax=200&retmode=json")
    try:
        r = NCBI_SESSION.get(url)
        r.raise_for_status()
        data = r.json()["esearchresult"]
        count = int(data.get("count", 0))
        ids   = data.get("idlist", [])
    except Exception:
        return {"sra_run_count": None, "sra_platforms": "", "sra_strategies": "",
                "sra_total_bases": None, "sra_date_range": ""}

    if count == 0 or not ids:
        return {"sra_run_count": 0, "sra_platforms": "", "sra_strategies": "",
                "sra_total_bases": None, "sra_date_range": ""}

    # Fetch summaries for up to 200 runs to collect platform/strategy stats
    uid_str = ",".join(ids[:200])
    sumurl  = (f"{NCBI_BASE}/esummary.fcgi"
               f"?db=sra&id={uid_str}&retmode=json")
    try:
        sr = NCBI_SESSION.get(sumurl)
        sr.raise_for_status()
        result = sr.json().get("result", {})
    except Exception:
        return {"sra_run_count": count, "sra_platforms": "", "sra_strategies": "",
                "sra_total_bases": None, "sra_date_range": ""}

    platforms  = set()
    strategies = set()
    total_bases = 0
    dates = []

    for uid, entry in result.items():
        if uid == "uids":
            continue
        expxml = entry.get("expxml", "")
        # platform
        plat_m = re.search(r'<Platform[^>]*instrument_model="([^"]*)"[^>]*>([^<]*)<', expxml)
        if plat_m:
            platforms.add(plat_m.group(2).strip() or plat_m.group(1).strip())
        # bases
        bases_m = re.search(r'total_bases="(\d+)"', expxml)
        if bases_m:
            total_bases += int(bases_m.group(1))
        # strategy
        strat_m = re.search(r'<LIBRARY_STRATEGY>([^<]+)</LIBRARY_STRATEGY>', expxml)
        if strat_m:
            strategies.add(strat_m.group(1).strip())
        # date
        cd = entry.get("createdate", "")
        if cd:
            dates.append(cd[:4])

    date_range = ""
    if dates:
        date_range = f"{min(dates)}-{max(dates)}"

    return {
        "sra_run_count":   count,
        "sra_platforms":   "; ".join(sorted(platforms)),
        "sra_strategies":  "; ".join(sorted(strategies)),
        "sra_total_bases": total_bases if total_bases else None,
        "sra_date_range":  date_range,
    }


def fetch_sra_stats(orgs: pd.DataFrame, limit: int | None = None) -> pd.DataFrame:
    """Fetch SRA stats for all organisms, skipping those already cached."""
    if SRA_CSV.exists():
        cached = pd.read_csv(SRA_CSV, dtype=str)
        done   = set(cached["organism"].tolist())
    else:
        cached = pd.DataFrame()
        done   = set()

    todo = orgs[~orgs["organism"].isin(done)].copy()
    if limit:
        todo = todo.head(limit)

    if todo.empty:
        print("  SRA: all organisms already cached.")
        return cached

    print(f"  SRA: fetching {len(todo)} organisms "
          f"(~{len(todo) * 0.5:.0f}s at 2.5 req/s)...")

    rows = []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        stats = _sra_search(row["organism"])
        stats["organism"] = row["organism"]
        rows.append(stats)
        if i % 50 == 0:
            print(f"    {i}/{len(todo)} done")

    new_df = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(SRA_CSV, index=False)
    print(f"  SRA: saved {len(combined)} rows to {SRA_CSV.name}")
    return combined


# ---------------------------------------------------------------------------
# BacDive fetch (public HTML scraping)
# ---------------------------------------------------------------------------
def _bacdive_search_ids(organism: str) -> list[str]:
    """Return list of BacDive strain IDs for an organism via the public search page."""
    url = f"https://bacdive.dsmz.de/?search={requests.utils.quote(organism)}"
    try:
        r = BACDIVE_SESSION.get(url)
        r.raise_for_status()
        ids = re.findall(r'href="/strain/(\d+)"', r.text)
        return list(dict.fromkeys(ids))[:3]
    except Exception:
        return []


def _extract_table_values(soup: BeautifulSoup, section_text: str) -> list[str]:
    """Find a table following a header containing section_text, return all cell values."""
    for h in soup.find_all(["h3", "h4"]):
        if section_text.lower() in h.get_text().lower():
            tbl = h.find_next("table")
            if tbl:
                return [td.get_text(strip=True) for td in tbl.find_all("td")]
    return []


def _bacdive_strain(strain_id: str) -> dict:
    """Scrape a BacDive strain page and return structured data."""
    url = f"https://bacdive.dsmz.de/strain/{strain_id}"
    try:
        r = BACDIVE_SESSION.get(url)
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "html.parser")
    except Exception:
        return {}

    data: dict = {"bacdive_strain_id": strain_id, "bacdive_url": url}

    # Taxonomy / name
    name_tag = soup.find("h1")
    if name_tag:
        data["bacdive_name"] = name_tag.get_text(strip=True)

    # Culture temperature
    temp_vals = _extract_table_values(soup, "Culture temp")
    temps = [v for v in temp_vals if re.match(r"^\d+\.?\d*$", v)]
    if temps:
        data["bacdive_temp_c"] = "; ".join(temps)

    # pH
    ph_vals = _extract_table_values(soup, "pH")
    phs = [v for v in ph_vals if re.match(r"^\d+\.?\d*$", v)]
    if phs:
        data["bacdive_ph"] = "; ".join(phs)

    # Oxygen tolerance
    oxy = _extract_table_values(soup, "Oxygen")
    oxy_clean = [v for v in oxy if len(v) > 2 and not v.startswith("@")]
    if oxy_clean:
        data["bacdive_oxygen"] = "; ".join(oxy_clean[:3])

    # Gram stain
    gram = _extract_table_values(soup, "Gram")
    gram_clean = [v for v in gram if v and not v.startswith("@")]
    if gram_clean:
        data["bacdive_gram"] = gram_clean[0]

    # Cell morphology
    morph = _extract_table_values(soup, "morphol")
    morph_clean = [v for v in morph if len(v) > 2 and not v.startswith("@")]
    if morph_clean:
        data["bacdive_morphology"] = "; ".join(morph_clean[:4])

    # Motility
    motility = _extract_table_values(soup, "Motil")
    mot_clean = [v for v in motility if v and not v.startswith("@")]
    if mot_clean:
        data["bacdive_motility"] = mot_clean[0]

    # Isolation
    isol_vals = _extract_table_values(soup, "Isolation")
    isol_clean = [v for v in isol_vals if len(v) > 2 and not v.startswith("@")]
    if isol_clean:
        data["bacdive_isolation"] = "; ".join(isol_clean[:6])

    return data


def fetch_bacdive_data(orgs: pd.DataFrame, limit: int | None = None) -> pd.DataFrame:
    """Fetch BacDive data for organisms not already cached."""
    if BACDIVE_CSV.exists():
        cached = pd.read_csv(BACDIVE_CSV, dtype=str)
        done   = set(cached["organism"].tolist())
    else:
        cached = pd.DataFrame()
        done   = set()

    # Prioritise organisms with more entries in PlasticDB
    todo = orgs[~orgs["organism"].isin(done)].copy()
    todo = todo.sort_values("n_entries", ascending=False)
    if limit:
        todo = todo.head(limit)

    if todo.empty:
        print("  BacDive: all organisms already cached.")
        return cached

    print(f"  BacDive: fetching {len(todo)} organisms "
          f"(~{len(todo) * 1.5:.0f}s at 1.5 req/s)...")

    rows = []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        org_name = row["organism"]
        strain_ids = _bacdive_search_ids(org_name)
        if not strain_ids:
            rows.append({"organism": org_name, "bacdive_found": "No"})
        else:
            strain_data = _bacdive_strain(strain_ids[0])
            strain_data["organism"]      = org_name
            strain_data["bacdive_found"] = "Yes"
            rows.append(strain_data)
        if i % 25 == 0:
            print(f"    {i}/{len(todo)} done")

    new_df   = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(BACDIVE_CSV, index=False)
    print(f"  BacDive: saved {len(combined)} rows to {BACDIVE_CSV.name}")
    return combined


# ---------------------------------------------------------------------------
# HTML page generation
# ---------------------------------------------------------------------------
def _slug(name: str) -> str:
    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       background: #f8f9fa; color: #1a1a2e; line-height: 1.6; }
.header { background: #1a1a2e; color: #fff; padding: 24px 40px; }
.header h1 { font-size: 1.9rem; font-weight: 700; letter-spacing: -0.5px; }
.header .meta { color: #9fb3cc; font-size: 0.9rem; margin-top: 6px; }
.breadcrumb { padding: 12px 40px; background: #eef1f5; font-size: 0.85rem; }
.breadcrumb a { color: #3a7bd5; text-decoration: none; }
.container { max-width: 1100px; margin: 0 auto; padding: 32px 40px; }
.section { background: #fff; border-radius: 8px; border: 1px solid #dde3ec;
           margin-bottom: 24px; overflow: hidden; }
.section-header { background: #f0f4fb; border-bottom: 1px solid #dde3ec;
                  padding: 14px 22px; font-weight: 600; font-size: 1rem;
                  color: #1a1a2e; }
.section-body { padding: 20px 22px; }
.grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
.kv-row { display: flex; gap: 12px; margin-bottom: 10px; font-size: 0.92rem; }
.kv-label { color: #6b7a99; min-width: 180px; flex-shrink: 0; }
.kv-value { font-weight: 500; }
.pill { display: inline-block; background: #e8f0fe; color: #3a5cc0;
        border-radius: 4px; padding: 2px 8px; font-size: 0.8rem;
        margin: 2px 2px 2px 0; }
.badge { display: inline-block; border-radius: 4px; padding: 2px 8px;
         font-size: 0.78rem; font-weight: 600; margin-right: 6px; }
.badge-yes { background: #d4edda; color: #155724; }
.badge-no  { background: #f8d7da; color: #721c24; }
.badge-na  { background: #e2e3e5; color: #383d41; }
table { width: 100%; border-collapse: collapse; font-size: 0.88rem; }
th { background: #f0f4fb; padding: 10px 14px; text-align: left;
     font-weight: 600; border-bottom: 2px solid #dde3ec; }
td { padding: 9px 14px; border-bottom: 1px solid #edf0f5; vertical-align: top; }
tr:last-child td { border-bottom: none; }
.stat-box { background: #f0f4fb; border-radius: 6px; padding: 16px;
            text-align: center; }
.stat-val { font-size: 1.8rem; font-weight: 700; color: #3a5cc0; }
.stat-lbl { font-size: 0.8rem; color: #6b7a99; margin-top: 4px; }
.no-data { color: #9aa3b5; font-style: italic; font-size: 0.9rem; }
a { color: #3a7bd5; }
@media (max-width: 680px) { .grid2 { grid-template-columns: 1fr; }
  .container { padding: 16px; } }
"""


def _badge(val) -> str:
    if val is True or str(val).lower() in ("true", "yes", "1"):
        return '<span class="badge badge-yes">Yes</span>'
    if val is False or str(val).lower() in ("false", "no", "0"):
        return '<span class="badge badge-no">No</span>'
    return '<span class="badge badge-na">N/A</span>'


def _fmt_bases(bases) -> str:
    try:
        b = float(bases)
        if b >= 1e12:
            return f"{b/1e12:.2f} Tbp"
        if b >= 1e9:
            return f"{b/1e9:.2f} Gbp"
        if b >= 1e6:
            return f"{b/1e6:.2f} Mbp"
        return f"{b:.0f} bp"
    except (TypeError, ValueError):
        return "N/A"


def _kv(label: str, value, pill: bool = False) -> str:
    if pill and value:
        pills = "".join(f'<span class="pill">{v.strip()}</span>'
                        for v in str(value).split(";") if v.strip())
        return (f'<div class="kv-row"><span class="kv-label">{label}</span>'
                f'<span class="kv-value">{pills}</span></div>')
    return (f'<div class="kv-row"><span class="kv-label">{label}</span>'
            f'<span class="kv-value">{value if value else "N/A"}</span></div>')


def _plastics_pills(plastics_list) -> str:
    if not plastics_list:
        return "<span class='no-data'>none recorded</span>"
    return "".join(f'<span class="pill">{p}</span>' for p in plastics_list)


def generate_organism_page(
    org_row: pd.Series,
    full_df: pd.DataFrame,
    sra_row: pd.Series | None,
    bd_row: pd.Series | None,
) -> str:
    name    = org_row["organism"]
    slug_id = _slug(name)

    # Extract entries for this organism
    ent = full_df[full_df["organism"] == name].copy()
    entries_html = ""
    if not ent.empty:
        rows_html = ""
        for _, e in ent.iterrows():
            seq_badge  = _badge(e.get("has_sequence", False))
            genb_badge = _badge(e.get("has_genbank",  False))
            enz_name   = e.get("enzyme_name", "") or ""
            doi_val    = e.get("doi", "") or ""
            doi_link   = (f'<a href="https://doi.org/{doi_val}" target="_blank">'
                          f'{doi_val[:40]}...</a>' if doi_val and len(doi_val) > 5 else doi_val)
            rows_html += (
                f"<tr>"
                f"<td>{e.get('plastic','')}</td>"
                f"<td>{int(e['year']) if pd.notna(e.get('year')) else 'N/A'}</td>"
                f"<td>{enz_name[:50] or 'N/A'}</td>"
                f"<td>{seq_badge}</td>"
                f"<td>{genb_badge}</td>"
                f"<td>{doi_link}</td>"
                f"</tr>"
            )
        entries_html = f"""
        <table>
          <thead><tr><th>Plastic</th><th>Year</th><th>Enzyme</th>
            <th>Sequence</th><th>GenBank</th><th>DOI</th></tr></thead>
          <tbody>{rows_html}</tbody>
        </table>"""

    # SRA section
    if sra_row is not None:
        run_count = sra_row.get("sra_run_count", "N/A")
        try:
            run_count = int(float(run_count))
        except (TypeError, ValueError):
            run_count = "N/A"
        sra_html = f"""
        <div class="grid2" style="margin-bottom:20px">
          <div class="stat-box"><div class="stat-val">{run_count}</div>
            <div class="stat-lbl">SRA runs deposited</div></div>
          <div class="stat-box"><div class="stat-val">{_fmt_bases(sra_row.get('sra_total_bases'))}</div>
            <div class="stat-lbl">Total sequenced bases (sampled)</div></div>
        </div>
        {_kv('Sequencing platforms', sra_row.get('sra_platforms','') or 'None in sample', pill=True)}
        {_kv('Library strategies',  sra_row.get('sra_strategies','') or 'None in sample', pill=True)}
        {_kv('Deposit year range',  sra_row.get('sra_date_range','') or 'N/A')}
        <p style="font-size:0.8rem;color:#9aa3b5;margin-top:10px">
          Counts via NCBI SRA E-utilities esearch/esummary.
          Platform and strategy data sampled from up to 200 runs.
        </p>"""
    else:
        sra_html = '<p class="no-data">SRA data not fetched for this organism.</p>'

    # BacDive section
    if bd_row is not None and str(bd_row.get("bacdive_found", "No")).lower() == "yes":
        bd_url = bd_row.get("bacdive_url", "")
        bd_link = (f'<a href="{bd_url}" target="_blank">'
                   f'BacDive strain {bd_row.get("bacdive_strain_id","")}</a>'
                   if bd_url else "N/A")
        bd_html = f"""
        {_kv('BacDive record', bd_link)}
        {_kv('Culture temperature (deg C)', bd_row.get('bacdive_temp_c','') or 'Not recorded')}
        {_kv('pH',                          bd_row.get('bacdive_ph','')     or 'Not recorded')}
        {_kv('Oxygen tolerance',            bd_row.get('bacdive_oxygen','') or 'Not recorded')}
        {_kv('Gram stain',                  bd_row.get('bacdive_gram','')   or 'Not recorded')}
        {_kv('Cell morphology',             bd_row.get('bacdive_morphology','') or 'Not recorded')}
        {_kv('Motility',                    bd_row.get('bacdive_motility','')   or 'Not recorded')}
        {_kv('Isolation source',            bd_row.get('bacdive_isolation','')  or 'Not recorded')}
        <p style="font-size:0.8rem;color:#9aa3b5;margin-top:10px">
          Data scraped from the public BacDive strain page.
        </p>"""
    elif bd_row is not None:
        bd_html = '<p class="no-data">Organism not found in BacDive public database.</p>'
    else:
        bd_html = '<p class="no-data">BacDive data not fetched for this organism.</p>'

    # Plastics list
    plastics_list = org_row.get("plastics", [])
    if isinstance(plastics_list, str):
        import ast
        try:
            plastics_list = ast.literal_eval(plastics_list)
        except Exception:
            plastics_list = [p.strip() for p in plastics_list.split(",")]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{name} | Organism Profile</title>
  <style>{CSS}</style>
</head>
<body>
<div class="header">
  <h1>{name}</h1>
  <div class="meta">
    Tax ID: {org_row.get('tax_id','N/A')} &nbsp;|&nbsp;
    {org_row.get('n_plastics', 0)} plastic type(s) &nbsp;|&nbsp;
    {org_row.get('n_entries', 0)} PlasticDB entries &nbsp;|&nbsp;
    First recorded: {int(org_row['first_year']) if pd.notna(org_row.get('first_year')) else 'N/A'} &nbsp;|&nbsp;
    Last recorded: {int(org_row['last_year']) if pd.notna(org_row.get('last_year')) else 'N/A'}
  </div>
</div>
<div class="breadcrumb">
  <a href="index.html">All organisms</a> / {name}
</div>
<div class="container">

  <!-- PlasticDB summary -->
  <div class="section">
    <div class="section-header">PlasticDB Summary</div>
    <div class="section-body">
      <div class="grid2">
        <div>
          {_kv('Plastics degraded', _plastics_pills(plastics_list))}
          {_kv('Thermophilic flag', org_row.get('thermophilic','') or 'Not recorded')}
          {_kv('Has linked sequence', _badge(org_row.get('has_sequence')))}
          {_kv('Has named enzyme',    _badge(org_row.get('has_enzyme')))}
          {_kv('Has GenBank ID',      _badge(org_row.get('has_genbank')))}
        </div>
        <div>
          {_kv('Isolation environments', org_row.get('isolation_envs','') or 'Not recorded', pill=True)}
          {_kv('Isolation locations',    org_row.get('isolation_locs','')  or 'Not recorded', pill=True)}
        </div>
      </div>
    </div>
  </div>

  <!-- Database entries table -->
  <div class="section">
    <div class="section-header">PlasticDB Entries ({org_row.get('n_entries',0)} total)</div>
    <div class="section-body">{entries_html}</div>
  </div>

  <!-- SRA stats -->
  <div class="section">
    <div class="section-header">NCBI SRA Sequencing Data</div>
    <div class="section-body">{sra_html}</div>
  </div>

  <!-- BacDive -->
  <div class="section">
    <div class="section-header">BacDive Physiological Data</div>
    <div class="section-body">{bd_html}</div>
  </div>

</div>
</body>
</html>"""
    return html


def generate_index(orgs: pd.DataFrame, sra_df: pd.DataFrame, bd_df: pd.DataFrame) -> str:
    sra_map = dict(zip(sra_df["organism"], sra_df.get("sra_run_count", pd.Series(dtype=str)))) \
        if not sra_df.empty else {}
    bd_map  = dict(zip(bd_df["organism"],  bd_df.get("bacdive_found",  pd.Series(dtype=str)))) \
        if not bd_df.empty else {}

    rows_html = ""
    for _, row in orgs.sort_values("n_entries", ascending=False).iterrows():
        org  = row["organism"]
        s    = _slug(org)
        rc   = sra_map.get(org, "")
        try:
            rc = int(float(rc))
        except (TypeError, ValueError):
            rc = "N/A"
        bd_found = str(bd_map.get(org, "")).lower() == "yes"
        bd_badge = _badge(bd_found)

        plastics_list = row.get("plastics", [])
        if isinstance(plastics_list, str):
            import ast
            try:
                plastics_list = ast.literal_eval(plastics_list)
            except Exception:
                plastics_list = [p.strip() for p in plastics_list.split(",")]

        pills = "".join(f'<span class="pill">{p}</span>' for p in plastics_list[:5])
        if len(plastics_list) > 5:
            pills += f'<span class="pill">+{len(plastics_list)-5}</span>'

        rows_html += (
            f'<tr>'
            f'<td><a href="{s}.html">{org}</a></td>'
            f'<td style="text-align:center">{row.get("n_entries",0)}</td>'
            f'<td>{pills}</td>'
            f'<td style="text-align:center">{rc}</td>'
            f'<td style="text-align:center">{bd_badge}</td>'
            f'<td style="text-align:center">{int(row["first_year"]) if pd.notna(row.get("first_year")) else "N/A"}</td>'
            f'</tr>'
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Organism Profiles | PlasticDB + SRA + BacDive</title>
  <style>{CSS}
  #search-input {{ width:100%;padding:10px 14px;font-size:1rem;border:1px solid #ccd3e0;
    border-radius:6px;margin-bottom:16px;outline:none; }}
  </style>
</head>
<body>
<div class="header">
  <h1>Organism Profiles</h1>
  <div class="meta">
    {len(orgs)} organisms from PlasticDB, cross-referenced with NCBI SRA and BacDive.
    Data fetched from live APIs. No values fabricated.
  </div>
</div>
<div class="container">
  <div class="section">
    <div class="section-header">All Organisms</div>
    <div class="section-body">
      <input id="search-input" placeholder="Filter by organism name..." oninput="filterTable(this.value)">
      <table id="org-table">
        <thead><tr>
          <th>Organism</th>
          <th style="text-align:center">PlasticDB entries</th>
          <th>Plastics degraded</th>
          <th style="text-align:center">SRA runs</th>
          <th style="text-align:center">BacDive</th>
          <th style="text-align:center">First year</th>
        </tr></thead>
        <tbody>{rows_html}</tbody>
      </table>
    </div>
  </div>
</div>
<script>
function filterTable(q) {{
  q = q.toLowerCase();
  document.querySelectorAll('#org-table tbody tr').forEach(r => {{
    r.style.display = r.cells[0].textContent.toLowerCase().includes(q) ? '' : 'none';
  }});
}}
</script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sra-only",   action="store_true")
    parser.add_argument("--pages-only", action="store_true")
    parser.add_argument("--limit",      type=int, default=None)
    args = parser.parse_args()

    print("Loading PlasticDB...")
    full_df, orgs = load_plasticdb()
    print(f"  {len(full_df):,} entries, {len(orgs):,} unique organisms.")

    if args.limit:
        orgs = orgs.sort_values("n_entries", ascending=False).head(args.limit)
        print(f"  Limiting to {args.limit} organisms.")

    if not args.pages_only:
        print("\nFetching SRA stats...")
        sra_df = fetch_sra_stats(orgs, limit=args.limit)

        if not args.sra_only:
            print("\nFetching BacDive data...")
            bd_df = fetch_bacdive_data(orgs, limit=args.limit)
        else:
            bd_df = pd.read_csv(BACDIVE_CSV) if BACDIVE_CSV.exists() else pd.DataFrame()
    else:
        sra_df = pd.read_csv(SRA_CSV)     if SRA_CSV.exists()     else pd.DataFrame()
        bd_df  = pd.read_csv(BACDIVE_CSV) if BACDIVE_CSV.exists() else pd.DataFrame()

    print(f"\nGenerating {len(orgs):,} organism pages...")
    sra_idx = sra_df.set_index("organism") if not sra_df.empty else pd.DataFrame()
    bd_idx  = bd_df.set_index("organism")  if not bd_df.empty  else pd.DataFrame()

    for i, (_, row) in enumerate(orgs.iterrows(), 1):
        name = row["organism"]
        s    = _slug(name)
        sra_row = sra_idx.loc[name] if (not sra_idx.empty and name in sra_idx.index) else None
        bd_row  = bd_idx.loc[name]  if (not bd_idx.empty  and name in bd_idx.index)  else None
        html = generate_organism_page(row, full_df, sra_row, bd_row)
        (PAGES_DIR / f"{s}.html").write_text(html, encoding="utf-8")
        if i % 100 == 0:
            print(f"  {i}/{len(orgs)} pages written")

    print(f"  All {len(orgs)} pages written.")

    print("Generating index.html...")
    index_html = generate_index(orgs, sra_df, bd_df)
    (PAGES_DIR / "index.html").write_text(index_html, encoding="utf-8")
    print(f"  index.html written.")
    print(f"\nDone. Open organism-profiles/pages/index.html to browse.")


if __name__ == "__main__":
    main()
