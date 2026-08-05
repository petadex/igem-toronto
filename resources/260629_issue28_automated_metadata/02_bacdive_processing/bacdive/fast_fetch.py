"""
Optimized data fetcher that processes one batch per run.
Each run picks up where the cache left off.

Usage:
    python fast_fetch.py --mode sra   --batch 300   # fetch SRA counts for next 300 organisms
    python fast_fetch.py --mode bd    --batch 50    # fetch BacDive for next 50 organisms
    python fast_fetch.py --mode pages               # generate HTML pages from cache
"""

import argparse
import re
import time
import unicodedata
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

ROOT        = Path(__file__).parent
DATA_DIR    = ROOT / "data"
PAGES_DIR   = ROOT / "pages"
PDB_TSV     = ROOT.parent / "plastic-biodegradation-analysis" / "data" / "plasticdb_microorganisms.tsv"
SRA_CSV     = DATA_DIR / "sra_stats.csv"
BACDIVE_CSV = DATA_DIR / "bacdive_data.csv"
DATA_DIR.mkdir(exist_ok=True)
PAGES_DIR.mkdir(exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "igem-toronto-research/1.0 (academic)"
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


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
        )
        .reset_index()
        .sort_values("n_entries", ascending=False)
    )
    return df, orgs


# ---------------------------------------------------------------------------
# SRA: count-only fetch (1 request per organism)
# ---------------------------------------------------------------------------
def fetch_sra_count(organism: str) -> dict:
    term = f'"{organism}"[Organism]'
    url  = (f"{NCBI_BASE}/esearch.fcgi"
            f"?db=sra&term={requests.utils.quote(term)}&retmax=5&retmode=json")
    time.sleep(0.34)
    try:
        r = session.get(url, timeout=20)
        r.raise_for_status()
        data = r.json()["esearchresult"]
        count = int(data.get("count", 0))
        ids   = data.get("idlist", [])
    except Exception as e:
        return {"organism": organism, "sra_run_count": None,
                "sra_platforms": "", "sra_strategies": "",
                "sra_total_bases": None, "sra_date_range": ""}

    if count == 0 or not ids:
        return {"organism": organism, "sra_run_count": 0,
                "sra_platforms": "", "sra_strategies": "",
                "sra_total_bases": None, "sra_date_range": ""}

    # Fetch summary for top 5 runs to get platform/strategy
    time.sleep(0.34)
    uid_str = ",".join(ids[:5])
    sumurl  = f"{NCBI_BASE}/esummary.fcgi?db=sra&id={uid_str}&retmode=json"
    try:
        sr = session.get(sumurl, timeout=20)
        sr.raise_for_status()
        result = sr.json().get("result", {})
    except Exception:
        return {"organism": organism, "sra_run_count": count,
                "sra_platforms": "", "sra_strategies": "",
                "sra_total_bases": None, "sra_date_range": ""}

    platforms, strategies, dates = set(), set(), []
    total_bases = 0
    for uid, entry in result.items():
        if uid == "uids":
            continue
        x = entry.get("expxml", "")
        m = re.search(r"<Platform[^>]*>([^<]+)<", x)
        if m:
            platforms.add(m.group(1).strip())
        m = re.search(r'total_bases="(\d+)"', x)
        if m:
            total_bases += int(m.group(1))
        m = re.search(r"<LIBRARY_STRATEGY>([^<]+)</LIBRARY_STRATEGY>", x)
        if m:
            strategies.add(m.group(1).strip())
        cd = entry.get("createdate", "")[:4]
        if cd:
            dates.append(cd)

    return {
        "organism":        organism,
        "sra_run_count":   count,
        "sra_platforms":   "; ".join(sorted(platforms)),
        "sra_strategies":  "; ".join(sorted(strategies)),
        "sra_total_bases": total_bases or None,
        "sra_date_range":  f"{min(dates)}-{max(dates)}" if dates else "",
    }


def run_sra_batch(orgs: pd.DataFrame, batch: int):
    cached = pd.read_csv(SRA_CSV, dtype=str) if SRA_CSV.exists() else pd.DataFrame()
    done   = set(cached["organism"].tolist()) if not cached.empty else set()
    todo   = orgs[~orgs["organism"].isin(done)].head(batch)

    if todo.empty:
        print(f"SRA: all {len(cached)} organisms already cached.")
        return

    print(f"SRA: fetching {len(todo)} organisms (cache has {len(done)})...")
    rows = []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        result = fetch_sra_count(row["organism"])
        rows.append(result)
        if i % 30 == 0:
            print(f"  {i}/{len(todo)}")

    new_df   = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(SRA_CSV, index=False)
    total = len(combined)
    found = int((pd.to_numeric(combined["sra_run_count"], errors="coerce") > 0).sum())
    print(f"SRA: saved. Total cached: {total}. With runs: {found}.")


# ---------------------------------------------------------------------------
# BacDive fetch
# ---------------------------------------------------------------------------
def _search_bd(organism: str) -> list:
    time.sleep(0.5)
    url = f"https://bacdive.dsmz.de/search?search={requests.utils.quote(organism)}"
    try:
        r = session.get(url, timeout=20)
        r.raise_for_status()
        # Direct redirect to strain page (exact match)
        if "/strain/" in r.url:
            m = re.search(r"/strain/(\d+)", r.url)
            if m:
                return [m.group(1)]
        # Multiple results page: extract strain IDs from HTML links
        ids = re.findall(r"/strain/(\d+)", r.text)
        return list(dict.fromkeys(ids))[:2]
    except Exception:
        return []


def _scrape_bd(strain_id: str) -> dict:
    time.sleep(0.5)
    url = f"https://bacdive.dsmz.de/strain/{strain_id}"
    try:
        r = session.get(url, timeout=20)
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "html.parser")
    except Exception:
        return {}

    data = {"bacdive_strain_id": strain_id, "bacdive_url": url}

    def tbl(section):
        for h in soup.find_all(["h3", "h4"]):
            if section.lower() in h.get_text().lower():
                t = h.find_next("table")
                if t:
                    return [td.get_text(strip=True) for td in t.find_all("td")]
        return []

    temps = [v for v in tbl("Culture temp") if re.match(r"^\d+\.?\d*$", v)]
    if temps:
        data["bacdive_temp_c"] = "; ".join(temps)
    phs = [v for v in tbl("pH") if re.match(r"^\d+\.?\d*$", v)]
    if phs:
        data["bacdive_ph"] = "; ".join(phs)
    oxy = [v for v in tbl("Oxygen") if len(v) > 2 and not v.startswith("@")]
    if oxy:
        data["bacdive_oxygen"] = "; ".join(oxy[:3])
    gram = [v for v in tbl("Gram") if v and not v.startswith("@")]
    if gram:
        data["bacdive_gram"] = gram[0]
    morph = [v for v in tbl("morphol") if len(v) > 2 and not v.startswith("@")]
    if morph:
        data["bacdive_morphology"] = "; ".join(morph[:4])
    mot = [v for v in tbl("Motil") if v and not v.startswith("@")]
    if mot:
        data["bacdive_motility"] = mot[0]
    isol = [v for v in tbl("Isolation") if len(v) > 2 and not v.startswith("@")]
    if isol:
        data["bacdive_isolation"] = "; ".join(isol[:6])
    return data


def run_bacdive_batch(orgs: pd.DataFrame, batch: int):
    cached = pd.read_csv(BACDIVE_CSV, dtype=str) if BACDIVE_CSV.exists() else pd.DataFrame()
    done   = set(cached["organism"].tolist()) if not cached.empty else set()
    todo   = orgs[~orgs["organism"].isin(done)].head(batch)

    if todo.empty:
        print(f"BacDive: all {len(cached)} organisms already cached.")
        return

    print(f"BacDive: fetching {len(todo)} organisms (cache has {len(done)})...")
    rows = []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        org = row["organism"]
        ids = _search_bd(org)
        if not ids:
            rows.append({"organism": org, "bacdive_found": "No"})
        else:
            d = _scrape_bd(ids[0])
            d.update({"organism": org, "bacdive_found": "Yes"})
            rows.append(d)
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")

    new_df   = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(BACDIVE_CSV, index=False)
    found = (combined.get("bacdive_found", pd.Series()) == "Yes").sum()
    print(f"BacDive: saved. Total cached: {len(combined)}. Found in BacDive: {found}.")


# ---------------------------------------------------------------------------
# Page generation
# ---------------------------------------------------------------------------
CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       background: #f8f9fa; color: #1a1a2e; line-height: 1.6; }
.header { background: #1a1a2e; color: #fff; padding: 24px 40px; }
.header h1 { font-size: 1.9rem; font-weight: 700; }
.header .meta { color: #9fb3cc; font-size: 0.9rem; margin-top: 6px; }
.breadcrumb { padding: 12px 40px; background: #eef1f5; font-size: 0.85rem; }
.breadcrumb a { color: #3a7bd5; text-decoration: none; }
.container { max-width: 1100px; margin: 0 auto; padding: 32px 40px; }
.section { background: #fff; border-radius: 8px; border: 1px solid #dde3ec;
           margin-bottom: 24px; overflow: hidden; }
.section-header { background: #f0f4fb; border-bottom: 1px solid #dde3ec;
                  padding: 14px 22px; font-weight: 600; font-size: 1rem; }
.section-body { padding: 20px 22px; }
.grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
.kv-row { display: flex; gap: 12px; margin-bottom: 10px; font-size: 0.92rem; }
.kv-label { color: #6b7a99; min-width: 180px; flex-shrink: 0; }
.kv-value { font-weight: 500; }
.pill { display: inline-block; background: #e8f0fe; color: #3a5cc0;
        border-radius: 4px; padding: 2px 8px; font-size: 0.8rem; margin: 2px; }
.badge { display: inline-block; border-radius: 4px; padding: 2px 8px;
         font-size: 0.78rem; font-weight: 600; }
.badge-yes { background: #d4edda; color: #155724; }
.badge-no  { background: #f8d7da; color: #721c24; }
.badge-na  { background: #e2e3e5; color: #383d41; }
table { width: 100%; border-collapse: collapse; font-size: 0.88rem; }
th { background: #f0f4fb; padding: 10px 14px; text-align: left;
     font-weight: 600; border-bottom: 2px solid #dde3ec; }
td { padding: 9px 14px; border-bottom: 1px solid #edf0f5; vertical-align: top; }
.stat-box { background: #f0f4fb; border-radius: 6px; padding: 16px; text-align: center; }
.stat-val { font-size: 1.8rem; font-weight: 700; color: #3a5cc0; }
.stat-lbl { font-size: 0.8rem; color: #6b7a99; margin-top: 4px; }
.no-data  { color: #9aa3b5; font-style: italic; font-size: 0.9rem; }
.note     { font-size: 0.8rem; color: #9aa3b5; margin-top: 10px; }
a { color: #3a7bd5; }
#search-input { width: 100%; padding: 10px 14px; font-size: 1rem;
  border: 1px solid #ccd3e0; border-radius: 6px; margin-bottom: 16px; }
"""


def _slug(name):
    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def _badge(val):
    s = str(val).lower()
    if s in ("true","yes","1"): return '<span class="badge badge-yes">Yes</span>'
    if s in ("false","no","0"): return '<span class="badge badge-no">No</span>'
    return '<span class="badge badge-na">N/A</span>'


def _kv(label, value, pill=False):
    if pill and value:
        pills = "".join(f'<span class="pill">{v.strip()}</span>'
                        for v in str(value).split(";") if v.strip())
        return (f'<div class="kv-row"><span class="kv-label">{label}</span>'
                f'<span class="kv-value">{pills}</span></div>')
    return (f'<div class="kv-row"><span class="kv-label">{label}</span>'
            f'<span class="kv-value">{value if value else "N/A"}</span></div>')


def _fmt_bases(b):
    try:
        b = float(b)
        if b >= 1e12: return f"{b/1e12:.2f} Tbp"
        if b >= 1e9:  return f"{b/1e9:.2f} Gbp"
        if b >= 1e6:  return f"{b/1e6:.2f} Mbp"
        return f"{b:.0f} bp"
    except Exception:
        return "N/A"


def _make_page(org_row, full_df, sra_row, bd_row):
    import ast
    name = org_row["organism"]
    ent  = full_df[full_df["organism"] == name]

    rows_html = ""
    for _, e in ent.iterrows():
        doi = str(e.get("doi","") or "")
        doi_link = (f'<a href="https://doi.org/{doi}" target="_blank">{doi[:35]}</a>'
                    if len(doi) > 5 else doi)
        yr = int(e["year"]) if pd.notna(e.get("year")) else "N/A"
        rows_html += (
            f"<tr><td>{e.get('plastic','')}</td><td>{yr}</td>"
            f"<td>{str(e.get('enzyme_name','') or '')[:50] or 'N/A'}</td>"
            f"<td>{_badge(e.get('has_sequence',False))}</td>"
            f"<td>{_badge(e.get('has_genbank',False))}</td>"
            f"<td style='font-size:0.8rem'>{doi_link}</td></tr>"
        )
    entries_html = (
        "<table><thead><tr><th>Plastic</th><th>Year</th><th>Enzyme</th>"
        "<th>Sequence</th><th>GenBank</th><th>DOI</th></tr></thead>"
        f"<tbody>{rows_html}</tbody></table>"
    )

    # SRA
    if sra_row is not None:
        try: rc = int(float(sra_row.get("sra_run_count","")))
        except: rc = "N/A"
        sra_html = (
            f'<div class="grid2" style="margin-bottom:18px">'
            f'<div class="stat-box"><div class="stat-val">{rc}</div>'
            f'<div class="stat-lbl">SRA runs deposited</div></div>'
            f'<div class="stat-box"><div class="stat-val">{_fmt_bases(sra_row.get("sra_total_bases"))}</div>'
            f'<div class="stat-lbl">Bases (sampled from top 5 runs)</div></div></div>'
            + _kv("Sequencing platforms", sra_row.get("sra_platforms","") or "None in sample", pill=True)
            + _kv("Library strategies",  sra_row.get("sra_strategies","") or "None in sample", pill=True)
            + _kv("Deposit year range",  sra_row.get("sra_date_range","") or "N/A")
            + '<p class="note">Source: NCBI SRA E-utilities (esearch + esummary). '
              'Run count is exact. Platform/strategy sampled from top 5 runs.</p>'
        )
    else:
        sra_html = '<p class="no-data">SRA data not yet fetched for this organism.</p>'

    # BacDive
    if bd_row is not None and str(bd_row.get("bacdive_found","")).lower() == "yes":
        bd_url  = bd_row.get("bacdive_url","")
        bd_link = (f'<a href="{bd_url}" target="_blank">BacDive strain '
                   f'{bd_row.get("bacdive_strain_id","")}</a>' if bd_url else "N/A")
        bd_html = (
            _kv("BacDive record",      bd_link)
            + _kv("Culture temp (C)", bd_row.get("bacdive_temp_c","") or "Not recorded")
            + _kv("pH",               bd_row.get("bacdive_ph","")     or "Not recorded")
            + _kv("Oxygen tolerance", bd_row.get("bacdive_oxygen","") or "Not recorded")
            + _kv("Gram stain",       bd_row.get("bacdive_gram","")   or "Not recorded")
            + _kv("Morphology",       bd_row.get("bacdive_morphology","") or "Not recorded")
            + _kv("Motility",         bd_row.get("bacdive_motility","")   or "Not recorded")
            + _kv("Isolation source", bd_row.get("bacdive_isolation","")  or "Not recorded")
            + '<p class="note">Source: BacDive public strain page (bacdive.dsmz.de).</p>'
        )
    elif bd_row is not None:
        bd_html = '<p class="no-data">Organism not found in BacDive public database.</p>'
    else:
        bd_html = '<p class="no-data">BacDive data not yet fetched for this organism.</p>'

    plastics = org_row.get("plastics", [])
    if isinstance(plastics, str):
        try:    plastics = ast.literal_eval(plastics)
        except: plastics = [p.strip() for p in plastics.split(",")]
    pills = "".join(f'<span class="pill">{p}</span>' for p in plastics)

    fy = int(org_row["first_year"]) if pd.notna(org_row.get("first_year")) else "N/A"
    ly = int(org_row["last_year"])  if pd.notna(org_row.get("last_year"))  else "N/A"

    return f"""<!DOCTYPE html><html lang="en">
<head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{name} | Organism Profile</title><style>{CSS}</style></head><body>
<div class="header"><h1>{name}</h1><div class="meta">
Tax ID: {org_row.get("tax_id","N/A")} | {org_row.get("n_plastics",0)} plastic type(s) |
{org_row.get("n_entries",0)} PlasticDB entries | Years: {fy} to {ly}
</div></div>
<div class="breadcrumb"><a href="index.html">All organisms</a> / {name}</div>
<div class="container">
<div class="section"><div class="section-header">PlasticDB Summary</div>
<div class="section-body"><div class="grid2"><div>
{_kv("Plastics degraded", pills)}
{_kv("Thermophilic flag", org_row.get("thermophilic","") or "Not recorded")}
{_kv("Has linked sequence", _badge(org_row.get("has_sequence")))}
{_kv("Has named enzyme",    _badge(org_row.get("has_enzyme")))}
{_kv("Has GenBank ID",      _badge(org_row.get("has_genbank")))}
</div><div>
{_kv("Isolation environments", org_row.get("isolation_envs","") or "Not recorded", pill=True)}
{_kv("Isolation locations",    org_row.get("isolation_locs","") or "Not recorded", pill=True)}
</div></div></div></div>
<div class="section"><div class="section-header">PlasticDB Entries ({org_row.get("n_entries",0)} total)</div>
<div class="section-body">{entries_html}</div></div>
<div class="section"><div class="section-header">NCBI SRA Sequencing Data</div>
<div class="section-body">{sra_html}</div></div>
<div class="section"><div class="section-header">BacDive Physiological Data</div>
<div class="section-body">{bd_html}</div></div>
</div></body></html>"""


def run_pages(full_df, orgs, sra_df, bd_df):
    sra_idx = sra_df.set_index("organism") if not sra_df.empty else pd.DataFrame()
    bd_idx  = bd_df.set_index("organism")  if not bd_df.empty  else pd.DataFrame()

    for i, (_, row) in enumerate(orgs.iterrows(), 1):
        name    = row["organism"]
        s       = _slug(name)
        sra_row = sra_idx.loc[name] if (not sra_idx.empty and name in sra_idx.index) else None
        bd_row  = bd_idx.loc[name]  if (not bd_idx.empty  and name in bd_idx.index)  else None
        html    = _make_page(row, full_df, sra_row, bd_row)
        (PAGES_DIR / f"{s}.html").write_text(html, encoding="utf-8")

    # Index
    import ast
    sra_map = dict(zip(sra_df["organism"], sra_df.get("sra_run_count", pd.Series()))) \
        if not sra_df.empty else {}
    bd_map  = dict(zip(bd_df["organism"],  bd_df.get("bacdive_found",  pd.Series()))) \
        if not bd_df.empty else {}

    rows_html = ""
    for _, row in orgs.sort_values("n_entries", ascending=False).iterrows():
        org = row["organism"]
        s   = _slug(org)
        try: rc = int(float(sra_map.get(org,"")))
        except: rc = "N/A"
        plastics = row.get("plastics",[])
        if isinstance(plastics, str):
            try:    plastics = ast.literal_eval(plastics)
            except: plastics = [p.strip() for p in plastics.split(",")]
        pills = "".join(f'<span class="pill">{p}</span>' for p in plastics[:5])
        if len(plastics) > 5: pills += f'<span class="pill">+{len(plastics)-5}</span>'
        bd_found = str(bd_map.get(org,"")).lower() == "yes"
        bd_badge = _badge(bd_found)
        fy = int(row["first_year"]) if pd.notna(row.get("first_year")) else "N/A"
        rows_html += (
            f'<tr><td><a href="{s}.html">{org}</a></td>'
            f'<td style="text-align:center">{row.get("n_entries",0)}</td>'
            f'<td>{pills}</td>'
            f'<td style="text-align:center">{rc}</td>'
            f'<td style="text-align:center">{bd_badge}</td>'
            f'<td style="text-align:center">{fy}</td></tr>'
        )

    index_html = f"""<!DOCTYPE html><html lang="en">
<head><meta charset="UTF-8"><title>Organism Profiles</title><style>{CSS}</style></head><body>
<div class="header"><h1>Organism Profiles</h1>
<div class="meta">{len(orgs)} organisms from PlasticDB, cross-referenced with NCBI SRA and BacDive.
All values fetched from live APIs. No values are estimated or interpolated.</div></div>
<div class="container"><div class="section"><div class="section-header">All Organisms</div>
<div class="section-body">
<input id="search-input" placeholder="Filter by organism name..." oninput="filterTable(this.value)">
<table id="org-table"><thead><tr>
<th>Organism</th><th style="text-align:center">PlasticDB entries</th>
<th>Plastics</th><th style="text-align:center">SRA runs</th>
<th style="text-align:center">BacDive</th><th style="text-align:center">First year</th>
</tr></thead><tbody>{rows_html}</tbody></table>
</div></div></div>
<script>function filterTable(q){{
  q=q.toLowerCase();
  document.querySelectorAll('#org-table tbody tr').forEach(r=>{{
    r.style.display=r.cells[0].textContent.toLowerCase().includes(q)?'':'none';
  }});
}}</script></body></html>"""

    (PAGES_DIR / "index.html").write_text(index_html, encoding="utf-8")
    total_pages = len(list(PAGES_DIR.glob("*.html")))
    print(f"Pages generated: {total_pages} (including index.html).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode",  choices=["sra","bd","pages","status"], required=True)
    parser.add_argument("--batch", type=int, default=200)
    args = parser.parse_args()

    full_df, orgs = load_plasticdb()
    print(f"Loaded {len(full_df):,} PlasticDB entries, {len(orgs):,} organisms.")

    if args.mode == "status":
        sra_cached = len(pd.read_csv(SRA_CSV)) if SRA_CSV.exists() else 0
        bd_cached  = len(pd.read_csv(BACDIVE_CSV)) if BACDIVE_CSV.exists() else 0
        pages      = len(list(PAGES_DIR.glob("*.html")))
        print(f"SRA cached:     {sra_cached}/{len(orgs)}")
        print(f"BacDive cached: {bd_cached}/{len(orgs)}")
        print(f"Pages built:    {pages}")

    elif args.mode == "sra":
        run_sra_batch(orgs, args.batch)

    elif args.mode == "bd":
        run_bacdive_batch(orgs, args.batch)

    elif args.mode == "pages":
        sra_df = pd.read_csv(SRA_CSV, dtype=str)     if SRA_CSV.exists()     else pd.DataFrame()
        bd_df  = pd.read_csv(BACDIVE_CSV, dtype=str) if BACDIVE_CSV.exists() else pd.DataFrame()
        run_pages(full_df, orgs, sra_df, bd_df)
