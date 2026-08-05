"""
Fetch BacDive physiological data for organisms from SAMN4.csv / unique_taxa.csv.
Robust column-based parser version.
"""

import os
import re
import time
import argparse
from pathlib import Path
import pandas as pd
import requests
from bs4 import BeautifulSoup

# Define paths
ROOT = Path(__file__).parent
DATA_DIR = ROOT.parent / "data"
UNIQUE_TAXA_CSV = DATA_DIR / "unique_taxa.csv"
BACDIVE_CSV = DATA_DIR / "bacdive_data.csv"

DATA_DIR.mkdir(exist_ok=True)

# Rate-limited HTTP session to be polite to the BacDive server
session = requests.Session()
session.headers["User-Agent"] = "igem-toronto-research/1.0 (academic, metadata-enrichment)"

def bacdive_search_ids(organism: str, sleep: float = 0.5) -> list:
    """Return BacDive strain IDs for an organism via public search page."""
    time.sleep(sleep)
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
    except Exception as e:
        print(f"Error searching for {organism}: {e}")
        return []

def extract_table_rows(soup: BeautifulSoup, section_title: str) -> list:
    """
    Finds the table under the section header and returns a list of dictionaries,
    where each dictionary maps header column names to their row values.
    """
    section_lower = section_title.lower()
    for h in soup.find_all(["h3", "h4"]):
        header_text = h.get_text().strip().lower()
        
        # Determine exact match condition to avoid partial match bugs
        if section_lower == "isolation":
            match = header_text == "isolation"
        elif section_lower == "ph":
            match = header_text == "culture ph"
        elif section_lower == "culture temp":
            match = header_text == "culture temp"
        elif section_lower == "oxygen":
            match = header_text == "oxygen tolerance"
        elif section_lower == "morphol":
            match = header_text == "cell morphology"
        else:
            match = section_lower in header_text
            
        if match:
            tbl = h.find_next("table")
            if tbl:
                headers = [th.get_text(strip=True) for th in tbl.find_all("th")]
                rows = []
                for tr in tbl.find_all("tr"):
                    tds = [td.get_text(strip=True) for td in tr.find_all("td")]
                    if tds and len(tds) <= len(headers):
                        tds += [""] * (len(headers) - len(tds))
                        row_dict = dict(zip(headers, tds))
                        rows.append(row_dict)
                return rows
    return []

def bacdive_strain(strain_id: str, sleep: float = 0.5) -> dict:
    """Scrape a BacDive public strain page and return structured fields."""
    time.sleep(sleep)
    url = f"https://bacdive.dsmz.de/strain/{strain_id}"
    try:
        r = session.get(url, timeout=20)
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "html.parser")
    except Exception as e:
        print(f"Error fetching strain {strain_id}: {e}")
        return {}

    data = {
        "bacdive_strain_id": strain_id,
        "bacdive_url": url,
        "bacdive_name": "",
        "bacdive_temp_c": "",
        "bacdive_ph": "",
        "bacdive_oxygen": "",
        "bacdive_gram": "",
        "bacdive_morphology": "",
        "bacdive_motility": "",
        "bacdive_isolation": ""
    }

    # Taxonomy / name
    name_tag = soup.find("h1")
    if name_tag:
        data["bacdive_name"] = name_tag.get_text(strip=True)

    # Temperature
    temp_rows = extract_table_rows(soup, "Culture temp")
    temps = []
    for r_dict in temp_rows:
        growth_val = next((v for k, v in r_dict.items() if "growth" in k.lower()), "")
        temp_val = next((v for k, v in r_dict.items() if "temp" in k.lower()), "")
        if "negative" not in growth_val.lower() and temp_val:
            temps.append(temp_val)
    if temps:
        data["bacdive_temp_c"] = "; ".join(dict.fromkeys(temps))

    # pH
    ph_rows = extract_table_rows(soup, "ph")
    phs = []
    for r_dict in ph_rows:
        ability_val = next((v for k, v in r_dict.items() if "ability" in k.lower() or "growth" in k.lower()), "")
        ph_val = next((v for k, v in r_dict.items() if k.lower() == "ph"), "")
        ph_range_val = next((v for k, v in r_dict.items() if "ph range" in k.lower()), "")
        val = ph_val or ph_range_val
        if "negative" not in ability_val.lower() and val and val.lower() not in ("alkaliphile", "acidophile"):
            phs.append(val)
    if phs:
        data["bacdive_ph"] = "; ".join(dict.fromkeys(phs))

    # Oxygen tolerance
    oxy_rows = extract_table_rows(soup, "Oxygen")
    oxys = []
    for r_dict in oxy_rows:
        oxy_val = next((v for k, v in r_dict.items() if "oxygen" in k.lower()), "")
        if oxy_val:
            oxys.append(oxy_val)
    if oxys:
        data["bacdive_oxygen"] = "; ".join(dict.fromkeys(oxys))

    # Gram stain / morphology
    morph_rows = extract_table_rows(soup, "morphol")
    gram_stains = []
    shapes = []
    motilities = []
    for r_dict in morph_rows:
        gram_val = next((v for k, v in r_dict.items() if "gram" in k.lower()), "")
        shape_val = next((v for k, v in r_dict.items() if "shape" in k.lower() or "morphology" in k.lower()), "")
        mot_val = next((v for k, v in r_dict.items() if "motil" in k.lower()), "")
        
        if gram_val: gram_stains.append(gram_val)
        if shape_val: shapes.append(shape_val)
        if mot_val: motilities.append(mot_val)
        
    if gram_stains:
        data["bacdive_gram"] = gram_stains[0]
    if shapes:
        data["bacdive_morphology"] = "; ".join(dict.fromkeys(shapes))
    if motilities:
        data["bacdive_motility"] = motilities[0]

    # Isolation
    isol_rows = extract_table_rows(soup, "Isolation")
    isolations = []
    for r_dict in isol_rows:
        parts = []
        source = next((v for k, v in r_dict.items() if "sample" in k.lower() or "isolate" in k.lower() or "source" in k.lower()), "")
        host = next((v for k, v in r_dict.items() if "host" in k.lower()), "")
        country = next((v for k, v in r_dict.items() if "country" in k.lower()), "")
        
        if source: parts.append(source)
        if host: parts.append(f"Host: {host}")
        if country: parts.append(country)
        
        if parts:
            isolations.append(", ".join(parts))
    if isolations:
        data["bacdive_isolation"] = "; ".join(dict.fromkeys(isolations[:6]))

    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int, default=50, help="Number of organisms to fetch in this run")
    parser.add_argument("--overwrite-cache", action="store_true", help="Delete the cache file and start fresh")
    args = parser.parse_args()

    if not UNIQUE_TAXA_CSV.exists():
        print(f"Error: {UNIQUE_TAXA_CSV} does not exist. Run process_athena_results.py first.")
        return

    # Load unique taxa
    taxa_df = pd.read_csv(UNIQUE_TAXA_CSV)
    taxa_df = taxa_df.dropna(subset=["taxon_name"])
    
    # Load cache if it exists and we're not overwriting
    if BACDIVE_CSV.exists() and not args.overwrite_cache:
        cached = pd.read_csv(BACDIVE_CSV, dtype=str)
        done = set(cached["organism"].tolist())
        print(f"Loaded cache: {len(cached)} organisms already processed.")
    else:
        cached = pd.DataFrame()
        done = set()
        print("Starting fresh (no cache used).")

    # Get list of organisms to process
    todo = taxa_df[~taxa_df["taxon_name"].isin(done)].copy()
    print(f"Total organisms in unique_taxa.csv: {len(taxa_df)}")
    print(f"Remaining to process: {len(todo)}")

    if todo.empty:
        print("All organisms already processed.")
        return

    todo_batch = todo.head(args.batch)
    print(f"Processing next batch of {len(todo_batch)} organisms...")

    rows = []
    for i, (_, row) in enumerate(todo_batch.iterrows(), 1):
        org = row["taxon_name"]
        taxon_id = row["taxon_id"]
        print(f"[{i}/{len(todo_batch)}] Searching '{org}' (Taxon ID: {taxon_id})...")
        ids = bacdive_search_ids(org)
        if not ids:
            print(f"  -> Not found in BacDive")
            rows.append({
                "organism": org,
                "taxon_id": taxon_id,
                "bacdive_found": "No",
                "bacdive_strain_id": "",
                "bacdive_url": "",
                "bacdive_name": "",
                "bacdive_temp_c": "",
                "bacdive_ph": "",
                "bacdive_oxygen": "",
                "bacdive_gram": "",
                "bacdive_morphology": "",
                "bacdive_motility": "",
                "bacdive_isolation": ""
            })
        else:
            print(f"  -> Found strain ID(s): {ids}. Scraping strain {ids[0]}...")
            d = bacdive_strain(ids[0])
            d["organism"] = org
            d["taxon_id"] = taxon_id
            d["bacdive_found"] = "Yes"
            rows.append(d)

    # Save progress
    new_df = pd.DataFrame(rows)
    combined = pd.concat([cached, new_df], ignore_index=True) if not cached.empty else new_df
    combined.to_csv(BACDIVE_CSV, index=False)
    
    found_count = (combined["bacdive_found"] == "Yes").sum()
    not_found_count = (combined["bacdive_found"] == "No").sum()
    print(f"\nSaved progress. Total cached: {len(combined)} (Yes: {found_count}, No: {not_found_count})")
    print(f"Output file: {BACDIVE_CSV}")

if __name__ == "__main__":
    main()
