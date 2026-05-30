"""
Phase 3: DOI Resolution Engine
Merges PMIDs from earlier pipeline phases, identifies unretrieved documents, and resolves PMIDs to DOIs using the Europe PMC API.
Outputs a master ledger of missing documents for downstream retrieval scripts.
(this feels out of place because it is, this was done at a later step when I had scraped some but not others and was unsure which papers I had)
"""
import os
import csv
import requests
import pandas as pd
from tqdm import tqdm

# Configuration & Paths
EMAIL = "aditya.ghosh@mail.utoronto.ca" 
PHASE1_PMIDS_FILE = "s3_pmid_mapping.csv"
PHASE2_PMIDS_FILE = "final_perfect_pmid_mapping.csv"
OUTPUT_DIR = "fulltext_articles"

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, "data")
MISSING_CSV_PATH = os.path.join(DATA_DIR, "still_missing.csv")

def main():
    """Executes the PMID merging and DOI resolution pipeline."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)

    phase1_pmids = set()
    try:
        df1 = pd.read_csv(PHASE1_PMIDS_FILE)
        success_df1 = df1[(df1['PMID'].notna()) & (df1['PMID'] != 'None')]
        all_pmids1 = success_df1['PMID'].astype(str).str.split(';').explode().str.strip()
        phase1_pmids = set(all_pmids1.unique())
    except Exception:
        pass

    phase2_pmids = set()
    try:
        df2 = pd.read_csv(PHASE2_PMIDS_FILE)
        df2['PMID_clean'] = df2['PMID'].astype(str).str.strip().str.replace('.0', '', regex=False)
        valid_df2 = df2[(df2['PMID_clean'].notna()) & (df2['PMID_clean'] != 'nan') & (df2['PMID_clean'] != 'None')]
        phase2_pmids = set(valid_df2['PMID_clean'].tolist())
    except Exception:
        pass

    master_pmids = phase1_pmids.union(phase2_pmids)
    completed = {f.replace("_fulltext", "").split(".")[0] for f in os.listdir(OUTPUT_DIR)}
    pmids_to_process = list(master_pmids - completed)

    if not pmids_to_process:
        return

    resolved_data = []
    batch_size = 40 
    batches = [pmids_to_process[i:i + batch_size] for i in range(0, len(pmids_to_process), batch_size)]

    generic_headers = {
        "User-Agent": f"SciBERT-Proteomics/1.0 (mailto:{EMAIL})"
    }

    for batch in tqdm(batches, desc="Phase 3: Resolving DOIs via EPMC"):
        query_string = " OR ".join([f"ext_id:{pmid}" for pmid in batch])
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        
        params = {
            "query": query_string,
            "format": "json",
            "resultType": "lite",
            "pageSize": 1000 
        }
        
        try:
            resp = requests.get(url, params=params, headers=generic_headers, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                results = data.get("resultList", {}).get("result", [])
                for res in results:
                    doi = res.get("doi")
                    pmid = res.get("pmid")
                    if doi and pmid:
                        resolved_data.append({"pmid": pmid, "doi": doi})
        except Exception:
            pass 
            
    with open(MISSING_CSV_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["pmid", "doi"])
        writer.writeheader()
        writer.writerows(resolved_data)

if __name__ == "__main__":
    main()
