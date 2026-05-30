"""
Phase 2: Structural Journal Mapping
Extracts unique citation strings from metadata and queries the Europe PMC API to resolve PubMed IDs.
Implements a two-pass system (Exact Match and Rescue Pass) to maximize resolution success rates.
(Theres a journal column in the dataset so might as well check it for new papers)
"""
import os
import time
import pandas as pd
import requests
from tqdm import tqdm

# Configuration
INPUT_CSV = "petadex_nr_metadata.csv"
INTERMEDIATE_CSV = "citation_to_pmid_mapping.csv"
FINAL_CSV = "final_perfect_pmid_mapping.csv"
BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

def initial_mapping_pass():
    """Extracts unique citations and performs the initial exact-match API search."""
    df = pd.read_csv(INPUT_CSV, usecols=['journal'])
    published_df = df.dropna(subset=['journal'])
    published_df = published_df[~published_df['journal'].str.contains('Unpublished', case=False)]
    unique_citations = published_df['journal'].unique()

    completed_citations = set()
    if os.path.exists(INTERMEDIATE_CSV):
        try:
            existing_df = pd.read_csv(INTERMEDIATE_CSV)
            completed_citations = set(existing_df['Citation'].astype(str))
        except Exception:
            pass
    else:
        with open(INTERMEDIATE_CSV, "w", encoding="utf-8") as f:
            f.write("Citation,PMID\n")

    citations_to_process = [c for c in unique_citations if c not in completed_citations]
    
    for citation in tqdm(citations_to_process, desc="Phase 2A: Exact Match Processing"):
        try:
            params = {
                "query": f'EXT_ID:"{citation}" OR "{citation}"',
                "format": "json",
                "resultType": "lite"
            }
            response = requests.get(BASE_URL, params=params, timeout=10)
            data = response.json()
            
            pmid = "None"
            if data.get("hitCount", 0) > 0:
                results = data.get("resultList", {}).get("result", [])
                for result in results:
                    if "pmid" in result:
                        pmid = result["pmid"]
                        break
            
            with open(INTERMEDIATE_CSV, "a", encoding="utf-8") as f:
                f.write(f'"{citation}",{pmid}\n')
                
        except Exception:
            time.sleep(2)
            continue
            
        time.sleep(0.2)

def rescue_mapping_pass():
    """Processes failed matches using a relaxed search query to maximize resolution."""
    df = pd.read_csv(INTERMEDIATE_CSV)
    df['PMID'] = df['PMID'].astype(str).str.replace('.0', '', regex=False)
    
    success_df = df[(df['PMID'].notna()) & (df['PMID'] != 'nan') & (df['PMID'] != 'None')]
    misses_df = df[(df['PMID'].isna()) | (df['PMID'] == 'nan') | (df['PMID'] == 'None')]

    success_df.to_csv(FINAL_CSV, index=False)

    for idx, row in tqdm(misses_df.iterrows(), total=len(misses_df), desc="Phase 2B: Rescue Pass Processing"):
        citation = str(row['Citation'])
        
        if "Submitted" in citation:
            with open(FINAL_CSV, "a", encoding="utf-8") as f:
                f.write(f'"{citation}",None\n')
            continue
            
        try:
            params = {
                "query": citation,
                "format": "json",
                "resultType": "lite"
            }
            response = requests.get(BASE_URL, params=params, timeout=10)
            data = response.json()
            
            pmid = "None"
            if data.get("hitCount", 0) > 0:
                results = data.get("resultList", {}).get("result", [])
                for result in results:
                    if "pmid" in result:
                        pmid = result["pmid"]
                        break
                        
            with open(FINAL_CSV, "a", encoding="utf-8") as f:
                f.write(f'"{citation}",{pmid}\n')
                
        except Exception:
            with open(FINAL_CSV, "a", encoding="utf-8") as f:
                f.write(f'"{citation}",None\n')
            time.sleep(2)
            
        time.sleep(0.15)

def main():
    """Executes the complete two-pass structural journal mapping pipeline."""
    initial_mapping_pass()
    rescue_mapping_pass()

if __name__ == "__main__":
    main()
