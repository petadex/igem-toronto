"""
Phase 3: Set Consolidation and API Retrieval
Merges PubMed IDs from Phase 1 and Phase 2 into a deduplicated master set.
Queries Unpaywall, Elsevier, and Springer APIs to fetch available full-text documents. (This paer did not work)
"""
import os
import time
import pandas as pd
import requests
from tqdm import tqdm

# Configuration
ELSEVIER_API_KEY = "*******************" #nonono
SPRINGER_API_KEY = "*******************" #nonono
EMAIL = "aditya.ghosh@mail.utoronto.ca"
OUTPUT_DIR = "fulltext_articles"

def load_and_merge_pmids():
    """Loads PMIDs from previous pipeline phases and returns a deduplicated master set."""
    phase1_pmids = set()
    try:
        with open("unique_pmids_to_download.txt", "r") as f:
            phase1_pmids = {line.strip() for line in f if line.strip()}
    except FileNotFoundError:
        pass

    phase2_pmids = set()
    try:
        df = pd.read_csv("final_perfect_pmid_mapping.csv")
        df['PMID_clean'] = df['PMID'].astype(str).str.strip().str.replace('.0', '', regex=False)
        valid_df = df[(df['PMID_clean'].notna()) & (df['PMID_clean'] != 'nan') & (df['PMID_clean'] != 'None')]
        phase2_pmids = set(valid_df['PMID_clean'].tolist())
    except FileNotFoundError:
        pass

    return phase1_pmids.union(phase2_pmids)

def main():
    """Executes the master consolidation and document retrieval process."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    master_pmids = load_and_merge_pmids()
    
    completed = {f.replace("_fulltext", "").split(".")[0] for f in os.listdir(OUTPUT_DIR)}
    pmids_to_process = [p for p in master_pmids if p not in completed]

    if not pmids_to_process:
        return

    generic_headers = {
        "User-Agent": f"Mozilla/5.0 (Windows NT 10.0; Win64; x64) ProteomicsPipeline/1.0 (mailto:{EMAIL})"
    }

    for pmid in tqdm(pmids_to_process, desc="Phase 3: Fetching Documents"):
        doi = None
        downloaded = False
        
        # Layer 1: Unpaywall Open Access Resolution
        try:
            up_url = f"https://api.unpaywall.org/v2/{pmid}?email={EMAIL}"
            up_resp = requests.get(up_url, timeout=10)
            if up_resp.status_code == 200:
                up_data = up_resp.json()
                doi = up_data.get("doi")
                if up_data.get("is_oa") and up_data.get("best_oa_location"):
                    oa_url = up_data["best_oa_location"].get("url_for_pdf") or up_data["best_oa_location"].get("url")
                    if oa_url:
                        oa_resp = requests.get(oa_url, headers=generic_headers, timeout=15)
                        if oa_resp.status_code == 200:
                            ext = "pdf" if "pdf" in oa_url.lower() else "html"
                            with open(os.path.join(OUTPUT_DIR, f"{pmid}.{ext}"), "wb") as f:
                                f.write(oa_resp.content)
                            downloaded = True
        except Exception:
            pass

        if downloaded:
            time.sleep(0.1)
            continue

        # Europe PMC Fallback for missing DOIs
        if not doi:
            try:
                epmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ext_id:{pmid} src:med&format=json&resultType=lite"
                epmc_data = requests.get(epmc_url, timeout=10).json()
                if epmc_data.get("hitCount", 0) > 0:
                    doi = epmc_data["resultList"]["result"][0].get("doi")
            except Exception:
                pass

        if not doi:
            continue

        # Layer 2: Elsevier XML Full-Text
        if "10.1016" in doi or "10.2139" in doi:
            try:
                el_url = f"https://api.elsevier.com/content/article/doi/{doi}"
                el_headers = {"X-ELS-APIKey": ELSEVIER_API_KEY, "Accept": "text/xml"}
                el_resp = requests.get(el_url, headers=el_headers, timeout=15)
                
                if el_resp.status_code == 200:
                    content = el_resp.content
                    if len(content) > 15000:
                        with open(os.path.join(OUTPUT_DIR, f"{pmid}.xml"), "wb") as f:
                            f.write(content)
                        downloaded = True
            except Exception:
                pass

        # Layer 3: Springer Nature Open Access JATS XML
        elif ("10.1007" in doi or "10.1038" in doi) and not downloaded:
            try:
                springer_url = f"https://api.springernature.com/openaccess/jats?q=doi:{doi}&api_key={SPRINGER_API_KEY}"
                sp_resp = requests.get(springer_url, timeout=15)
                
                if sp_resp.status_code == 200:
                    content = sp_resp.content
                    if len(content) > 15000:
                        with open(os.path.join(OUTPUT_DIR, f"{pmid}.xml"), "wb") as f:
                            f.write(content)
                        downloaded = True
            except Exception:
                pass

        # Layer 4: Universal Publisher HTML Link-Out
        if not downloaded:
            try:
                doi_url = f"https://doi.org/{doi}"
                pub_resp = requests.get(doi_url, headers=generic_headers, timeout=15)
                if pub_resp.status_code == 200:
                    if b"login" not in pub_resp.url.lower():
                        with open(os.path.join(OUTPUT_DIR, f"{pmid}.html"), "wb") as f:
                            f.write(pub_resp.content)
            except Exception:
                pass

        time.sleep(0.2)

if __name__ == "__main__":
    main()
