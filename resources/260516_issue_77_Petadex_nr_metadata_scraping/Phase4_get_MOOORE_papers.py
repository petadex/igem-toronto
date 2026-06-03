"""
Phase 4: Europe PMC Dual-Engine Fetcher
Queries the Europe PMC REST API for Open Access PDF endpoints.
Utilizes a stealth HTTP/2 client with an automatic fallback to standard HTTP/1.1 requests to ensure connection stability.
"""
import os
import csv
import re
import time
from curl_cffi import requests as stealth_requests
import requests as std_requests

# Configuration & Paths
MASTER_DIR = r"C:\Users\somna\source\repos\why do i have so many directories\master audit"
INPUT_CSV = r"C:\Users\somna\source\repos\scraping logan attempt 16\scraping logan attempt 16\data\ABSOLUTE_FINAL_MISSING_DOIS.csv"
CLOUDFLARE_DIR = os.path.join(MASTER_DIR, "CLOUDFLARE_PDFS")
EPMC_DIR = os.path.join(MASTER_DIR, "EUROPE_PMC_PDFS")

def doi_to_safe_filename(doi):
    """Sanitizes DOI strings for local filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe)

def main():
    """Executes the Europe PMC retrieval pipeline."""
    os.makedirs(EPMC_DIR, exist_ok=True)
    
    secured_safe_dois = set()
    if os.path.exists(CLOUDFLARE_DIR):
        for filename in os.listdir(CLOUDFLARE_DIR):
            if filename.lower().endswith(".pdf"):
                stem = filename.rsplit(".", 1)[0].lower()
                secured_safe_dois.add(stem)
                
    targets = []
    if os.path.exists(INPUT_CSV):
        with open(INPUT_CSV, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                doi = row.get("doi")
                if doi:
                    safe_doi = doi_to_safe_filename(doi)
                    if safe_doi not in secured_safe_dois:
                        targets.append(doi.strip())
                        
    successful_downloads = 0
    session = stealth_requests.Session(impersonate="chrome110")
    
    for i, doi in enumerate(targets, 1):
        safe_name = doi_to_safe_filename(doi)
        pdf_path = os.path.join(EPMC_DIR, f"{safe_name}.pdf")
        
        if os.path.exists(pdf_path):
            continue

        api_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        params = {
            "query": f'DOI:"{doi}"',
            "resultType": "core",
            "format": "json"
        }
        
        try:
            response = session.get(api_url, params=params, timeout=15)
            data = response.json()
            
            hit_count = data.get("hitCount", 0)
            if hit_count > 0:
                result = data.get("resultList", {}).get("result", [])[0]
                has_pdf = result.get("hasPDF", "N")
                pmcid = result.get("pmcid")
                
                if has_pdf == "Y" and pmcid:
                    download_url = f"https://europepmc.org/articles/{pmcid}?pdf=render"
                    headers = {
                        "Referer": f"https://europepmc.org/articles/{pmcid}",
                        "Accept": "application/pdf,application/xhtml+xml,text/html,*/*"
                    }
                    
                    try:
                        # Primary Engine: Stealth HTTP/2
                        pdf_response = session.get(download_url, headers=headers, timeout=30, allow_redirects=True)
                        
                        if pdf_response.status_code == 200 and b'%PDF' in pdf_response.content[:10]:
                            with open(pdf_path, "wb") as pdf_file:
                                pdf_file.write(pdf_response.content)
                            successful_downloads += 1
                            
                    except Exception as e:
                        if "STREAM_CLOSED" in str(e) or "HTTP/2" in str(e):
                            try:
                                # Fallback Engine: Standard HTTP/1.1
                                fallback_headers = {
                                    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Chrome/120.0.0.0 Safari/537.36",
                                    **headers
                                }
                                fb_response = std_requests.get(download_url, headers=fallback_headers, timeout=30, allow_redirects=True)
                                
                                if fb_response.status_code == 200 and b'%PDF' in fb_response.content[:10]:
                                    with open(pdf_path, "wb") as pdf_file:
                                        pdf_file.write(fb_response.content)
                                    successful_downloads += 1
                            except Exception:
                                pass
        except Exception:
            pass
            
        time.sleep(1.5)

if __name__ == "__main__":
    main()
