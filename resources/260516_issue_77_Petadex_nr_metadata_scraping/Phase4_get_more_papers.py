"""
Phase 4: Wiley TDM API Fetcher
Retrieves full-text PDFs from Wiley Online Library utilizing their Text and Data Mining API.
Implements specific user-agent headers and authentication flows for institutional access.
"""
import csv
import os
import re
import time
import urllib.error
import urllib.request

# Configuration & Paths
ROOT = os.path.dirname(os.path.abspath(__file__))
MISSING = os.path.join(ROOT, "data", "still_missing.csv")
OUT_DIR = os.path.join(ROOT, "data", "wiley_pdf")
LOG_CSV = os.path.join(ROOT, "data", "wiley_fetch_log.csv")

WILEY_TOKEN = "******************" #nononnononono
API_URL = "https://api.wiley.com/onlinelibrary/tdm/v1/articles/"
MAILTO = "aditya.ghosh@mail.utoronto.ca"
DELAY = 0.5  
RETRY_BACKOFF = [2, 5] 

def doi_to_filename(doi):
    """Sanitizes DOI strings for secure local filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe) + ".pdf"

def fetch_one(doi):
    """Executes a single authenticated API request with backoff logic."""
    clean_doi = doi.strip().replace("http://dx.doi.org/", "").replace("https://doi.org/", "")
    url = API_URL + clean_doi
    
    headers = {
        "Wiley-TDM-Client-OS-Version": "1.0",
        "Wiley-TDM-Client-Token": WILEY_TOKEN,
        "Accept": "application/pdf",
        "User-Agent": f"SciBERT-Proteomics/1.0 (mailto:{MAILTO})"
    }
    req = urllib.request.Request(url, headers=headers)

    for attempt in range(len(RETRY_BACKOFF) + 1):
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = resp.read()
                return resp.status, data, None
        except urllib.error.HTTPError as e:
            if e.code in (401, 403, 404):
                return e.code, None, f"HTTPError: {e.code} {e.reason}"
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return e.code, None, f"HTTPError: {e.code} {e.reason}"
        except Exception as e:
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return 0, None, f"System Error: {str(e)}"
    return 0, None, "exhausted retries"

def main():
    """Executes the batch retrieval pipeline for Wiley DOIs."""
    if WILEY_TOKEN == "YOUR_WILEY_TDM_TOKEN":
        return

    os.makedirs(OUT_DIR, exist_ok=True)

    dois = []
    try:
        with open(MISSING, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                doi = r.get("doi")
                if doi and (doi.startswith("10.1002/") or doi.startswith("10.1111/")):
                    dois.append(doi)
    except FileNotFoundError:
        return

    existing = set(os.listdir(OUT_DIR)) if os.path.isdir(OUT_DIR) else set()
    todo = [d for d in dois if doi_to_filename(d) not in existing]

    if not todo:
        return

    ok = not_entitled = not_found = phantom_pdfs = errors = 0
    log_rows = []

    for doi in todo:
        status, data, err = fetch_one(doi)
        outcome = ""

        if data and status == 200:
            if len(data) < 15000:
                phantom_pdfs += 1
                outcome = "phantom_rejected"
            else:
                fname = doi_to_filename(doi)
                with open(os.path.join(OUT_DIR, fname), "wb") as f:
                    f.write(data)
                ok += 1
                outcome = "ok"
        elif status in (401, 403):
            not_entitled += 1
            outcome = "not_entitled"
        elif status == 404:
            not_found += 1
            outcome = "not_found"
        else:
            errors += 1
            outcome = f"error_{status}"

        log_rows.append({"doi": doi, "outcome": outcome, "error": err or ""})
        time.sleep(DELAY)

    with open(LOG_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["doi", "outcome", "error"])
        w.writeheader()
        w.writerows(log_rows)

if __name__ == "__main__":
    main()
