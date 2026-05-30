"""
Phase 4: Elsevier TDM API Fetcher
Retrieves full-text XML documents from ScienceDirect utilizing the Elsevier Article Retrieval API.
Includes a byte-size filter that rejects full text downloads that clearly don't have the full text.
"""
import csv
import os
import re
import tempfile
import time
import urllib.error
import urllib.request

# Configuration & Paths
ROOT = os.path.dirname(os.path.abspath(__file__))
MISSING = os.path.join(ROOT, "data", "still_missing.csv")
OUT_DIR = os.path.join(ROOT, "data", "elsevier_xml")
LOG_CSV = os.path.join(ROOT, "data", "elsevier_fetch_log.csv")

API_KEY = "*****************" #nononono
API_URL = "https://api.elsevier.com/content/article/doi/"
MAILTO = "aditya.ghosh@mail.utoronto.ca"
DELAY = 0.12          
RETRY_BACKOFF = [2, 4, 8]

def doi_to_filename(doi):
    """Sanitizes DOI strings for secure local filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe) + ".xml"

def fetch_one(doi):
    """Executes a single API request with exponential backoff for rate limiting."""
    url = API_URL + doi
    headers = {
        "X-ELS-APIKey": API_KEY,
        "Accept": "text/xml",
        "User-Agent": f"SciBERT-Proteomics/1.0 (mailto:{MAILTO})",
    }
    req = urllib.request.Request(url, headers=headers)

    for attempt in range(len(RETRY_BACKOFF) + 1):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                return resp.status, data, None
        except urllib.error.HTTPError as e:
            code = e.code
            if code in (401, 403, 404):
                return code, None, str(e)
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return code, None, str(e)
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return 0, None, str(e)
    return 0, None, "exhausted retries"

def main():
    """Executes the batch retrieval pipeline for Elsevier DOIs."""
    if API_KEY == "YOUR_ELSEVIER_API_KEY":
        return

    os.makedirs(OUT_DIR, exist_ok=True)

    dois = []
    try:
        with open(MISSING, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                doi = r.get("doi")
                if doi and doi.startswith("10.1016/"):
                    dois.append(doi)
    except FileNotFoundError:
        return

    existing = set(os.listdir(OUT_DIR)) if os.path.isdir(OUT_DIR) else set()
    todo = [d for d in dois if doi_to_filename(d) not in existing]

    if not todo:
        return

    ok = not_entitled = not_found = phantom_abstracts = errors = 0
    log_rows = []

    for doi in todo:
        status, data, err = fetch_one(doi)

        if data and status == 200:
            if len(data) < 15000:
                phantom_abstracts += 1
                outcome = "metadata_only_rejected"
            else:
                fname = doi_to_filename(doi)
                fd, tmp = tempfile.mkstemp(dir=OUT_DIR, suffix=".tmp")
                try:
                    os.write(fd, data)
                    os.close(fd)
                    os.replace(tmp, os.path.join(OUT_DIR, fname))
                    ok += 1
                    outcome = "ok"
                except Exception:
                    os.close(fd) if not os.get_inheritable(fd) else None
                    if os.path.exists(tmp):
                        os.remove(tmp)
                    raise
        elif status in (401, 403):
            not_entitled += 1
            outcome = "not_entitled"
        elif status == 404:
            not_found += 1
            outcome = "not_found"
        else:
            errors += 1
            outcome = f"error_{status}"

        log_rows.append({
            "doi": doi, "http_status": status,
            "outcome": outcome, "error": err or "",
        })
        time.sleep(DELAY)

    with open(LOG_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["doi", "http_status", "outcome", "error"])
        w.writeheader()
        w.writerows(log_rows)

if __name__ == "__main__":
    main()
