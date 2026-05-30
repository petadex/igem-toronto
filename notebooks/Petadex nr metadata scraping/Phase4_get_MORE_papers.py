"""
Phase 4: CrossRef TDM Link Fetcher
Queries CrossRef for Text and Data Mining (TDM) links associated with target DOIs.
Downloads available full-text XML or PDF documents and filters out files that are too small.
"""
import csv
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

# Configuration & Paths
ROOT = os.path.dirname(os.path.abspath(__file__))
MISSING = os.path.join(ROOT, "data", "still_missing.csv")
OUT_DIR = os.path.join(ROOT, "data", "crossref_tdm")
LOG_CSV = os.path.join(ROOT, "data", "crossref_tdm_log.csv")

MAILTO = "aditya.ghosh@mail.utoronto.ca"
CROSSREF_DELAY = 0.5   
FETCH_DELAY = 0.5      
TIMEOUT = 30
RETRY_BACKOFF = [2, 5]

def doi_to_safe(doi):
    """Sanitizes DOI strings for local filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe)

def already_fetched_dirs():
    """Compiles a set of previously downloaded files to prevent redundant network requests."""
    dirs = [
        os.path.join(ROOT, "data", "elsevier_xml"),
        os.path.join(ROOT, "data", "wiley_pdf"),
    ]
    names = set()
    for d in dirs:
        if os.path.isdir(d):
            for f in os.listdir(d):
                stem = f.rsplit(".", 1)[0] if "." in f else f
                names.add(stem.lower())
    return names

def get_tdm_links(doi):
    """Retrieves TDM links from CrossRef metadata prioritizing XML formats."""
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi, safe='')}"
    headers = {
        "User-Agent": f"SciBERT-Proteomics/1.0 (mailto:{MAILTO})",
        "Accept": "application/json",
    }
    req = urllib.request.Request(url, headers=headers)

    for attempt in range(len(RETRY_BACKOFF) + 1):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                data = json.loads(resp.read())
                break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, None
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return None, None
        except Exception:
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return None, None

    links = data.get("message", {}).get("link", [])
    xml_link = pdf_link = other_link = None
    for link in links:
        if link.get("intended-application") != "text-mining":
            continue
        ct = link.get("content-type", "")
        u = link.get("URL", "")
        if not u:
            continue
        if "xml" in ct:
            xml_link = (u, ct)
        elif "pdf" in ct:
            pdf_link = (u, ct)
        elif not other_link:
            other_link = (u, ct)

    return xml_link or pdf_link or other_link or (None, None)

def fetch_fulltext(url, accept):
    """Downloads the document payload from the resolved TDM URL."""
    headers = {
        "User-Agent": f"SciBERT-Proteomics/1.0 (mailto:{MAILTO})",
        "Accept": accept,
    }
    req = urllib.request.Request(url, headers=headers)

    for attempt in range(len(RETRY_BACKOFF) + 1):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                data = resp.read()
                return resp.status, data, None
        except urllib.error.HTTPError as e:
            if e.code in (401, 403, 404):
                return e.code, None, f"HTTPError: {e.code}"
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return e.code, None, f"HTTPError: {e.code}"
        except Exception as e:
            if attempt < len(RETRY_BACKOFF):
                time.sleep(RETRY_BACKOFF[attempt])
                continue
            return 0, None, f"System Error: {str(e)}"
    return 0, None, "exhausted retries"

def detect_format(data):
    """Verifies the integrity and format of the downloaded byte sequence."""
    if data[:4] == b"%PDF":
        return "pdf"
    head = data[:200].lstrip()
    if head[:5] == b"<?xml" or head[:1] == b"<":
        return "xml"
    return None

def main():
    """Executes the CrossRef TDM batch retrieval process."""
    os.makedirs(OUT_DIR, exist_ok=True)

    all_dois = []
    try:
        with open(MISSING, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("doi"):
                    all_dois.append(r["doi"])
    except FileNotFoundError:
        return

    already_secured = already_fetched_dirs()
    own_files = set()
    if os.path.isdir(OUT_DIR):
        for f in os.listdir(OUT_DIR):
            stem = f.rsplit(".", 1)[0] if "." in f else f
            own_files.add(stem.lower())

    todo = [doi for doi in all_dois if doi_to_safe(doi) not in own_files and doi_to_safe(doi) not in already_secured]
    if not todo:
        return

    ok = no_link = not_entitled = phantoms = bad_content = errors = 0
    log_rows = []

    for doi in todo:
        tdm_url, content_type = get_tdm_links(doi)
        time.sleep(CROSSREF_DELAY)

        outcome = ""
        if not tdm_url:
            no_link += 1
            outcome = "no_tdm_link"
        else:
            accept = content_type if content_type else "*/*"
            status, data, err = fetch_fulltext(tdm_url, accept)
            time.sleep(FETCH_DELAY)

            if data and status == 200:
                # Reject metadata-only files masquerading as full text (< 15KB)
                if len(data) < 15000:
                    phantoms += 1
                    outcome = "phantom_rejected"
                else:
                    fmt = detect_format(data)
                    if fmt:
                        ext = ".xml" if fmt == "xml" else ".pdf"
                        fname = doi_to_safe(doi) + ext
                        with open(os.path.join(OUT_DIR, fname), "wb") as f:
                            f.write(data)
                        ok += 1
                        outcome = "ok"
                    else:
                        bad_content += 1
                        outcome = "unrecognized_format"
            elif status in (401, 403):
                not_entitled += 1
                outcome = "not_entitled"
            else:
                errors += 1
                outcome = f"error_{status}"

        log_rows.append({"doi": doi, "outcome": outcome})

    with open(LOG_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["doi", "outcome"])
        w.writeheader()
        w.writerows(log_rows)

if __name__ == "__main__":
    main()
