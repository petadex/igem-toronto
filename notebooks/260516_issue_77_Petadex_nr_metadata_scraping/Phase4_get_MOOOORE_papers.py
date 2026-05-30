"""
Phase 4: Unpaywall Resolution & Splash Page Piercer
Queries the Unpaywall API for Open Access endpoints.
Implements a DOM parser to extract underlying PDF assets when endpoints route to HTML landing pages.
"""
import os
import csv
import time
import re
import requests as std_requests
from curl_cffi import requests as stealth_requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin

# Configuration & Paths
MASTER_DIR = r"C:\Users\somna\source\repos\why do i have so many directories\master audit"
INPUT_CSV = r"C:\Users\somna\source\repos\scraping logan attempt 16\scraping logan attempt 16\data\ABSOLUTE_FINAL_MISSING_DOIS.csv"
CLOUDFLARE_DIR = os.path.join(MASTER_DIR, "CLOUDFLARE_PDFS")
UNPAYWALL_DIR = os.path.join(MASTER_DIR, "UNPAYWALL_PDFS")
EMAIL = "aditya.ghosh@mail.utoronto.ca"

def doi_to_safe_filename(doi):
    """Sanitizes DOI strings for local filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe)

def extract_pdf_link(html_content, base_url):
    """Parses HTML DOM to locate standard academic PDF metadata tags or actionable download links."""
    soup = BeautifulSoup(html_content, 'html.parser')
    
    meta_tag = soup.find('meta', attrs={'name': 'citation_pdf_url'})
    if meta_tag and meta_tag.get('content'):
        return meta_tag['content']
    
    for a_tag in soup.find_all('a', href=True):
        if 'pdf' in a_tag['href'].lower() and ('download' in a_tag.text.lower() or 'pdf' in a_tag.text.lower()):
            return urljoin(base_url, a_tag['href'])
            
    return None

def main():
    """Executes the Unpaywall retrieval and HTML parsing pipeline."""
    os.makedirs(UNPAYWALL_DIR, exist_ok=True)
    
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
                    doi = doi.strip()
                    if doi.lower().startswith("10.2210/pdb"):
                        continue
                    
                    safe_doi = doi_to_safe_filename(doi)
                    if safe_doi not in secured_safe_dois:
                        targets.append(doi)
                        
    successful_downloads = 0
    session = stealth_requests.Session(impersonate="chrome110")
    
    for doi in targets:
        safe_name = doi_to_safe_filename(doi)
        pdf_path = os.path.join(UNPAYWALL_DIR, f"{safe_name}.pdf")
        
        if os.path.exists(pdf_path):
            continue

        api_url = f"https://api.unpaywall.org/v2/{doi}?email={EMAIL}"
        
        try:
            response = std_requests.get(api_url, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                is_oa = data.get("is_oa", False)
                best_oa_location = data.get("best_oa_location")
                
                if is_oa and best_oa_location:
                    pdf_url = best_oa_location.get("url_for_pdf") or best_oa_location.get("url_for_landing_page")
                    
                    if pdf_url:
                        try:
                            pdf_response = session.get(pdf_url, timeout=30, allow_redirects=True)
                            
                            if pdf_response.status_code == 200:
                                if b'%PDF' in pdf_response.content[:10]:
                                    with open(pdf_path, "wb") as pdf_file:
                                        pdf_file.write(pdf_response.content)
                                    successful_downloads += 1
                                else:
                                    real_pdf_link = extract_pdf_link(pdf_response.text, pdf_response.url)
                                    if real_pdf_link:
                                        final_response = session.get(real_pdf_link, timeout=30)
                                        if final_response.status_code == 200 and b'%PDF' in final_response.content[:10]:
                                            with open(pdf_path, "wb") as pdf_file:
                                                pdf_file.write(final_response.content)
                                            successful_downloads += 1
                        except Exception:
                            pass
        except Exception:
            pass
            
        time.sleep(1)

if __name__ == "__main__":
    main()
