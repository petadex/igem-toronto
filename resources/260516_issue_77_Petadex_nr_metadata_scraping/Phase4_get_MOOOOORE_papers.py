"""
Phase 4: Targeted Cloudflare Spoofer & Queue Manager
Pre-processes the target list to isolate specific publishers (e.g. ScienceDirect) for authenticated Selenium processing.
Executes automated Cloudflare evasion on the remaining corpus to extract embedded PDF assets.
"""
import os
import csv
import re
import time
import random
from curl_cffi import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse

# Configuration & Paths
MASTER_DIR = r"C:\Users\somna\source\repos\why do i have so many directories\master audit"
INPUT_CSV = os.path.join(MASTER_DIR, "SCIENCEDIRECT_DOIS_FOR_SELENIUM.csv")
OUTPUT_DIR = os.path.join(MASTER_DIR, "CLOUDFLARE_PDFS")
SELENIUM_QUEUE_CSV = os.path.join(MASTER_DIR, "SCIENCEDIRECT_DOIS_FOR_SELENIUM.csv")

def doi_to_safe_filename(doi):
    """Sanitizes DOI strings for standard filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe)

def extract_pdf_link(html_content, base_url):
    """Parses publisher HTML to locate embedded PDF metadata tags."""
    soup = BeautifulSoup(html_content, 'html.parser')
    
    meta_tag = soup.find('meta', attrs={'name': 'citation_pdf_url'})
    if meta_tag and meta_tag.get('content'):
        return meta_tag['content']
    
    for a_tag in soup.find_all('a', href=True):
        if 'pdf' in a_tag['href'].lower() and ('download' in a_tag.text.lower() or 'pdf' in a_tag.text.lower()):
            href = a_tag['href']
            if href.startswith('/'):
                parsed_uri = urlparse(base_url)
                root = '{uri.scheme}://{uri.netloc}'.format(uri=parsed_uri)
                return root + href
            return href
            
    return None

def main():
    """Executes the pre-filtering and Cloudflare evasion pipeline."""
    if not os.path.exists(INPUT_CSV):
        return
        
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    targets = []
    skipped_for_selenium = []
    
    with open(INPUT_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            doi = row.get("doi")
            pmid = row.get("pmid", "")
            
            if not doi:
                continue
                
            doi = doi.strip()
            
            if doi.lower().startswith("10.2210/pdb"):
                continue
                
            targets.append(doi)
                
    if skipped_for_selenium:
        with open(SELENIUM_QUEUE_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["pmid", "doi"])
            writer.writeheader()
            writer.writerows(skipped_for_selenium)
    
    session = requests.Session(impersonate="chrome110")
    
    for doi in targets:
        safe_name = doi_to_safe_filename(doi)
        pdf_path = os.path.join(OUTPUT_DIR, f"{safe_name}.pdf")
        
        if os.path.exists(pdf_path):
            continue
            
        doi_url = f"https://doi.org/{doi}"
        
        try:
            response = session.get(doi_url, allow_redirects=True, timeout=30)
            
            if response.status_code == 200:
                pdf_link = extract_pdf_link(response.text, response.url)
                
                if pdf_link:
                    pdf_response = session.get(pdf_link, timeout=60)
                    
                    if pdf_response.status_code == 200 and b'%PDF' in pdf_response.content[:10]:
                        with open(pdf_path, "wb") as pdf_file:
                            pdf_file.write(pdf_response.content)
                            
        except Exception:
            pass
            
        sleep_time = random.uniform(3.0, 7.0)
        time.sleep(sleep_time)

if __name__ == "__main__":
    main()
