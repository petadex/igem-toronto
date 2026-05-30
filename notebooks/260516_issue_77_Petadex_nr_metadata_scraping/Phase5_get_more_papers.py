"""
Phase 5: Authenticated Institutional Proxy Fetcher
Automates the retrieval of full-text PDFs through an institutional EZProxy portal using an undetected Chrome WebDriver.
Features manual authentication handoff, adaptive human pacing (jitter), and scheduled pauses (coffee breaks) to evade automated bot detection rate limits.
"""
import os
import time
import csv
import re
import random
import pyautogui
import undetected_chromedriver as uc
from selenium.webdriver.common.by import By

# Configuration & Paths
INPUT_CSV = r"C:\Users\somna\source\repos\scraping logan attempt 16\scraping logan attempt 16\data\ABSOLUTE_FINAL_MISSING_DOIS.csv"
DOWNLOAD_DIR = r"C:\Users\somna\source\repos\why do i have so many directories\master audit\FINAL_PROXY_PDFS"
ERROR_LOG_CSV = r"C:\Users\somna\source\repos\why do i have so many directories\master audit\PROXY_ERROR_LOG.csv"

# Institutional EZProxy Prefix
PROXY_PREFIX = "http://myaccess.library.utoronto.ca/login?url=https://doi.org/"

def doi_to_safe_filename(doi):
    """Sanitizes DOI strings for local filesystem storage."""
    safe = doi.strip().lower().replace("/", "_")
    return re.sub(r'[<>:"|?*\\]', "_", safe)

def log_failed_doi(doi, reason):
    """Logs unretrievable DOIs and associated error messages to a persistent CSV."""
    file_exists = os.path.isfile(ERROR_LOG_CSV)
    with open(ERROR_LOG_CSV, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(["doi", "error_reason"])
        writer.writerow([doi, reason])

def main():
    """Executes the authenticated Selenium fetching pipeline."""
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)

    failed_dois = set()
    if os.path.exists(ERROR_LOG_CSV):
        with open(ERROR_LOG_CSV, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("doi"):
                    failed_dois.add(row["doi"].strip())

    targets = []
    if os.path.exists(INPUT_CSV):
        with open(INPUT_CSV, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                doi = row.get("doi")
                if doi and not doi.lower().startswith("10.2210/pdb"):
                    doi = doi.strip()
                    if doi not in failed_dois:
                        targets.append(doi)
    else:
        return
                    
    prefs = {
        "download.default_directory": os.path.abspath(DOWNLOAD_DIR),
        "download.prompt_for_download": False,
        "plugins.always_open_pdf_externally": True 
    }
    options = uc.ChromeOptions()
    options.add_experimental_option("prefs", prefs)
    
    driver = uc.Chrome(version_main=148, options=options)
    
    # find unsuspecting party to sign in (jkjk)
    if targets:
        driver.get(f"{PROXY_PREFIX}{targets[0]}")
        
        print("WAITING FOR MANUAL AUTHENTICATION...")
        print("1. Log in with your institutional credentials.")
        print("2. Approve the 2FA push.")
        print("3. Wait for the publisher page to finish loading.")
        print("4. Tell me to come back and return to this terminal and press ENTER.")
        print("I swear I don't have your login info lol (thank you so much)")
        input("Press ENTER when authenticated and ready to begin batch processing... ")
    
    successful = 0
    papers_processed = 0

    for i, doi in enumerate(targets, 1):
        safe_name = doi_to_safe_filename(doi)
        expected_file = os.path.join(DOWNLOAD_DIR, f"{safe_name}.pdf")

        if os.path.exists(expected_file):
            continue 
            
        papers_processed += 1
        
        try:
            # Context refresh: Lock focus to the primary page frame
            driver.switch_to.window(driver.window_handles[0])
            driver.get(f"{PROXY_PREFIX}{doi}")
            time.sleep(8) 
            
            # Check for embedded PDF metadata URL
            pdf_url = driver.execute_script("return document.querySelector('meta[name=\"citation_pdf_url\"]')?.content;")
            if pdf_url:
                driver.get(pdf_url)
                time.sleep(1)
                pyautogui.moveTo(974, 632, duration=0.2)
                pyautogui.click()
            else:
                # Fallback to text based PDF links
                pdf_links = driver.find_elements(By.XPATH, "//a[contains(translate(., 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'pdf')]")
                if pdf_links:
                    driver.execute_script("arguments[0].click();", pdf_links[0])
                else:
                    log_failed_doi(doi, "No PDF link found")
                    continue
                    
            # Tab cleanup logic
            time.sleep(2) 
            if len(driver.window_handles) > 1:
                driver.switch_to.window(driver.window_handles[1])
                driver.close()
                driver.switch_to.window(driver.window_handles[0])
                driver.switch_to.default_content()
                    
            # Dynamic download verification loop
            timeout = time.time() + 12 
            success_found = False
        
            while time.time() < timeout:
                files = [os.path.join(DOWNLOAD_DIR, f) for f in os.listdir(DOWNLOAD_DIR) if f.endswith('.pdf')]
                if files:
                    latest_file = max(files, key=os.path.getctime)
                    if time.time() - os.path.getctime(latest_file) < 12:
                        try:
                            os.rename(latest_file, expected_file)
                            successful += 1
                            success_found = True
                            break
                        except Exception:
                            pass
                time.sleep(1) 
            
            # Final fallback utilizing dynamic UI querying for obscured download buttons
            if not success_found:
                alt_buttons = driver.find_elements(By.XPATH, "//*[self::button or self::a][@id='new-download-btn' or contains(translate(@aria-label, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download') or contains(translate(@title, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download')]")
                
                if alt_buttons:
                    pyautogui.moveTo(1880, 123, duration=0.2)
                    pyautogui.click()
                    
                    pyautogui.moveTo(1797, 250, duration=0.2)
                    pyautogui.click()

                    timeout = time.time() + 12
                    while time.time() < timeout:
                        files = [os.path.join(DOWNLOAD_DIR, f) for f in os.listdir(DOWNLOAD_DIR) if f.endswith('.pdf')]
                        if files:
                            latest_file = max(files, key=os.path.getctime)
                            if time.time() - os.path.getctime(latest_file) < 12:
                                try:
                                    os.rename(latest_file, expected_file)
                                    successful += 1
                                    success_found = True
                                    break
                                except Exception:
                                    pass
                        time.sleep(1)
                        
                    pyautogui.moveTo(555, 256, duration=0.2)
                    pyautogui.click()
                else:
                    log_failed_doi(doi, "No PDF link found")
                    continue
                    
                if not success_found:
                    log_failed_doi(doi, "Timeout or download failed")
        
        except Exception as e:
            error_msg = str(e)[:100].replace('\n', ' ')
            log_failed_doi(doi, error_msg)
            
        # Pacing engine: Randomized delay to mimic human reading patterns
        jitter = random.randint(8, 20)
        time.sleep(jitter)

        # Extended pause to reset library rate limiters
        if papers_processed > 0 and papers_processed % 75 == 0:
            break_time = random.randint(300, 600) 
            time.sleep(break_time)

    driver.quit()

if __name__ == "__main__":
    main()
