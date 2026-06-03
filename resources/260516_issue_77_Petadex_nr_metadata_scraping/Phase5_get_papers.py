"""
Phase 5: Scihub scraper #DaGoat
Utilizes undetected_chromedriver and PyAutoGUI to navigate target repositories and simulate user interactions for PDF retrieval.
Implements robust error handling and browser state recovery for uninterrupted batch processing.
"""
import os
import time
import csv
import pyautogui #basically made my computer unusable for like 2 days
from selenium.webdriver.common.by import By
import undetected_chromedriver as uc
import re

# Configuration & Paths
DOWNLOAD_DIR = r"C:\Users\somna\Downloads\New_PDFs"
CSV_PATH = r"C:\Users\somna\source\repos\scraping logan attempt 16\scraping logan attempt 16\data\custom_scraper_untried_dois.csv"
LOG_CSV_PATH = r"C:\Users\somna\source\repos\scraping logan attempt 16\scraping logan attempt 16\data\missing_dois_log1.csv"
TARGET_HOST = "sci-hub.box"

def main():
    """Executes the automated UI scraping pipeline."""
    os.makedirs(DOWNLOAD_DIR, exist_ok=True)

    if not os.path.exists(LOG_CSV_PATH):
        with open(LOG_CSV_PATH, mode='w', newline='', encoding='utf-8') as log_file:
            log_writer = csv.writer(log_file)
            log_writer.writerow(["Missing_DOI"])

    dois_to_scrape = []
    with open(CSV_PATH, mode='r', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        next(reader, None)  
        for row in reader:
            if len(row) > 1:
                dois_to_scrape.append(row[1].strip())

    doi_retry_counts = {}

    while True:
        already_downloaded = set()
        for filename in os.listdir(DOWNLOAD_DIR):
            if filename.endswith(".pdf"):
                already_downloaded.add(filename)

        already_failed = set()
        if os.path.exists(LOG_CSV_PATH):
            with open(LOG_CSV_PATH, mode='r', encoding='utf-8') as file:
                reader = csv.reader(file)
                next(reader, None)
                for row in reader:
                    if row:
                        clean_doi = row[0].split(" (ERROR")[0].strip()
                        already_failed.add(clean_doi)

        remaining_dois = [doi for doi in dois_to_scrape if re.sub(r'[\\/*?:"<>|]', '_', doi) + ".pdf" not in already_downloaded and doi not in already_failed]
        
        if not remaining_dois:
            break

        prefs = {
            "download.default_directory": DOWNLOAD_DIR,
            "download.prompt_for_download": False,
            "download.directory_upgrade": True,
            "plugins.always_open_pdf_externally": False
        }

        options = uc.ChromeOptions()
        options.add_experimental_option("prefs", prefs)
        
        current_doi = None  

        try:
            driver = uc.Chrome(version_main=148, options=options)

            screen_width, screen_height = pyautogui.size()
            click_x = int(screen_width * 0.95)
            click_y = int(screen_height * 0.12)
            neutral_x = int(screen_width * 0.10)
            neutral_y = int(screen_height * 0.50)

            for doi in remaining_dois:
                current_doi = doi 
                
                safe_doi = re.sub(r'[\\/*?:"<>|]', '_', doi)
                expected_filename = f"{safe_doi}.pdf"
                save_path = os.path.join(DOWNLOAD_DIR, expected_filename)
                
                target_url = f"https://{TARGET_HOST}/{doi}"
                driver.get(target_url)

                missing_paper = False
                timeout = 30
                start_time = time.time()
                page_loaded = False
                captcha_clicked = False

                while time.time() - start_time < timeout:
                    # Handle specific CAPTCHA events
                    if "robot" in driver.title.lower() and not captcha_clicked:
                        pyautogui.moveTo(192, 421, duration=0.5) 
                        pyautogui.click()
                        captcha_clicked = True 
                        time.sleep(2) 
                        continue 
                    
                    # Verify document availability in the repository
                    error_elements = driver.find_elements(By.CSS_SELECTOR, "body > fixed-width > column-split > block-rounded")
                    if error_elements and "Alas, the following paper is not yet available" in error_elements[0].text:
                        missing_paper = True
                        page_loaded = True
                        break
                    
                    if "no articles found" in driver.title:
                        missing_paper = True
                        page_loaded = True
                        break

                    if "Sci" in driver.title:
                        page_loaded = True
                        break
                    
                    time.sleep(1)

                if not page_loaded:
                    raise Exception("Browser security check timed out.")
                
                if not missing_paper:
                    time.sleep(2)

                if missing_paper:
                    with open(LOG_CSV_PATH, mode='a', newline='', encoding='utf-8') as log_file:
                        log_writer = csv.writer(log_file)
                        log_writer.writerow([doi])
                    continue

                # Execute UI interaction sequence to initiate download
                pyautogui.moveTo(click_x, click_y, duration=0.5)
                pyautogui.click()
                time.sleep(2) 

                pyautogui.write(save_path, interval=0.02)
                time.sleep(1)
                pyautogui.press('enter')

                file_saved = False
                i = 1
                start_time = time.time()

                for attempt in range(4):
                    time.sleep(4) 
                    if os.path.exists(save_path):
                        file_saved = True
                        break
                    
                    if i == 2:
                        pyautogui.hotkey('ctrl', 'r')
                        while time.time() - start_time < timeout:
                            if "Sci" in driver.title:
                                page_loaded = True
                                time.sleep(5)
                                break
                    
                    pyautogui.moveTo(click_x, click_y, duration=0.5)
                    pyautogui.click()
                    time.sleep(2)
                    
                    pyautogui.moveTo(1124, 535, duration=0.2)
                    pyautogui.click()
                    time.sleep(0.5)
                    pyautogui.hotkey('ctrl', 'a')
                    time.sleep(0.5)
                    pyautogui.press('backspace')
                    time.sleep(0.5)
                    
                    pyautogui.write(save_path, interval=0.05)
                    time.sleep(1)
                    pyautogui.press('enter')
         
                    time.sleep(4)
                    i += 1
                    
                if not file_saved:
                    pyautogui.press('esc')
                    raise Exception("File resolution timeout.")

                pyautogui.moveTo(neutral_x, neutral_y, duration=0.2)
                pyautogui.click()
                time.sleep(1)

            driver.quit()
            break 

        except Exception as e:
            if current_doi:
                doi_retry_counts[current_doi] = doi_retry_counts.get(current_doi, 0) + 1
                
                if doi_retry_counts[current_doi] >= 2:
                    with open(LOG_CSV_PATH, mode='a', newline='', encoding='utf-8') as log_file:
                        log_writer = csv.writer(log_file)
                        log_writer.writerow([f"{current_doi} (ERROR: Failed 2 attempts)"])
            
            try:
                driver.quit()
            except:
                pass
            
            time.sleep(5) 
            continue

if __name__ == "__main__":
    main()
