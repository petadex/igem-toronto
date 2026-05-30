"""
Phase 1: Robust Accession Mapper (Google Colab Optimized)
Maps protein accessions to PubMed IDs utilizing NCBI Entrez API JSON payloads.
Includes connection stability enforcement, offline pre-filtering, and auto-healing for invalid IDs.
(I just needed to check if the database really only have 6.6k papers from the accessions causes thats crazy)
"""
import os
import time
import http.client
import re  
import json
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
from google.colab import drive

# Enforce HTTP/1.0 to prevent NCBI mid-transfer connection drops
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

# Configure Entrez credentials
Entrez.email = "aditya.ghosh@mail.utoronto.ca"

def setup_environment():
    """Mounts Google Drive and initializes file paths."""
    drive.mount('/content/drive', force_remount=True)
    base_path = "/content/drive/MyDrive/petadex"
    os.makedirs(base_path, exist_ok=True)
    
    return {
        "input": f"{base_path}/unique_s3_projects.txt",
        "output": f"{base_path}/s3_pmid_mapping.csv",
        "error_log": f"{base_path}/s3_error_log.txt",
        "chunk_size": 200,
        "max_retries": 3
    }

def main():
    """Executes the robust batch mapping pipeline."""
    config = setup_environment()
    
    with open(config["input"], "r") as f:
        raw_accessions = [line.strip() for line in f if line.strip()]

    completed_accessions = set()
    if os.path.exists(config["output"]):
        try:
            existing_df = pd.read_csv(config["output"], on_bad_lines='skip')
            completed_accessions = set(existing_df['Accession'].astype(str))
        except Exception:
            pass
    else:
        with open(config["output"], "w") as f:
            f.write("Accession,PMID\n")

    accessions_to_process = [acc for acc in raw_accessions if acc not in completed_accessions]

    # Filter incompatible biological identifiers prior to API submission
    bad_pattern = re.compile(r"^[MHNE][A-Z]{2}\d+(\.\d+)?$|^[A-Z0-9]{10}$|^[A-Z0-9]{6}$")
    clean_accessions = []
    
    with open(config["output"], "a") as f:
        for acc in accessions_to_process:
            if bad_pattern.match(acc):
                f.write(f"{acc},INVALID_ID_FORMAT\n")
            else:
                clean_accessions.append(acc)

    chunks = [clean_accessions[i:i + config["chunk_size"]] for i in range(0, len(clean_accessions), config["chunk_size"])]

    for chunk in tqdm(chunks, desc="Processing Batches"):
        current_chunk = chunk.copy() 
        retries = 0  
        
        while current_chunk:
            try:
                # Fetch JSON Summary Data
                summary_handle = Entrez.esummary(db="protein", id=",".join(current_chunk), retmode="json")
                summary_data = json.load(summary_handle)
                summary_handle.close()
                
                uid_to_acc = {}
                uids = summary_data.get("result", {}).get("uids", [])
                
                if not uids:
                    with open(config["output"], "a") as f:
                        for acc in current_chunk:
                            f.write(f"{acc},None\n")
                    break 

                for uid in uids:
                    doc = summary_data["result"][uid]
                    acc_ver = doc.get("accessionversion", "Unknown")
                    acc_base = acc_ver.split('.')[0] if '.' in acc_ver else acc_ver
                    uid_to_acc[uid] = (acc_ver, acc_base)

                # Map UIDs to PubMed IDs via JSON structure
                link_handle = Entrez.elink(dbfrom="protein", db="pubmed", id=",".join(uids), retmode="json")
                link_data = json.load(link_handle)
                link_handle.close()

                results_to_save = {}
                for linkset in link_data.get("linksets", []):
                    if not linkset.get("ids"):
                        continue
                        
                    uid = str(linkset["ids"][0])
                    acc_tuple = uid_to_acc.get(uid)
                    
                    if acc_tuple:
                        acc_ver, acc_base = acc_tuple
                        pmids = set() 
                        
                        if linkset.get("linksetdbs"):
                            for linkset_db in linkset["linksetdbs"]:
                                if linkset_db.get("linkname") == "protein_pubmed":
                                    for pmid in linkset_db.get("links", []):
                                        pmids.add(str(pmid))
                        
                        results_to_save[acc_ver] = list(pmids)
                        results_to_save[acc_base] = list(pmids)
                        
                with open(config["output"], "a") as f:
                    for acc in current_chunk:
                        if acc in results_to_save and results_to_save[acc]:
                            f.write(f"{acc},{';'.join(results_to_save[acc])}\n")
                        else:
                            f.write(f"{acc},None\n")
                        
                break

            except RuntimeError as e:
                error_msg = str(e)
                # Trap and isolate individual rejected sequence IDs
                if "Invalid uid" in error_msg:
                    match = re.search(r"Invalid uid ([^\s]+)", error_msg)
                    if match:
                        bad_id = match.group(1)
                        if bad_id in current_chunk:
                            current_chunk.remove(bad_id)
                            with open(config["output"], "a") as f:
                                f.write(f"{bad_id},INVALID_ID\n")
                        continue 
                raise e
                
            except Exception as e:
                # Standard error handler for network timeouts
                if retries < config["max_retries"]:
                    retries += 1
                    time.sleep(5)
                    continue
                else:
                    with open(config["output"], "a") as f:
                        for acc in current_chunk:
                            f.write(f"{acc},API_ERROR\n")
                            
                    with open(config["error_log"], "a") as err_file:
                        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
                        err_file.write(f"[{timestamp}] ERROR: {type(e).__name__} {str(e)}\n")
                        err_file.write(f"Affected Accessions: {current_chunk[0]} to {current_chunk[-1]}\n\n")
                    
                    time.sleep(5)
                    break
            
        time.sleep(0.1)

if __name__ == "__main__":
    main()
