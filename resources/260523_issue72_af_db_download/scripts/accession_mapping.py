import json
import time
import sys
import requests

UNIPROT_API_URL = "https://rest.uniprot.org/idmapping"
DEFAULT_INPUT_FILE = "accessions/clean_pdb_only.csv"
DEFAULT_FROM_DB = "PDB"
DEFAULT_OUTPUT_FILE = "accessions/pdb_to_uniprot.tsv"
DEFAULT_LOG_PREFIX = "mapping_run"

def log_message(message, log_path, also_print=False):
    with open(log_path, "a", encoding="utf-8") as log_f:
        log_f.write(message + "\n")
    if also_print:
        print(message)

def status_progress_summary(status_data):
    progress = status_data.get("progress")
    processed = status_data.get("processed") or status_data.get("processedCount")
    total = status_data.get("total") or status_data.get("totalCount")

    if progress is not None:
        return f"progress={progress}%"

    if processed is not None and total:
        try:
            percent = (float(processed) / float(total)) * 100
            return f"progress={percent:.1f}% ({processed}/{total})"
        except (TypeError, ValueError, ZeroDivisionError):
            return f"processed={processed}, total={total}"

    return None

def submit_mapping_job(from_db, ids_list, log_path):
    """Submits a list of IDs to the UniProt mapping queue and returns the Job ID."""
    payload = {
        "from": from_db,
        "to": "UniProtKB",
        "ids": ",".join(ids_list)
    }
    response = requests.post(f"{UNIPROT_API_URL}/run", data=payload)
    if response.status_code == 200:
        return response.json()["jobId"]
    else:
        error_path = f"{log_path}.submit_error.txt"
        with open(error_path, "w", encoding="utf-8") as err_f:
            err_f.write(response.text)
        log_message(f"[-] Error submitting job: HTTP {response.status_code}. Saved body to {error_path}", log_path, True)
        return None

def check_job_status(job_id, log_path):
    """Polls the UniProt server until the mapping job finishes."""
    while True:
        response = requests.get(f"{UNIPROT_API_URL}/status/{job_id}")
        status_data = response.json()
        progress_summary = status_progress_summary(status_data)
        
        if status_data.get("jobStatus") == "RUNNING":
            if progress_summary:
                log_message(f"    -> Job still processing ({progress_summary}), waiting 60 seconds...", log_path, False)
            else:
                log_message("    -> Job still processing, waiting 60 seconds...", log_path, False)
            time.sleep(60)
        elif status_data.get("jobStatus") == "FINISHED" or "results" in response.links:
            return status_data
        elif "results" in status_data:
            return status_data
        else:
            summary = {key: status_data.get(key) for key in ("jobStatus", "warnings", "failedIds") if key in status_data}
            log_message(f"[-] Job failed or encountered an error: {summary}", log_path, True)
            return False

def extract_uniprot_accession(to_value):
    if isinstance(to_value, dict):
        for key in ("primaryAccession", "uniProtkbId", "id", "accession"):
            value = to_value.get(key)
            if value:
                return value
        return None
    if isinstance(to_value, str):
        return to_value
    return None

def write_results_pairs(results, output_handle):
    mapped_count = 0
    for pair in results:
        from_id = pair.get("from")
        to_id = extract_uniprot_accession(pair.get("to"))
        if from_id and to_id:
            output_handle.write(f"{from_id}\t{to_id}\n")
            mapped_count += 1
    return mapped_count

def stream_mapping_pairs(job_id, output_handle, debug_dir=".", log_path=None):
    """Streams mapping pairs from UniProt and writes 'From' -> 'Entry' rows."""
    response = requests.get(
        f"{UNIPROT_API_URL}/stream/{job_id}",
        params={"format": "tsv"},
        stream=True,
        timeout=120,
    )

    if response.status_code != 200:
        if log_path:
            log_message(f"[-] Error retrieving results: HTTP {response.status_code}", log_path, True)
        else:
            print(f"[-] Error retrieving results: HTTP {response.status_code}")
        error_path = f"{debug_dir}/mapping_error_{job_id}.txt"
        with open(error_path, "w", encoding="utf-8") as err_f:
            err_f.write(response.text)
        if log_path:
            log_message(f"    -> Saved error body to {error_path}", log_path, True)
        else:
            print(f"    -> Saved error body to {error_path}")
        return 0

    mapped_count = 0
    header = None
    from_idx = None
    to_idx = None
    debug_path = f"{debug_dir}/mapping_response_{job_id}.tsv"
    debug_f = open(debug_path, "w", encoding="utf-8")
    for raw_line in response.iter_lines(decode_unicode=True):
        if not raw_line:
            continue

        debug_f.write(raw_line + "\n")

        if header is None:
            header = [col.strip().lower() for col in raw_line.split("\t")]
            if "from" in header:
                from_idx = header.index("from")
            for candidate in ("entry", "to", "uniprotkb", "uniprotkb_ac", "accession", "primaryaccession"):
                if candidate in header:
                    to_idx = header.index(candidate)
                    break

            # If header row doesn't include columns, keep scanning
            if from_idx is not None and to_idx is not None:
                continue

        parts = raw_line.split("\t")
        if from_idx is None or to_idx is None:
            continue
        if len(parts) <= max(from_idx, to_idx):
            continue

        from_id = parts[from_idx].strip()
        to_id = parts[to_idx].strip()
        if from_id and to_id:
            output_handle.write(f"{from_id}\t{to_id}\n")
            mapped_count += 1
    debug_f.close()
    if log_path:
        log_message(f"    -> Saved raw response to {debug_path}", log_path, True)
    else:
        print(f"    -> Saved raw response to {debug_path}")

    return mapped_count


def run_bulk_mapping(input_csv, from_db, output_tsv, batch_size=50000):
    log_path = f"{DEFAULT_LOG_PREFIX}_{int(time.time())}.log"
    log_message(f"[*] Reading IDs from {input_csv}...", log_path, True)
    
    # Read IDs, accounting for commas, tabs, and stripping whitespace
    ids = []
    with open(input_csv, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Extract first token before a comma or tab
            first_val = line.split(',')[0].split('\t')[0].strip()
            # Skip header lines safely
            if first_val.lower() not in ["accession", "id", "from", "entry"]:
                ids.append(first_val)
                
    log_message(f"[+] Successfully extracted {len(ids):,} total IDs from your input file.", log_path, True)
    
    if len(ids) == 0:
        log_message("[-] Stopping: No valid source IDs found to process. Check your input file format.", log_path, True)
        return
    
    log_message(f"[+] Found {len(ids):,} total IDs to map.", log_path, True)
    
    # Open file with .tsv extension and write the matching headers
    with open(output_tsv, 'w') as out_f:
        out_f.write("From\tEntry\n")  # Matches web layout exactly
        
        for i in range(0, len(ids), batch_size):
            chunk = ids[i:i + batch_size]
            log_message(f"[*] Submitting batch {(i//batch_size)+1} ({len(chunk):,} IDs)...", log_path, True)
            
            job_id = submit_mapping_job(from_db, chunk, log_path)
            if not job_id:
                continue
                
            log_message(f"    -> Job submitted. ID: {job_id}. Polling status...", log_path, True)
            status_data = check_job_status(job_id, log_path)
            if status_data:
                details = requests.get(f"{UNIPROT_API_URL}/details/{job_id}")
                if details.status_code == 200:
                    details_path = f"mapping_details_{job_id}.json"
                    with open(details_path, "w", encoding="utf-8") as details_f:
                        json.dump(details.json(), details_f, indent=2)
                    log_message(f"    -> Saved job details to {details_path}", log_path, True)

                mapped_count = stream_mapping_pairs(job_id, out_f, log_path=log_path)
                log_message(f"[+] Mapped and saved {mapped_count:,} pairs from this batch.", log_path, True)
            
            time.sleep(5)

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT_FILE
    from_db = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_FROM_DB
    output_file = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_OUTPUT_FILE
    batch_size = int(sys.argv[4]) if len(sys.argv) > 4 else 50000

    # Usage: python accession_mapping.py [input_csv] [from_db] [output_tsv] [batch_size]

    run_bulk_mapping(input_file, from_db, output_file, batch_size=batch_size)