"""
Phase 1: Basic Accession Mapper
Maps protein accessions to PubMed IDs (PMIDs) utilizing the NCBI Entrez API.
Implements chunked API requests to handle large input datasets efficiently.
"""
import os
import time
import pandas as pd
from Bio import Entrez
from tqdm import tqdm

# Configure Entrez credentials
Entrez.email = "aditya.ghosh@mail.utoronto.ca"

INPUT_FILE = "unique_s3_projects.txt"
OUTPUT_FILE = "s3_pmid_mapping.csv"
CHUNK_SIZE = 250 #fk this took forever

def main():
    """Executes the accession to PMID mapping pipeline."""
    # Load unique accessions from input file
    with open(INPUT_FILE, "r") as f:
        accessions = [line.strip() for line in f if line.strip()]

    # Check for existing progress to enable resuming interrupted runs
    completed_accessions = set()
    if os.path.exists(OUTPUT_FILE):
        try:
            existing_df = pd.read_csv(OUTPUT_FILE, on_bad_lines='skip')
            completed_accessions = set(existing_df['Accession'].astype(str))
        except Exception:
            pass
    else:
        with open(OUTPUT_FILE, "w") as f:
            f.write("Accession,PMID\n")

    accessions_to_process = [acc for acc in accessions if acc not in completed_accessions]
    chunks = [accessions_to_process[i:i + CHUNK_SIZE] for i in range(0, len(accessions_to_process), CHUNK_SIZE)]

    for chunk in tqdm(chunks, desc="Processing Batches"):
        try:
            # Retrieve UIDs for the current accession chunk
            search_term = " OR ".join([f"{acc}[accn]" for acc in chunk])
            search_handle = Entrez.esearch(db="protein", term=search_term, retmax=CHUNK_SIZE)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            uids = search_results.get("IdList", [])
            
            if not uids:
                with open(OUTPUT_FILE, "a") as f:
                    for acc in chunk:
                        f.write(f"{acc},None\n")
                continue

            # Map protein UIDs to PubMed UIDs
            link_handle = Entrez.elink(dbfrom="protein", db="pubmed", id=",".join(uids))
            link_results = Entrez.read(link_handle)
            link_handle.close()

            # Map UIDs back to original accessions
            summary_handle = Entrez.esummary(db="protein", id=",".join(uids))
            summary_results = Entrez.read(summary_handle)
            summary_handle.close()

            uid_to_acc = {str(doc["Id"]): doc.get("AccessionVersion", "Unknown") for doc in summary_results}
            results_to_save = {}

            # Extract mappings
            for linkset in link_results:
                uid = str(linkset["IdList"][0])
                acc = uid_to_acc.get(uid)
                
                if acc:
                    pmids = []
                    if linkset.get("LinkSetDb"):
                        for link in linkset["LinkSetDb"][0]["Link"]:
                            pmids.append(link["Id"])
                    results_to_save[acc] = pmids
                    
            # Append results to CSV
            with open(OUTPUT_FILE, "a") as f:
                for acc in chunk:
                    if acc in results_to_save and results_to_save[acc]:
                        f.write(f"{acc},{';'.join(results_to_save[acc])}\n")
                    else:
                        f.write(f"{acc},None\n")

        except Exception:
            time.sleep(5)
            
        time.sleep(0.35)

if __name__ == "__main__":
    main()
