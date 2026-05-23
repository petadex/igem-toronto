import re
import sys
import time

# Define regex patterns optimized for fast compilation

# Matches standard 6-character and 10-character UniProtKB primary accessions
UNIPROT_PAT = re.compile(r'^([A-N_R-Z][0-9][A-Z0-9]{3}[0-9]|[O,P,Q][0-9][A-Z0-9]{4})(?:[A-Z0-9]{4})?$')
# Matches RefSeq prefixes (WP_, NP_, XP_, YP_) followed by any number of digits
REFSEQ_PAT = re.compile(r'^(WP|NP|XP|YP)_\d+$')
# Matches GenBank protein accessions (3 letters followed by 5 or 7 digits)
GENBANK_PAT = re.compile(r'^[A-Z]{3}\d{5,7}$')

def split_accessions(input_csv_path):
    print(f"[*] Initializing processing for: {input_csv_path}")
    start_time = time.time()
    
    counters = {"uniprot": 0, "refseq": 0, "genbank": 0, "other": 0, "total": 0}
    
    try:
        with open(input_csv_path, 'r', encoding='utf-8') as infile, \
            open('accessions/uniprot_accessions.csv', 'w', encoding='utf-8') as f_uni, \
             open('accessions/refseq_accessions.csv', 'w', encoding='utf-8') as f_ref, \
             open('accessions/genbank_accessions.csv', 'w', encoding='utf-8') as f_gen, \
             open('accessions/other_accessions.csv', 'w', encoding='utf-8') as f_oth:
            # Read first line to check for header
            first_line = infile.readline()
            if not first_line:
                print("[-] Error: The input file is empty.")
                return
            
            # Detect header row
            if any(x in first_line.lower() for x in ["accession", "id", "protein", "query"]):
                f_uni.write(first_line)
                f_ref.write(first_line)
                f_gen.write(first_line)
                f_oth.write(first_line)
                print("[*] Detected header row. Copied header across all output files.")
            else:
                # Process the first line as a data row if it's not a header
                process_row(first_line, f_uni, f_ref, f_gen, f_oth, counters)
                
            # Iterate through the rest of the 4.7 million rows using low-memory streaming
            for line in infile:
                process_row(line, f_uni, f_ref, f_gen, f_oth, counters)
                
                # Progress logging loop
                if counters["total"] % 500000 == 0:
                    print(f"    -> Processed {counters['total']:,} rows...")
                    
    except FileNotFoundError:
        print(f"[-] Error: File '{input_csv_path}' not found. Please verify the filename.")
        return

    elapsed_time = time.time() - start_time
    print("[+] Processing Completed Successfully!")
    print(f"Time Taken: {elapsed_time:.2f} seconds")
    print("-" * 40)
    print(f"Total Rows Parsed:     {counters['total']:,}")
    print(f"├── UniProt IDs:       {counters['uniprot']:,} -> written to 'uniprot_accessions.csv'")
    print(f"├── RefSeq IDs:        {counters['refseq']:,} -> written to 'refseq_accessions.csv'")
    print(f"├── GenBank IDs:       {counters['genbank']:,} -> written to 'genbank_accessions.csv'")
    print(f"└── Unclassified/Other: {counters['other']:,} -> written to 'other_accessions.csv'")
    print("-" * 40)

def process_row(line, f_uni, f_ref, f_gen, f_oth, counters):
    raw_val = line.strip()
    if not raw_val:
        return
        
    # Isolate first column if line contains multiple comma-separated fields
    base_col = raw_val.split(',')[0]
    
    # Clean the ID: Strip out version suffixes (e.g., .1, .2) and whitespaces
    clean_id = base_col.split('.')[0].strip()
    
    counters["total"] += 1
    
    # Route the original line to the matching bucket using precompiled regex
    if UNIPROT_PAT.match(clean_id):
        f_uni.write(line)
        counters["uniprot"] += 1
    elif REFSEQ_PAT.match(clean_id):
        f_ref.write(line)
        counters["refseq"] += 1
    elif GENBANK_PAT.match(clean_id):
        f_gen.write(line)
        counters["genbank"] += 1
    else:
        f_oth.write(line)
        counters["other"] += 1

if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "clean_accessions_no_version.csv"
    split_accessions(input_file)

