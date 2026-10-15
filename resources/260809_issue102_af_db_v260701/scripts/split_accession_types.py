#!/usr/bin/env python3
"""Split a flat accession list into UniProt / RefSeq / GenBank / PDB / other buckets.

Copied from resources/260523_issue72_af_db_download/scripts/split_accession_types.py (Run 1)
and changed in three ways for Run 2 -- the Run 1 copy is left untouched so its results stay
reproducible.

  1. --outdir. Run 1 hardcoded 'accessions/*.csv' relative to the repo root, which would
     silently overwrite the Run 1 outputs still sitting there.
  2. An explicit PDB bucket. Run 1 had no PDB pattern, so PDB IDs fell into 'other' and the
     mapping step then treated all of 'other' as PDB. That was true in Run 1 (all 1,448 were
     PDB) but is not guaranteed for nr, and a silent misclassification there produces wrong
     mappings rather than an error. 'other' now means genuinely unclassified and must be
     inspected before Step 3.
  3. --header / --no-header instead of sniffing the first line. The Run 1 sniff fired when the
     first line contained any of "accession"/"id"/"protein"/"query", which a real accession
     beginning "MID..." would trigger, dropping a real row.

Usage:
    python split_accession_types.py unique_accessions.csv --outdir accessions/ --header
"""

import argparse
import json
import os
import re
import sys
import time

# Standard 6-character and 10-character UniProtKB primary accessions.
UNIPROT_PAT = re.compile(r"^([A-N_R-Z][0-9][A-Z0-9]{3}[0-9]|[O,P,Q][0-9][A-Z0-9]{4})(?:[A-Z0-9]{4})?$")
# RefSeq protein prefixes followed by digits (WP_305403363, XP_067660346, NP_..., YP_...).
REFSEQ_PAT = re.compile(r"^(WP|NP|XP|YP|AP|ZP)_\d+$")
# GenBank/EMBL protein accessions: 3 letters + 5 or 7 digits (AVI23960, HEX8796638).
GENBANK_PAT = re.compile(r"^[A-Z]{3}\d{5,7}$")
# Legacy 2-letter GenBank/EMBL protein accessions (CAA12345-style predecessors).
GENBANK_LEGACY_PAT = re.compile(r"^[A-Z]{2}\d{5,6}$")
# PDB entry IDs: digit + 3 alphanumerics (chain suffix already stripped upstream).
PDB_PAT = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

BUCKETS = ["uniprot", "refseq", "genbank", "pdb", "other"]


def classify(clean_id):
    if UNIPROT_PAT.match(clean_id):
        return "uniprot"
    if REFSEQ_PAT.match(clean_id):
        return "refseq"
    if GENBANK_PAT.match(clean_id) or GENBANK_LEGACY_PAT.match(clean_id):
        return "genbank"
    if PDB_PAT.match(clean_id):
        return "pdb"
    return "other"


def split_accessions(input_csv_path, outdir, has_header):
    print(f"[*] Initializing processing for: {input_csv_path}")
    os.makedirs(outdir, exist_ok=True)
    start_time = time.time()

    counters = {b: 0 for b in BUCKETS}
    counters["total"] = 0

    handles = {
        b: open(os.path.join(outdir, f"{b}_accessions.csv"), "w", encoding="utf-8", newline="\n")
        for b in BUCKETS
    }
    try:
        with open(input_csv_path, "r", encoding="utf-8") as infile:
            if has_header:
                header = infile.readline()
                for h in handles.values():
                    h.write(header)
                print("[*] Skipped header row (copied into every output file).")

            for line in infile:
                raw_val = line.strip()
                if not raw_val:
                    continue
                # Isolate the first column, then strip any residual version suffix.
                clean_id = raw_val.split(",")[0].split("\t")[0].split(".")[0].strip()
                counters["total"] += 1
                bucket = classify(clean_id)
                handles[bucket].write(clean_id + "\n")
                counters[bucket] += 1

                if counters["total"] % 1_000_000 == 0:
                    print(f"    -> Processed {counters['total']:,} rows...", flush=True)
    except FileNotFoundError:
        sys.exit(f"[-] Error: File '{input_csv_path}' not found.")
    finally:
        for h in handles.values():
            h.close()

    elapsed_time = time.time() - start_time
    print("[+] Processing Completed Successfully!")
    print(f"Time Taken: {elapsed_time:.2f} seconds")
    print("-" * 46)
    print(f"Total Rows Parsed:      {counters['total']:,}")
    for b in BUCKETS:
        print(f"  {b:<10} {counters[b]:>12,}  -> {outdir}/{b}_accessions.csv")
    print("-" * 46)
    print("[!] Inspect other_accessions.csv before running Step 3 -- unclassified accessions")
    print("    are NOT mapped and would silently disappear from the manifest.")

    with open(os.path.join(outdir, "split_stats.json"), "w", encoding="utf-8") as f:
        json.dump({**counters, "seconds": round(elapsed_time, 1)}, f, indent=2)
    return counters


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_csv", help="flat accession list (one per line), e.g. unique_accessions.csv")
    ap.add_argument("--outdir", required=True, help="directory for the per-type CSVs")
    ap.add_argument("--header", dest="has_header", action="store_true", help="input has a header row")
    ap.add_argument("--no-header", dest="has_header", action="store_false")
    ap.set_defaults(has_header=True)
    args = ap.parse_args()
    split_accessions(args.input_csv, args.outdir, args.has_header)


if __name__ == "__main__":
    main()
