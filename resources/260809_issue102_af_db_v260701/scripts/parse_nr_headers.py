#!/usr/bin/env python3
"""Parse PETadex v260701 nr significant-hits FASTA headers into record-linked tables.

Run 1 (issue 72, May 2026) could get away with a one-line awk extraction because every
FASTA record carried exactly one accession. The nr hits do not: a single sequence in nr is
shared by several accessions, so a header looks like

    >WP_453139269.1,MGW6465901.1|3.06e-50|42.9|181.0|PL127|C014
     ^-- accessions --------------^ e_value  pid bits  query   component

and "does this sequence already have an AlphaFold structure?" becomes an OR over every
accession on the record. That requires keeping the record -> accession link, which is what
this script produces.

Outputs (all under --outdir):

  records.tsv           rec_id, query_petadex_id, component_id, e_value, pid, bitscore,
                        n_accessions
  record_accessions.tsv rec_id, accession_raw, accession_clean   (one row per pair)
  unique_accessions.csv accession_clean, deduplicated -- the input to split_accession_types.py
  parse_stats.json      counts for the notebook

rec_id is the 1-based ordinal of the '>' line in the file. It is deliberately NOT
query_petadex_id|component_id: many nr hits share the same query, so that pair is not unique.
Every downstream script keys off rec_id and re-derives it by streaming the same FASTA in the
same order, so the input file must not be re-sorted between steps.

Usage:
    python parse_nr_headers.py INPUT.fasta --outdir out/ [--sep ,] [--limit N]
"""

import argparse
import gzip
import io
import json
import os
import re
import sys
import time

# Strip a trailing ".1" / ".2" version suffix (NCBI accessions always carry one).
VERSION_RE = re.compile(r"\.\d+$")
# Strip a trailing single-letter PDB chain suffix ("2ZPQ_A" -> "2ZPQ"). Deliberately a single
# letter: a broader "_[A-Za-z0-9]+$" would eat the numeric body of every RefSeq accession
# ("WP_305403363" -> "WP"), which is exactly the bug carried by the Run 1 mapping script.
CHAIN_RE = re.compile(r"_[A-Za-z]$")


def open_fasta(path):
    """Open a FASTA as text, transparently handling .zst / .gz / plain."""
    if path.endswith(".zst"):
        import zstandard

        raw = open(path, "rb")
        reader = zstandard.ZstdDecompressor().stream_reader(raw)
        return io.TextIOWrapper(reader, encoding="utf-8", errors="replace")
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace", buffering=1 << 22)


def clean_accession(acc):
    """Normalize an accession to the form used as a join key everywhere downstream."""
    acc = acc.strip()
    acc = VERSION_RE.sub("", acc)
    acc = CHAIN_RE.sub("", acc)
    return acc


def parse(fasta_path, outdir, sep=",", limit=None, progress_every=1_000_000):
    os.makedirs(outdir, exist_ok=True)
    start = time.time()

    rec_id = 0
    n_pairs = 0
    n_malformed = 0
    max_accs = 0
    accs_hist = {}
    unique = set()

    f_rec = open(os.path.join(outdir, "records.tsv"), "w", encoding="utf-8", newline="\n")
    f_pair = open(os.path.join(outdir, "record_accessions.tsv"), "w", encoding="utf-8", newline="\n")
    f_rec.write("rec_id\tquery_petadex_id\tcomponent_id\te_value\tpid\tbitscore\tn_accessions\n")
    f_pair.write("rec_id\taccession_raw\taccession_clean\n")

    with open_fasta(fasta_path) as fin, f_rec, f_pair:
        for line in fin:
            if not line.startswith(">"):
                continue
            rec_id += 1

            fields = line[1:].rstrip("\n").rstrip("\r").split("|")
            acc_field = fields[0]
            # Pad so a short/odd header still produces a row rather than crashing the run.
            if len(fields) < 6:
                n_malformed += 1
                fields = fields + [""] * (6 - len(fields))
            e_value, pid, bitscore, query_id, component_id = fields[1:6]

            raw_accs = [a for a in acc_field.split(sep) if a.strip()]
            n = len(raw_accs)
            max_accs = max(max_accs, n)
            accs_hist[n] = accs_hist.get(n, 0) + 1

            f_rec.write(
                f"{rec_id}\t{query_id}\t{component_id}\t{e_value}\t{pid}\t{bitscore}\t{n}\n"
            )
            for raw in raw_accs:
                raw = raw.strip()
                cleaned = clean_accession(raw)
                if not cleaned:
                    continue
                f_pair.write(f"{rec_id}\t{raw}\t{cleaned}\n")
                unique.add(cleaned)
                n_pairs += 1

            if rec_id % progress_every == 0:
                rate = rec_id / (time.time() - start)
                print(f"    -> {rec_id:,} records, {n_pairs:,} pairs ({rate:,.0f} rec/s)", flush=True)
            if limit and rec_id >= limit:
                break

    uniq_path = os.path.join(outdir, "unique_accessions.csv")
    with open(uniq_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\n")
        for acc in sorted(unique):
            f.write(acc + "\n")

    elapsed = time.time() - start
    stats = {
        "input": fasta_path,
        "separator": sep,
        "records": rec_id,
        "accession_pairs": n_pairs,
        "unique_accessions": len(unique),
        "max_accessions_per_record": max_accs,
        "mean_accessions_per_record": round(n_pairs / rec_id, 4) if rec_id else 0,
        "accessions_per_record_hist": {str(k): v for k, v in sorted(accs_hist.items())},
        "malformed_headers": n_malformed,
        "seconds": round(elapsed, 1),
    }
    with open(os.path.join(outdir, "parse_stats.json"), "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print("[+] Parsing complete.")
    print(f"    records                  : {rec_id:,}")
    print(f"    (record, accession) pairs: {n_pairs:,}")
    print(f"    unique accessions        : {len(unique):,}")
    print(f"    max accessions / record  : {max_accs}")
    print(f"    malformed headers        : {n_malformed:,}")
    print(f"    elapsed                  : {elapsed:,.1f}s")
    print(f"[+] Wrote records.tsv, record_accessions.tsv, unique_accessions.csv to {outdir}/")
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="input FASTA (.fastaa / .fasta / .fa / .gz / .zst)")
    ap.add_argument("--outdir", required=True, help="directory for the output tables")
    ap.add_argument(
        "--sep",
        default=",",
        help="separator between accessions inside header field 1 (default: ',' -- verified "
        "against the v260701 file; pass ':' if a future build changes it)",
    )
    ap.add_argument("--limit", type=int, default=None, help="stop after N records (for pilots)")
    args = ap.parse_args()

    if not os.path.exists(args.fasta):
        sys.exit(f"[-] Input not found: {args.fasta}")
    parse(args.fasta, args.outdir, sep=args.sep, limit=args.limit)


if __name__ == "__main__":
    main()
