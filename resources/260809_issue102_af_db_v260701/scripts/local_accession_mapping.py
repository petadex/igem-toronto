#!/usr/bin/env python3
"""Map RefSeq / GenBank / PDB accessions to UniProtKB ACs via UniProt's idmapping table.

Copied from resources/260523_issue72_af_db_download/scripts/local_accession_mapping.py (Run 1)
and changed in four ways. The Run 1 copy is left untouched so its results stay reproducible.

  1. FIXED: type-dependent normalization. Run 1 applied `_[A-Za-z0-9]+$` -> "" to the master
     table's source column for every type, with the comment "only applies for pdb". It does
     not: polars replaces the leftmost match, so "WP_054022242" becomes "WP" and every RefSeq
     row in the master table collapses to three distinct keys. The subsequent membership
     filter against real WP_/NP_/XP_ IDs then matches nothing. Chain stripping is now applied
     only for --type pdb, and uses `_[A-Za-z]$` (a single trailing letter) so it cannot reach
     into a numeric accession body.
  2. --type {refseq,genbank,pdb} replaces sniffing the accession class out of the input
     filename. Renaming an input file silently changed which master column was joined on.
  3. --master for the idmapping path (Run 1 hardcoded it to the current directory).
  4. The membership filter is a semi-join on a LazyFrame rather than `is_in(list(...))`.
     Run 2 has an order of magnitude more accessions than Run 1 and materializing a
     multi-million-element Python list into the query plan is the slow path.

Master table columns (0-indexed) in idmapping_selected.tab.gz:
    0 UniProtKB-AC   3 RefSeq   5 PDB   17 EMBL-CDS (GenBank protein)

Usage:
    python local_accession_mapping.py IN.csv OUT.tsv --type refseq \
        --master /mnt/afdb/idmapping_selected.tab.gz
"""

import argparse
import os
import sys

import polars as pl

# Master-table column holding each accession class.
SOURCE_COL = {"refseq": 3, "genbank": 17, "pdb": 5}
TARGET_COL = 0  # UniProtKB-AC


def build_mapping_dict(master_mapping_gz, acc_type, target_ids):
    source_col = SOURCE_COL[acc_type]
    print(f"[*] Scanning master mapping file (source col {source_col} = {acc_type})...")

    source_name = f"source_col_{source_col}"
    target_name = f"target_col_{TARGET_COL}"

    q = (
        pl.scan_csv(
            master_mapping_gz,
            separator="\t",
            has_header=False,
            with_column_names=lambda cols: [
                target_name if i == TARGET_COL else source_name if i == source_col else f"col_{i}"
                for i in range(len(cols))
            ],
            infer_schema_length=0,
            rechunk=False,
        )
        .select([source_name, target_name])
        .filter(pl.col(source_name).is_not_null() & (pl.col(source_name) != ""))
    )

    print("[*] Exploding and normalizing accessions...")
    normalized = pl.col(source_name).str.strip_chars().str.replace(r"\.\d+$", "")
    if acc_type == "pdb":
        # Only PDB carries a chain suffix. Single trailing letter only -- see docstring note 1.
        normalized = normalized.str.replace(r"_[A-Za-z]$", "")

    q = (
        q.with_columns(pl.col(source_name).str.split(";"))
        .explode(source_name)
        .with_columns(normalized)
        .filter(pl.col(source_name).is_not_null() & (pl.col(source_name) != ""))
    )

    # Keep only the accessions we actually asked about.
    targets = pl.LazyFrame({source_name: sorted(target_ids)})
    q = q.join(targets, on=source_name, how="semi")

    # Prefer reviewed (Swiss-Prot, 6-char) ACs over unreviewed (TrEMBL, 10-char) ones.
    mapping_df = (
        q.sort(pl.col(target_name).str.len_chars(), descending=False)
        .unique(subset=[source_name], keep="first")
        .collect(engine="streaming")
    )

    print("[*] Building lookup dictionary...")
    lookup_dict = dict(zip(mapping_df[source_name], mapping_df[target_name]))
    print(f"[+] Built {len(lookup_dict):,} unique mappings.")
    return lookup_dict


def read_target_ids(input_csv):
    print(f"[*] Reading target accessions from {input_csv}...")
    header_labels = {"accession", "id", "from", "entry"}
    target_ids = []
    with open(input_csv, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            token = line.strip().split(",")[0].split("\t")[0]
            if token.lower() in header_labels:
                continue
            target_ids.append(token)
    print(f"[+] Read {len(target_ids):,} accessions.")
    return target_ids


def convert_accessions(target_ids, lookup_dict, output_tsv):
    print(f"[*] Looking up {len(target_ids):,} accessions...")
    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)
    mapped_count = 0
    with open(output_tsv, "w", encoding="utf-8", newline="\n") as out_f:
        out_f.write("From\tEntry\n")
        for acc in target_ids:
            uniprot_ac = lookup_dict.get(acc)
            if uniprot_ac:
                out_f.write(f"{acc}\t{uniprot_ac}\n")
                mapped_count += 1

    pct = 100.0 * mapped_count / len(target_ids) if target_ids else 0.0
    print(f"[+] Mapped {mapped_count:,} / {len(target_ids):,} ({pct:.1f}%) -> {output_tsv}")
    if mapped_count == 0 and target_ids:
        print("[!] Zero mappings. Check --type against the input file's accession class.")
    return mapped_count


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_csv", help="accession list from split_accession_types.py")
    ap.add_argument("output_tsv", help="output TSV: From<TAB>Entry")
    ap.add_argument("--type", required=True, choices=sorted(SOURCE_COL), help="accession class")
    ap.add_argument(
        "--master",
        default="idmapping_selected.tab.gz",
        help="path to UniProt idmapping_selected.tab.gz",
    )
    args = ap.parse_args()

    if not os.path.exists(args.master):
        sys.exit(
            f"[-] Master table not found: {args.master}\n"
            "    Download from https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/idmapping/"
        )

    target_ids = read_target_ids(args.input_csv)
    if not target_ids:
        print("[!] Input is empty; writing header-only output.")
        convert_accessions([], {}, args.output_tsv)
        return

    lookup = build_mapping_dict(args.master, args.type, set(target_ids))
    convert_accessions(target_ids, lookup, args.output_tsv)


if __name__ == "__main__":
    main()
