#!/usr/bin/env python3
"""Build the AlphaFold DB download manifest and the UniProt-AC provenance table.

Takes the per-type mapping TSVs from local_accession_mapping.py (plus any accessions that were
already UniProt ACs and needed no mapping), unions and deduplicates the UniProt ACs, and emits
one gs:// URI per AC.

Two things differ from Run 1's inline notebook cell:

  * Deduplication across source types. An nr sequence carries several accessions, and a RefSeq
    and a GenBank accession on the same record routinely resolve to the same UniProt AC. Run 1
    kept a separate manifest and output directory per source type, which would now download the
    same structure two or three times. Here the manifest is one deduplicated list, provenance
    moves into ac_to_source.tsv, and structures/ stays flat.
  * Optional existence prefilter. Run 1 requested 830,217 files and 372,043 came back missing.
    If AFDB's accession index is available, intersecting against it up front turns those into
    zero wasted requests and gives an exact expected count before the download starts.

    For v4 the index lives at gs://public-datasets-deepmind-alphafold-v4/accession_ids.csv and
    needs an authenticated gcloud (it is not anonymously readable). Fetch it once with
        gcloud storage cp gs://public-datasets-deepmind-alphafold-v4/accession_ids.csv .
    and pass --afdb-accessions. If it is unavailable, omit the flag: fetch_afdb.py runs
    `gcloud storage cp -c` and records 404s as misses, which is Run 1's behaviour.

    Do NOT substitute EBI's anonymous accession_ids.csv here -- that index is v6, and v6 removed
    39.3M entries relative to v4, so it would wrongly prune ACs that v4 does have.

Usage:
    python build_manifest.py --outdir out/ \
        --mapping refseq:map/refseq_to_uniprot.tsv \
        --mapping genbank:map/genbank_to_uniprot.tsv \
        --mapping pdb:map/pdb_to_uniprot.tsv \
        --direct-uniprot accessions/uniprot_accessions.csv \
        [--afdb-accessions accession_ids.csv]
"""

import argparse
import json
import os
import re
import sys
import time

UNIPROT_AC_RE = re.compile(r"^[A-Z0-9]{6}(?:[A-Z0-9]{4})?$")


def read_mapping(path):
    """Yield (source_accession, uniprot_ac) from a From<TAB>Entry TSV."""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, ac = parts[0].strip(), parts[1].strip()
            if not src or not ac or src == "From":
                continue
            yield src, ac


def read_direct_uniprot(path):
    """Yield accessions that were already UniProt ACs and needed no mapping."""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            token = line.strip().split(",")[0].split("\t")[0]
            if not token or token.lower() in {"accession", "id", "from", "entry"}:
                continue
            yield token


def load_afdb_index(path):
    """Load the set of UniProt ACs present in AlphaFold DB.

    accession_ids.csv rows look like:
        A0A7Y8APW1,1,271,AF-A0A7Y8APW1-F1,6
    A bare one-column list is also accepted.
    """
    print(f"[*] Loading AFDB accession index from {path} (this file is several GB)...")
    start = time.time()
    index = set()
    with open(path, "r", encoding="utf-8", errors="replace", buffering=1 << 22) as f:
        for i, line in enumerate(f, 1):
            ac = line.split(",", 1)[0].strip()
            if ac:
                index.add(ac)
            if i % 20_000_000 == 0:
                print(f"    -> {i:,} rows...", flush=True)
    print(f"[+] Index holds {len(index):,} accessions ({time.time() - start:,.0f}s).")
    return index


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--mapping",
        action="append",
        default=[],
        metavar="TYPE:PATH",
        help="repeatable, e.g. --mapping refseq:mappings/refseq_to_uniprot.tsv",
    )
    ap.add_argument("--direct-uniprot", default=None, help="uniprot_accessions.csv from the split step")
    ap.add_argument("--afdb-accessions", default=None, help="AFDB accession index for the existence prefilter")
    ap.add_argument("--version", default="v4", help="AFDB model version (default: v4, matching Run 1)")
    ap.add_argument(
        "--gcs-bucket",
        default="public-datasets-deepmind-alphafold-v4",
        help="source bucket (default: the v4 public-datasets bucket)",
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    manifest_dir = os.path.join(args.outdir, "manifests")
    mapping_dir = os.path.join(args.outdir, "mappings")
    os.makedirs(manifest_dir, exist_ok=True)
    os.makedirs(mapping_dir, exist_ok=True)

    # uniprot_ac -> list of (source_type, source_accession)
    provenance = {}
    per_type_counts = {}

    for spec in args.mapping:
        if ":" not in spec:
            sys.exit(f"[-] --mapping expects TYPE:PATH, got '{spec}'")
        acc_type, path = spec.split(":", 1)
        if not os.path.exists(path):
            print(f"[!] Skipping missing mapping file: {path}")
            continue
        n = 0
        for src, ac in read_mapping(path):
            provenance.setdefault(ac, []).append((acc_type, src))
            n += 1
        per_type_counts[acc_type] = n
        print(f"[+] {acc_type:<8} {n:>12,} mapped pairs from {path}")

    if args.direct_uniprot and os.path.exists(args.direct_uniprot):
        n = 0
        for ac in read_direct_uniprot(args.direct_uniprot):
            provenance.setdefault(ac, []).append(("uniprot", ac))
            n += 1
        per_type_counts["uniprot"] = n
        print(f"[+] {'uniprot':<8} {n:>12,} direct accessions from {args.direct_uniprot}")

    all_acs = set(provenance)
    print(f"[+] {len(all_acs):,} unique UniProt ACs before filtering.")

    malformed = {ac for ac in all_acs if not UNIPROT_AC_RE.match(ac)}
    if malformed:
        print(f"[!] Dropping {len(malformed):,} ACs that do not look like UniProt accessions, e.g. "
              f"{sorted(malformed)[:5]}")
        all_acs -= malformed

    prefiltered = None
    if args.afdb_accessions:
        index = load_afdb_index(args.afdb_accessions)
        before = len(all_acs)
        all_acs &= index
        prefiltered = before - len(all_acs)
        print(f"[+] Prefilter removed {prefiltered:,} ACs absent from AFDB; {len(all_acs):,} remain.")
    else:
        print("[!] No --afdb-accessions index supplied. The manifest will include ACs that AFDB")
        print("    may not have; fetch_afdb.py records those as misses (Run 1 behaviour).")

    manifest_path = os.path.join(manifest_dir, "afdb_manifest.txt")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as f:
        for ac in sorted(all_acs):
            f.write(f"gs://{args.gcs_bucket}/AF-{ac}-F1-model_{args.version}.cif\n")

    prov_path = os.path.join(mapping_dir, "ac_to_source.tsv")
    with open(prov_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("uniprot_ac\tsource_type\tsource_accession\n")
        for ac in sorted(all_acs):
            for acc_type, src in provenance[ac]:
                f.write(f"{ac}\t{acc_type}\t{src}\n")

    stats = {
        "version": args.version,
        "gcs_bucket": args.gcs_bucket,
        "mapped_pairs_by_type": per_type_counts,
        "unique_acs_before_filter": len(provenance),
        "malformed_acs_dropped": len(malformed),
        "prefilter_applied": bool(args.afdb_accessions),
        "prefilter_removed": prefiltered,
        "manifest_entries": len(all_acs),
    }
    with open(os.path.join(args.outdir, "manifest_stats.json"), "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[+] Wrote {len(all_acs):,} entries -> {manifest_path}")
    print(f"[+] Wrote provenance          -> {prov_path}")


if __name__ == "__main__":
    main()
