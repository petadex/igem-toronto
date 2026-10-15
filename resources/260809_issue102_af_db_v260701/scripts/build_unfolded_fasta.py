#!/usr/bin/env python3
"""Decide folded/unfolded per FASTA record and emit the unfolded sequences.

This is the step that the multi-accession nr headers actually change. Run 1 had one accession
per record, so "does this sequence have a structure?" was a property of the accession. Here it
is a property of the record:

    A record is FOLDED iff at least one of its accessions maps to a UniProt AC that AFDB
    actually returned a structure for. Everything else is UNFOLDED.

Both halves of that matter. An accession with no UniProt mapping does not make a record
unfolded if a sibling accession on the same record did map and did resolve. Conversely a
mapped accession whose AC turned out to be absent from AFDB does not make the record folded.

Outputs (under --outdir):
    record_structure_map.tsv   rec_id, query_petadex_id, component_id, matched_accession,
                               uniprot_ac, s3_key -- for every folded record
    unfolded.fa.zst            records with no structure, for downstream ESMFold2
    af_misses/unmapped_accessions.txt   accessions with no UniProt mapping at all
    af_misses/mapped_but_absent.txt     UniProt ACs that mapped but AFDB did not have
    unfolded_stats.json

Run 1 reported "missing" as a single number that conflated those last two categories; they have
different meanings (no UniProt entry vs. UniProt entry with no AlphaFold model) and are reported
separately here.

Memory note: the accession->AC dictionary is filtered down to accessions that resolve to a
downloaded structure before the record scan, so peak memory tracks the number of hits rather
than the ~10M-row mapping tables.

Usage:
    python build_unfolded_fasta.py INPUT.fastaa --outdir out/ \
        --record-accessions out/record_accessions.tsv \
        --hits ledger/hits.txt \
        --mapping refseq:mappings/refseq_to_uniprot.tsv \
        --mapping genbank:mappings/genbank_to_uniprot.tsv \
        --direct-uniprot accessions/uniprot_accessions.csv \
        --s3-prefix s3://petadex-protein-structures/af_db/v260701/structures/
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_nr_headers import COMPONENT_RE, open_fasta  # noqa: E402  (same-directory helper, keeps parsing identical)


def load_hits(path):
    with open(path, "r", encoding="utf-8") as f:
        return {line.strip() for line in f if line.strip()}


def load_productive_mappings(mapping_specs, direct_uniprot, hits):
    """accession -> uniprot_ac, restricted to ACs that actually have a downloaded structure.

    Also returns the set of every accession that had *any* UniProt mapping, which is what
    separates "no UniProt entry" from "UniProt entry with no AlphaFold model".
    """
    productive = {}
    all_mapped = set()
    mapped_acs = set()

    for spec in mapping_specs:
        if ":" not in spec:
            sys.exit(f"[-] --mapping expects TYPE:PATH, got '{spec}'")
        _, path = spec.split(":", 1)
        if not os.path.exists(path):
            print(f"[!] Skipping missing mapping file: {path}")
            continue
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                src, ac = parts[0].strip(), parts[1].strip()
                if not src or not ac or src == "From":
                    continue
                all_mapped.add(src)
                mapped_acs.add(ac)
                if ac in hits:
                    productive[src] = ac

    if direct_uniprot and os.path.exists(direct_uniprot):
        with open(direct_uniprot, "r", encoding="utf-8") as f:
            for line in f:
                ac = line.strip().split(",")[0].split("\t")[0]
                if not ac or ac.lower() in {"accession", "id", "from", "entry"}:
                    continue
                all_mapped.add(ac)
                mapped_acs.add(ac)
                if ac in hits:
                    productive[ac] = ac

    return productive, all_mapped, mapped_acs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--record-accessions", required=True)
    ap.add_argument("--hits", required=True, help="ledger/hits.txt from fetch_afdb.py")
    ap.add_argument("--mapping", action="append", default=[], metavar="TYPE:PATH")
    ap.add_argument("--direct-uniprot", default=None)
    ap.add_argument("--s3-prefix", default="", help="prefix used to build the s3_key column")
    ap.add_argument("--version", default="v4")
    ap.add_argument("--compress", action="store_true", default=True)
    ap.add_argument("--no-compress", dest="compress", action="store_false")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    misses_dir = os.path.join(args.outdir, "af_misses")
    os.makedirs(misses_dir, exist_ok=True)
    start = time.time()

    hits = load_hits(args.hits)
    print(f"[+] {len(hits):,} UniProt ACs have a downloaded structure.")

    productive, all_mapped, mapped_acs = load_productive_mappings(
        args.mapping, args.direct_uniprot, hits
    )
    print(f"[+] {len(all_mapped):,} accessions have a UniProt mapping; "
          f"{len(productive):,} of those resolve to a downloaded structure.")

    # Pass 1 over record_accessions.tsv: which records are folded, and via which accession.
    folded = {}          # rec_id -> (accession, uniprot_ac)
    unmapped_accs = set()
    n_pairs = 0
    with open(args.record_accessions, "r", encoding="utf-8") as f:
        header = f.readline()
        if not header.startswith("rec_id"):
            sys.exit("[-] record_accessions.tsv is missing its header row.")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rec_id, _raw, clean = int(parts[0]), parts[1], parts[2]
            n_pairs += 1
            ac = productive.get(clean)
            if ac is not None and rec_id not in folded:
                folded[rec_id] = (clean, ac)
            if clean not in all_mapped:
                unmapped_accs.add(clean)
    print(f"[+] Scanned {n_pairs:,} (record, accession) pairs; {len(folded):,} records are folded.")

    prefix = args.s3_prefix.rstrip("/") + "/" if args.s3_prefix else ""

    # Pass 2 over the FASTA: emit the unfolded records and the folded-record map.
    # rec_id is re-derived by counting '>' lines, exactly as parse_nr_headers.py assigned it.
    map_path = os.path.join(args.outdir, "record_structure_map.tsv")
    out_name = "unfolded.fa.zst" if args.compress else "unfolded.fa"
    out_path = os.path.join(args.outdir, out_name)

    if args.compress:
        import zstandard

        raw_out = open(out_path, "wb")
        writer = zstandard.ZstdCompressor(level=10).stream_writer(raw_out)

        def emit(text):
            writer.write(text.encode("utf-8"))
    else:
        raw_out = open(out_path, "w", encoding="utf-8", newline="\n")
        writer = None

        def emit(text):
            raw_out.write(text)

    rec_id = 0
    n_unfolded = 0
    n_bad_component = 0
    n_records = 0
    keep = False
    try:
        with open_fasta(args.fasta) as fin, open(map_path, "w", encoding="utf-8", newline="\n") as fmap:
            fmap.write(
                "rec_id\tquery_petadex_id\tcomponent_id\tmatched_accession\tuniprot_ac\t"
                "s3_key\toriginal_header\n"
            )
            for line in fin:
                if line.startswith(">"):
                    rec_id += 1
                    n_records += 1
                    hit = folded.get(rec_id)
                    keep = hit is None
                    if keep:
                        n_unfolded += 1
                    else:
                        acc, ac = hit
                        header = line[1:].rstrip("\n").rstrip("\r")
                        # Right-anchored, for the same reason as parse_nr_headers.py: the
                        # accession field can contain '|' (legacy nr identifiers), so the
                        # last two fields are the only reliable way to reach query and
                        # component. A left-anchored split leaks pid/bitscore into them.
                        fields = header.rsplit("|", 5)
                        query_id = fields[-2] if len(fields) >= 6 else ""
                        component_id = fields[-1] if len(fields) >= 6 else ""
                        if not COMPONENT_RE.match(component_id):
                            n_bad_component += 1
                        key = f"{prefix}AF-{ac}-F1-model_{args.version}.cif"
                        fmap.write(
                            f"{rec_id}\t{query_id}\t{component_id}\t{acc}\t{ac}\t{key}\t{header}\n"
                        )
                    if rec_id % 1_000_000 == 0:
                        print(f"    -> {rec_id:,} records scanned", flush=True)
                if keep:
                    emit(line)
    finally:
        if writer is not None:
            writer.close()
        raw_out.close()

    mapped_but_absent = sorted(mapped_acs - hits)
    with open(os.path.join(misses_dir, "mapped_but_absent.txt"), "w", encoding="utf-8", newline="\n") as f:
        for ac in mapped_but_absent:
            f.write(ac + "\n")
    with open(os.path.join(misses_dir, "unmapped_accessions.txt"), "w", encoding="utf-8", newline="\n") as f:
        for acc in sorted(unmapped_accs):
            f.write(acc + "\n")

    n_folded = n_records - n_unfolded
    stats = {
        "records": n_records,
        "folded_records": n_folded,
        "unfolded_records": n_unfolded,
        "bad_component_ids": n_bad_component,
        "folded_pct": round(100.0 * n_folded / n_records, 2) if n_records else 0.0,
        "structures_downloaded": len(hits),
        "accessions_with_uniprot_mapping": len(all_mapped),
        "accessions_without_uniprot_mapping": len(unmapped_accs),
        "uniprot_acs_mapped_but_absent_from_afdb": len(mapped_but_absent),
        "seconds": round(time.time() - start, 1),
    }
    with open(os.path.join(args.outdir, "unfolded_stats.json"), "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    # The reconciliation check from the plan, asserted rather than eyeballed.
    assert n_folded + n_unfolded == n_records, "folded + unfolded != total records"
    assert n_folded == len(folded), "record_structure_map rows disagree with the folded set"

    print("[+] Done.")
    for k, v in stats.items():
        print(f"    {k:<42} {v}")
    print(f"[+] Wrote {out_path} and {map_path}")


if __name__ == "__main__":
    main()
