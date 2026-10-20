#!/usr/bin/env python
"""Complete the 60% structure set on S3: copy-job lists for the 19,430 centroids with their own AFDB
structure and the 31 exact duplicates, a per-ORF structure index, and a sequence check of the AFDB
files against the manifest md5. Stdlib only. Used by fold60_complete_s3.sbatch.

    python fold60_complete.py lists OUTDIR
    python fold60_complete.py check-afdb STAGEDIR OUTDIR
"""
import csv
import glob
import gzip
import hashlib
import os
import sys

MAN = "data/fold60pid_input/manifest.tsv"
OVER = "data/fold60pid_input/over_cap_2500.tsv"
# 27 AFDB keys whose model is not the centroid's exact sequence (qc/afdb_inexact_refolded.tsv): folded instead
REFOLD = "data/fold60pid_input/afdb_refold.fa"
DST = "s3://petadex-protein-structures/petadex-v1.1-60pid-centroids-structures"
SOURCE = {"short_full": "esmfold2_full", "mid_fast": "esmfold2_fast", "tail_h100": "esmfold2_fast"}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def rel(tier, kind, orf):
    ext = "cif.gz" if kind == "structures" else "json"
    return f"{tier}/{kind}/{int(orf) % 1000:03d}/orf{orf}.{ext}"


def lists(out):
    os.makedirs(out, exist_ok=True)
    over = {r["orf_id"].removeprefix("orf") for r in csv.DictReader(open(OVER), delimiter="\t")}
    rows = list(csv.DictReader(open(MAN), delimiter="\t"))
    by_id = {r["orf_id"]: r for r in rows}
    refold = {l[4:].strip() for l in open(REFOLD) if l.startswith(">orf")}
    remove = open(f"{out}/afdb_remove.tsv", "w")  # inexact AFDB copies to delete from S3
    afdb = open(f"{out}/afdb_jobs.tsv", "w")      # src_s3, local_rel (uncompressed), orf_id, md5
    dups = open(f"{out}/dup_jobs.tsv", "w")       # src_s3, dst_s3
    idx = open(f"{out}/structure_index.tsv", "w")
    idx.write("orf_id\tsource\tstructure_key\tmetrics_key\tcopied_from\n")
    n = {}
    for r in rows:
        o, t, s = r["orf_id"], r["tier"], r["skip_reason"]
        if s == "afdb" and o in refold:
            remove.write(f"{DST}/{rel('afdb', 'structures', o)}\n")
            line, src = (o, SOURCE[t], rel(t, "structures", o), rel(t, "metrics", o), ""), "afdb_refolded"
        elif s == "afdb":
            k = rel("afdb", "structures", o)
            acc = r["afdb_s3_key"].rsplit("/", 1)[-1].removesuffix("-model_v4.cif")
            afdb.write(f"{r['afdb_s3_key']}\t{k.removesuffix('.gz')}\t{o}\t{r['md5']}\n")
            line, src = (o, "afdb", k, "", acc), "afdb"
        elif s == "duplicate" and by_id[r["folded_as"]]["skip_reason"] == "afdb":  # twin was not folded
            dups.write(f"{DST}/{rel('afdb', 'structures', r['folded_as'])}\t{DST}/{rel('afdb', 'structures', o)}\n")
            line, src = (o, "afdb", rel("afdb", "structures", o), "", f"orf{r['folded_as']}"), "duplicate"
        elif s == "duplicate":
            for kind in ("structures", "metrics"):
                dups.write(f"{DST}/{rel(t, kind, r['folded_as'])}\t{DST}/{rel(t, kind, o)}\n")
            line, src = (o, SOURCE[t], rel(t, "structures", o), rel(t, "metrics", o), f"orf{r['folded_as']}"), "duplicate"
        elif o in over:
            line, src = (o, "not_folded_over_2500aa", "", "", ""), "over_cap"
        else:
            line, src = (o, SOURCE[t], rel(t, "structures", o), rel(t, "metrics", o), ""), SOURCE[t]
        idx.write("\t".join(line) + "\n")
        n[src] = n.get(src, 0) + 1
    print(f"{len(rows):,} ORFs:", ", ".join(f"{k} {v:,}" for k, v in sorted(n.items())))


def cif_seq(path):
    """Sequence from the CA atoms of _atom_site (first model, chain order)."""
    cols, seq, last = [], [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.split(".", 1)[1].strip())
            elif cols and line.startswith(("ATOM", "HETATM")):
                f = dict(zip(cols, line.split()))
                if f["label_atom_id"] != "CA" or f.get("pdbx_PDB_model_num", "1") != "1":
                    continue
                key = (f["label_asym_id"], f["label_seq_id"])
                if key != last:
                    seq.append(AA3.get(f["label_comp_id"], "X"))
                    last = key
            elif cols and seq and line.startswith("#"):
                break
    return "".join(seq)


def check_afdb(stage, out):
    jobs = [l.rstrip("\n").split("\t") for l in open(f"{out}/afdb_jobs.tsv")]
    ok, bad, missing = 0, [], []
    for src, local, orf, md5 in jobs:
        p = os.path.join(stage, local)
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            missing.append(orf)
            continue
        if hashlib.md5(cif_seq(p).encode()).hexdigest() == md5:
            ok += 1
        else:
            bad.append(orf)
    with open(f"{out}/afdb_seq_check.tsv", "w") as fh:
        fh.write("orf_id\tproblem\n" + "".join(f"{o}\tmissing\n" for o in missing)
                 + "".join(f"{o}\tsequence_md5_mismatch\n" for o in bad))
    print(f"AFDB sequence check: {ok:,}/{len(jobs):,} match the manifest md5; "
          f"{len(bad)} mismatched, {len(missing)} missing")
    return 1 if missing else 0


if __name__ == "__main__":
    if sys.argv[1] == "lists":
        lists(sys.argv[2])
    elif sys.argv[1] == "check-afdb":
        sys.exit(check_afdb(sys.argv[2], sys.argv[3]))
