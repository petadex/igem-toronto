#!/usr/bin/env python
"""Ground-truth set, step 2: from data/gt/pdb_accession_orfs.tsv (gt_scan_pdb_accessions.py),
pull RCSB release date / method / resolution, keep experimental structures that are
  * released after the ESMFold2 training cutoff (--cutoff, default 2021-09-30) -> "clean"
  * or before it, as length-matched leakage controls (the metric's ceiling)    -> "control"
  * X-ray or EM at <= --max-res A
then download those mmCIFs from RCSB and write the fold input (sequences from the HF orfs spine).

Public endpoints only (RCSB GraphQL + files.rcsb.org, HF parquet over HTTPS); no AWS needed.

    python gt_build_set.py                      # writes data/gt/{gt_map.tsv, gt.fa, pdb/*.cif}
"""
import argparse
import csv
import json
import os
import random
import time
import urllib.request

GQL = "https://data.rcsb.org/graphql"
QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
  rcsb_accession_info{initial_release_date} exptl{method} rcsb_entry_info{resolution_combined}}}"""


def post(url, payload, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, json.dumps(payload).encode(), {"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 ** k)


def fetch(url, dest, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r, open(dest, "wb") as fh:
                fh.write(r.read())
            return True
        except Exception:
            time.sleep(2 ** k)
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orfs", default="data/gt/pdb_accession_orfs.tsv")
    ap.add_argument("--outdir", default="data/gt")
    ap.add_argument("--cutoff", default="2021-09-30")
    ap.add_argument("--max-res", type=float, default=3.0)
    ap.add_argument("--max-len", type=int, default=2048)
    ap.add_argument("--n-control", type=int, default=40, help="pre-cutoff controls (length-matched)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.orfs), delimiter="\t"))
    for r in rows:
        r["pdb_id"] = r["accession"].split("_")[0].upper()
        r["chain"] = r["accession"].split("_")[1] if "_" in r["accession"] else ""
    ids = sorted({r["pdb_id"] for r in rows})

    info = {}
    for i in range(0, len(ids), 200):
        res = post(GQL, {"query": QUERY, "variables": {"ids": ids[i:i + 200]}})
        for e in (res.get("data") or {}).get("entries") or []:
            if not e:
                continue
            resol = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
            info[e["rcsb_id"]] = {
                "release": (e["rcsb_accession_info"]["initial_release_date"] or "")[:10],
                "method": ";".join(x["method"] for x in e.get("exptl") or []),
                "resolution": resol[0],
            }
    print(f"RCSB metadata for {len(info)}/{len(ids)} PDB entries")

    ok = []
    for r in rows:
        m = info.get(r["pdb_id"])
        if not m or int(r["length"]) > args.max_len:
            continue
        good = (("X-RAY" in m["method"] or "ELECTRON MICROSCOPY" in m["method"])
                and m["resolution"] is not None and m["resolution"] <= args.max_res)
        if not good:
            continue
        r.update(m)
        r["role"] = "clean" if m["release"] > args.cutoff else "pre_cutoff"
        ok.append(r)
    clean = [r for r in ok if r["role"] == "clean"]
    pre = [r for r in ok if r["role"] == "pre_cutoff"]
    # length-matched controls: for a random subset of clean rows, the pre-cutoff row nearest in length
    random.Random(args.seed).shuffle(clean_copy := list(clean))
    controls, used = [], set()
    for c in clean_copy:
        if len(controls) >= args.n_control or not pre:
            break
        best = min((p for p in pre if p["orf_id"] not in used),
                   key=lambda p: abs(int(p["length"]) - int(c["length"])), default=None)
        if best:
            used.add(best["orf_id"])
            best = dict(best, role="control")
            controls.append(best)
    sel = clean + controls
    print(f"selected {len(clean)} clean (post-{args.cutoff}) + {len(controls)} pre-cutoff controls "
          f"(X-ray/EM <= {args.max_res} A, <= {args.max_len} aa)")

    pdb_dir = os.path.join(args.outdir, "pdb")
    os.makedirs(pdb_dir, exist_ok=True)
    miss = 0
    for pid in sorted({r["pdb_id"] for r in sel}):
        dest = os.path.join(pdb_dir, f"{pid}.cif")
        if not os.path.exists(dest) and not fetch(f"https://files.rcsb.org/download/{pid}.cif", dest):
            miss += 1
    print(f"downloaded mmCIFs -> {pdb_dir} ({miss} failed)")

    # sequences for the selected ORFs from the HF spine (row-group pushdown on sorted orf_id)
    import duckdb
    c = duckdb.connect()
    c.sql("INSTALL httpfs; LOAD httpfs;")
    B = "https://huggingface.co/datasets/petadex/petadex-corpus-full/resolve/main/orfs/part-"
    files = [f"{B}{i}.parquet" for i in range(131)]
    want = ",".join(str(int(r["orf_id"])) for r in sel)
    seqs = dict(c.sql(f"select orf_id, rtrim(upper(sequence),'*') from read_parquet({files}) "
                      f"where orf_id in ({want})").fetchall())
    with open(os.path.join(args.outdir, "gt.fa"), "w") as fa:
        for r in sorted(sel, key=lambda r: int(r["length"])):
            s = seqs.get(int(r["orf_id"]))
            if s:
                fa.write(f">orf{r['orf_id']}\n{s}\n")
    cols = ["orf_id", "role", "pdb_id", "chain", "accession", "release", "method", "resolution",
            "length", "origin", "c90_id", "is_c90_centroid", "md5"]
    with open(os.path.join(args.outdir, "gt_map.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(sel, key=lambda r: (r["role"], int(r["length"]))))
    print(f"wrote {args.outdir}/gt.fa ({len(seqs)} seqs) and gt_map.tsv")


if __name__ == "__main__":
    main()
