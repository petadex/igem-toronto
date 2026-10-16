#!/usr/bin/env python
"""Pairwise agreement between two fold.py output dirs on their shared ORFs.

For each ORF folded in both runs: TM-score between the two predicted structures (tmtools,
CA atoms, normalised by length), RMSD, and the deltas in mean pLDDT / pTM, plus the speed
and peak-VRAM ratio. Used to check that a cheaper config (new backend, kernel, precision,
fewer steps/loops) reproduces a reference config before adopting it.

    python compare_runs.py data/sample_fast data/stack_native_cueq --out cmp.tsv
"""
import argparse
import csv
import glob
import json
import os
import statistics as st


def load_metrics(d):
    out = {}
    for p in glob.glob(os.path.join(d, "metrics", "*.json")):
        with open(p) as fh:
            r = json.load(fh)
        if r.get("status", "ok") == "ok":
            out[r["id"]] = r
    return out


def ca_coords(cif_path):
    """CA coordinates + one-letter sequence of the first protein chain in an mmCIF."""
    import numpy as np
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import three_to_index, index_to_one

    s = MMCIFParser(QUIET=True).get_structure("x", cif_path)
    xyz, seq = [], []
    for res in next(s[0].get_chains()):
        if "CA" in res:
            xyz.append(res["CA"].coord)
            try:
                seq.append(index_to_one(three_to_index(res.get_resname())))
            except Exception:
                seq.append("X")
    return np.asarray(xyz, dtype=float), "".join(seq)


def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * (len(xs) - 1) + 0.5))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref", help="reference fold.py outdir")
    ap.add_argument("test", help="candidate fold.py outdir")
    ap.add_argument("--out", help="per-ORF TSV")
    args = ap.parse_args()

    from tmtools import tm_align

    ref, test = load_metrics(args.ref), load_metrics(args.test)
    ids = sorted(set(ref) & set(test), key=lambda i: ref[i]["seq_len"])
    if not ids:
        raise SystemExit(f"no shared ok ORFs between {args.ref} and {args.test}")

    rows = []
    for i in ids:
        a, sa = ca_coords(os.path.join(args.ref, "structures", f"{i}.cif"))
        b, sb = ca_coords(os.path.join(args.test, "structures", f"{i}.cif"))
        res = tm_align(a, b, sa, sb)
        r, t = ref[i], test[i]
        rows.append({
            "id": i, "seq_len": r["seq_len"],
            "tm": round(res.tm_norm_chain1, 4), "rmsd": round(res.rmsd, 3),
            "plddt_ref": r["mean_plddt"], "plddt_test": t["mean_plddt"],
            "d_plddt": round(t["mean_plddt"] - r["mean_plddt"], 3),
            "d_ptm": round((t["ptm"] or 0) - (r["ptm"] or 0), 4),
            "wall_ref": r["wall_s"], "wall_test": t["wall_s"],
            "speedup": round(r["wall_s"] / t["wall_s"], 3),
            "vram_ref": r.get("peak_vram_gb"), "vram_test": t.get("peak_vram_gb"),
        })

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    tm = [x["tm"] for x in rows]
    dp = [x["d_plddt"] for x in rows]
    dt = [x["d_ptm"] for x in rows]
    sp = [x["speedup"] for x in rows]
    print(f"{args.ref}  vs  {args.test}   n={len(rows)} shared ORFs "
          f"({rows[0]['seq_len']}-{rows[-1]['seq_len']}aa)")
    print(f"  TM(ref,test)  median {st.median(tm):.3f}  p10 {pct(tm, .10):.3f}  min {min(tm):.3f}"
          f"   (<0.9: {sum(x < 0.9 for x in tm)}, <0.5: {sum(x < 0.5 for x in tm)})")
    print(f"  d_pLDDT       median {st.median(dp):+.2f}  mean {st.mean(dp):+.2f}  worst {min(dp):+.2f}")
    print(f"  d_pTM         median {st.median(dt):+.4f}  mean {st.mean(dt):+.4f}  worst {min(dt):+.4f}")
    tot = sum(x["wall_ref"] for x in rows) / sum(x["wall_test"] for x in rows)
    print(f"  speedup       median {st.median(sp):.2f}x  total-time {tot:.2f}x")
    vr = [x["vram_test"] for x in rows if x["vram_test"] is not None]
    if vr:
        print(f"  peak VRAM     test max {max(vr):.1f} GB at {rows[max(range(len(rows)), key=lambda k: rows[k]['vram_test'] or 0)]['seq_len']}aa")


if __name__ == "__main__":
    main()
