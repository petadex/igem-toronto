#!/usr/bin/env python
"""Score the Phase 2 sweep (sweep.sbatch) against its reference config and the seed-noise floor.

For every config dir under --sweep-dir and every ORF it shares with the reference:
  * TM-score of the config's structure vs the reference structure (tmtools, CA, ref-normalised)
  * delta mean pLDDT / pTM vs the reference
  * geometry proxy from issue 111 (cheap_geom.score_structure: clashes/1k atoms, % Rama favoured)
  * wall time, which is power-law fitted (wall = a*L^b) and integrated over the real 90pid
    length histogram to project GPU-hours for all 18.17M centroids
The seed replicates (<cfg>_seed1/_seed2 vs the seed-0 run) give the noise floor: a config is only
"worse" than the reference by more than two seeds of the same config disagree with each other.

    python score_sweep.py --sweep-dir data/sweep --ref ref_full_s100_l20 \
        --map data/bench/bench_B_map.tsv --length-hist data/length_hist.tsv --out data/sweep_scores
"""
import argparse
import csv
import glob
import importlib.util
import json
import math
import os
import statistics as st
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
GEOM_PATH = os.path.join(HERE, "..", "260729_issue111_folddistributions", "cheap_geom.py")
LEN_BINS = [(29, 100), (100, 150), (150, 250), (250, 350), (350, 500), (500, 800), (800, 1281), (1281, 2049)]
_geom = None


def geom():
    global _geom
    if _geom is None:
        spec = importlib.util.spec_from_file_location("cheap_geom", GEOM_PATH)
        _geom = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_geom)
    return _geom


def load_metrics(d):
    out = {}
    for p in glob.glob(os.path.join(d, "metrics", "*.json")):
        with open(p) as fh:
            r = json.load(fh)
        if r.get("status", "ok") == "ok":
            out[r["id"]] = r
    return out


def ca(cif):
    import numpy as np
    xyz, seq = [], []
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    for res in next(MMCIFParser(QUIET=True).get_structure("x", cif)[0].get_chains()):
        if "CA" in res:
            xyz.append(res["CA"].coord)
            try:
                seq.append(index_to_one(three_to_index(res.get_resname())))
            except Exception:
                seq.append("X")
    return np.asarray(xyz, float), "".join(seq)


def tm_pair(job):
    """(cfg, id, cif_test, cif_ref) -> TM normalised by the reference length."""
    from tmtools import tm_align
    cfg, oid, a_path, b_path = job
    try:
        a, sa = ca(a_path)
        b, sb = ca(b_path)
        r = tm_align(a, b, sa, sb)
        return cfg, oid, round(r.tm_norm_chain2, 4), round(r.rmsd, 3)
    except Exception:
        return cfg, oid, None, None


def geom_one(job):
    cfg, oid, path = job
    try:
        with open(path) as fh:
            g = geom().score_structure(fh.read())
        return cfg, oid, g["clash_per_1k"], g["pct_rama_favored"], g["cheap_score"]
    except Exception:
        return cfg, oid, None, None, None


def len_bin(L):
    for lo, hi in LEN_BINS:
        if lo <= L < hi:
            return f"{lo}-{hi - 1}"
    return "other"


def fit_project(rows, hist, max_len):
    """Power-law fit of wall time vs length, integrated over the length histogram -> GPU-hours."""
    import numpy as np
    pts = [(r["seq_len"], r["wall_s"]) for r in rows if r.get("wall_s")]
    if len(pts) < 5:
        return None, None, None
    L = np.array([p[0] for p in pts], float)
    W = np.array([p[1] for p in pts], float)
    b, loga = np.polyfit(np.log(L), np.log(W), 1)
    a = math.exp(loga)
    secs = sum(c * a * l ** b for l, c in hist if l <= max_len)
    return secs / 3600.0, a, b


def med(xs):
    xs = [x for x in xs if x is not None]
    return st.median(xs) if xs else None


def frac(xs, f):
    xs = [x for x in xs if x is not None]
    return sum(1 for x in xs if f(x)) / len(xs) if xs else None


def fmt(x, nd=3):
    return "" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-dir", default="data/sweep")
    ap.add_argument("--ref", default="ref_full_s100_l20")
    ap.add_argument("--map", default="data/bench/bench_B_map.tsv")
    ap.add_argument("--length-hist", default="data/length_hist.tsv")
    ap.add_argument("--max-len", type=int, default=2048, help="projection cap (longer seqs excluded)")
    ap.add_argument("--out", default="data/sweep_scores")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 8)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    meta = {}
    with open(args.map) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[f"orf{r['centroid_orf_id']}"] = r
    hist = []
    with open(args.length_hist) as fh:
        for line in fh:
            l, c = line.split()
            hist.append((int(l), int(c)))

    cfgs = sorted(os.path.basename(d.rstrip("/")) for d in glob.glob(os.path.join(args.sweep_dir, "*/"))
                  if os.path.isdir(os.path.join(d, "metrics")))
    if args.ref not in cfgs:
        raise SystemExit(f"reference {args.ref} not found in {args.sweep_dir} (have {cfgs})")
    M = {c: load_metrics(os.path.join(args.sweep_dir, c)) for c in cfgs}
    ref = M[args.ref]
    sd = lambda c, i: os.path.join(args.sweep_dir, c, "structures", f"{i}.cif")

    # comparison pairs: every config vs the reference, plus seed replicates vs their seed-0 run
    pairs = [(c, args.ref) for c in cfgs if c != args.ref]
    for c in cfgs:
        for s in ("_seed1", "_seed2"):
            if c.endswith(s) and c[: -len(s)] in M:
                pairs.append((c, c[: -len(s)]))
        if c.endswith("_seed2") and c.replace("_seed2", "_seed1") in M:
            pairs.append((c, c.replace("_seed2", "_seed1")))

    tm_jobs = []
    for c, r in pairs:
        for i in sorted(set(M[c]) & set(M[r])):
            tm_jobs.append((f"{c}|{r}", i, sd(c, i), sd(r, i)))
    geom_jobs = [(c, i, sd(c, i)) for c in cfgs for i in M[c]]
    with ProcessPoolExecutor(args.workers) as ex:
        tms = list(ex.map(tm_pair, tm_jobs, chunksize=16))
        geoms = list(ex.map(geom_one, geom_jobs, chunksize=16))
    TM = {}
    for key, i, tm, rmsd in tms:
        TM.setdefault(key, {})[i] = (tm, rmsd)
    G = {(c, i): (cl, rf, cs) for c, i, cl, rf, cs in geoms}

    # per-ORF long table
    with open(os.path.join(args.out, "per_orf.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["config", "vs", "id", "seq_len", "len_bin", "origin_class", "arm", "tm", "rmsd",
                    "plddt", "plddt_vs", "ptm", "ptm_vs", "wall_s", "peak_vram_gb",
                    "clash_per_1k", "pct_rama_favored", "cheap_score"])
        for key, d in TM.items():
            c, r = key.split("|")
            for i, (tm, rmsd) in d.items():
                m, mr, mm = M[c][i], M[r][i], meta.get(i, {})
                w.writerow([c, r, i, m["seq_len"], len_bin(m["seq_len"]), mm.get("origin_class", ""),
                            mm.get("arm", ""), tm, rmsd, m["mean_plddt"], mr["mean_plddt"], m["ptm"],
                            mr["ptm"], m["wall_s"], m.get("peak_vram_gb"), *G.get((c, i), (None,) * 3)])

    # summary per comparison
    summ = []
    for key, d in TM.items():
        c, r = key.split("|")
        ids = list(d)
        tm = [d[i][0] for i in ids]
        dp = [M[c][i]["mean_plddt"] - M[r][i]["mean_plddt"] for i in ids]
        conf = [i for i in ids if M[r][i]["mean_plddt"] >= 70]
        gpu_h, a, b = fit_project(list(M[c].values()), hist, args.max_len)
        row = {
            "config": c, "vs": r, "n": len(ids),
            "gpu_h_18M": round(gpu_h) if gpu_h else None,
            "fit_b": round(b, 3) if b else None,
            "tm_med": med(tm), "tm_p10": sorted(x for x in tm if x is not None)[len(tm) // 10] if tm else None,
            "frac_tm_lt_0.5": frac(tm, lambda x: x < 0.5),
            "tm_med_ref_conf70": med([d[i][0] for i in conf]),
            "frac_tm_lt_0.5_ref_conf70": frac([d[i][0] for i in conf], lambda x: x < 0.5),
            "d_plddt_med": med(dp), "d_plddt_mean": st.mean(dp) if dp else None,
            "frac_plddt_lt70": frac([M[c][i]["mean_plddt"] for i in ids], lambda x: x < 70),
            # confidence loss: the comparison run was confident (>=70) and this run fell >15 pLDDT.
            # Low-confidence TM disagreement is seed noise; a reproducible collapse is not.
            "frac_conf_loss": frac([M[c][i]["mean_plddt"] - M[r][i]["mean_plddt"] for i in conf],
                                   lambda x: x < -15),
            "n_conf_loss": sum(1 for i in conf if M[c][i]["mean_plddt"] - M[r][i]["mean_plddt"] < -15),
            "cheap_score_med": med([G.get((c, i), (None,) * 3)[2] for i in ids]),
            "clash_per_1k_med": med([G.get((c, i), (None,) * 3)[0] for i in ids]),
        }
        for lo, hi in LEN_BINS:
            sel = [d[i][0] for i in ids if lo <= M[c][i]["seq_len"] < hi]
            row[f"tm_med_{lo}-{hi - 1}"] = med(sel)
        summ.append(row)

    ref_rows = list(ref.values())
    gpu_h, _, b = fit_project(ref_rows, hist, args.max_len)
    summ.append({"config": args.ref, "vs": "-", "n": len(ref_rows), "gpu_h_18M": round(gpu_h) if gpu_h else None,
                 "fit_b": round(b, 3) if b else None,
                 "frac_plddt_lt70": frac([x["mean_plddt"] for x in ref_rows], lambda x: x < 70),
                 "cheap_score_med": med([G.get((args.ref, x["id"]), (None,) * 3)[2] for x in ref_rows]),
                 "clash_per_1k_med": med([G.get((args.ref, x["id"]), (None,) * 3)[0] for x in ref_rows])})
    summ.sort(key=lambda r: (r["vs"] != args.ref, r.get("gpu_h_18M") or 1e12))
    cols = sorted({k for r in summ for k in r}, key=lambda k: list(summ[0]).index(k) if k in summ[0] else 99)
    with open(os.path.join(args.out, "summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(summ)

    show = ["config", "vs", "n", "gpu_h_18M", "tm_med_ref_conf70", "frac_tm_lt_0.5_ref_conf70",
            "n_conf_loss", "frac_conf_loss", "d_plddt_med", "frac_plddt_lt70", "cheap_score_med"]
    print("\t".join(show))
    for r in summ:
        print("\t".join(fmt(r.get(k)) for k in show))
    print(f"\nwrote {args.out}/summary.tsv and per_orf.tsv")


if __name__ == "__main__":
    main()
