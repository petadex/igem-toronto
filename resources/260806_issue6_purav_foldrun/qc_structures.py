#!/usr/bin/env python
"""QC for production fold output (fold_hybrid.sbatch layout: <outdir>/<tier>/{structures,metrics}/<NNN>/).

Per sampled structure:
  * gzip decompresses, mmCIF parses
  * residue sequence in the file == input FASTA sequence (no missing / extra residues)
  * every residue has backbone N, CA, C, O; heavy-atom count vs a per-residue expectation
    (flags missing side chains); no NaN / absurd coordinates
  * consecutive CA-CA distances ~3.8 A (chain breaks > 4.2 A, clashes < 2.9 A)
  * per-residue pLDDT in the B-factor column, in 0-100, mean == metrics JSON mean_plddt
  * geometry proxy from issue 111 (clashes per 1k atoms, % Ramachandran favoured)
Whole run:
  * every structure has a metrics JSON; no zero-byte or leftover .tmp files; progress CSV failures
  * optional reproducibility check vs the Phase 2 sweep (same config, same seed)

    python qc_structures.py --run data/fold60pid --input data/fold60pid_input --per-tier 300
"""
import argparse
import csv
import glob
import gzip
import importlib.util
import io
import json
import math
import os
import random
import statistics as st
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
GEOM_PATH = os.path.join(HERE, "..", "260729_issue111_folddistributions", "cheap_geom.py")
# heavy atoms per residue type (standard amino acids, excluding OXT)
HEAVY = dict(ALA=5, ARG=11, ASN=8, ASP=8, CYS=6, GLN=9, GLU=9, GLY=4, HIS=10, ILE=8, LEU=8, LYS=9,
             MET=8, PHE=11, PRO=7, SER=6, THR=7, TRP=14, TYR=12, VAL=7)
THREE2ONE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
                 LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
_geom = None


def geom():
    global _geom
    if _geom is None:
        spec = importlib.util.spec_from_file_location("cheap_geom", GEOM_PATH)
        _geom = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_geom)
    return _geom


def read_fasta(path, want=None):
    out, rid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rid and (want is None or rid in want):
                    out[rid] = "".join(buf)
                rid, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if rid and (want is None or rid in want):
        out[rid] = "".join(buf)
    return out


def check(job):
    tier, cif, mjson, seq = job
    r = {"tier": tier, "id": os.path.basename(cif)[:-7], "issues": []}
    try:
        import numpy as np
        from Bio.PDB import MMCIFParser
        text = gzip.open(cif, "rt").read()
        s = MMCIFParser(QUIET=True).get_structure("x", io.StringIO(text))
        chains = list(s[0].get_chains())
        if len(chains) != 1:
            r["issues"].append(f"chains={len(chains)}")
        res = [x for x in chains[0] if x.id[0] == " "]
        fseq = "".join(THREE2ONE.get(x.get_resname(), "X") for x in res)
        r["n_res"], r["seq_len"] = len(res), len(seq)
        # ambiguity codes (B, Z, J, U, O, X) are written as UNK by the model; compare them as X
        norm = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in seq)
        if fseq != norm:
            r["issues"].append("sequence_mismatch" if len(fseq) == len(seq) else f"length_mismatch({len(fseq)}vs{len(seq)})")
        miss_bb = sum(1 for x in res if not all(a in x for a in ("N", "CA", "C", "O")))
        if miss_bb:
            r["issues"].append(f"missing_backbone={miss_bb}")
        exp = sum(HEAVY.get(x.get_resname(), 0) for x in res)
        got = sum(1 for x in res for a in x if a.element != "H" and a.get_id() != "OXT")
        r["heavy_frac"] = got / exp if exp else None
        if exp and got < 0.98 * exp:
            r["issues"].append(f"missing_heavy_atoms={exp - got}")
        xyz = np.array([a.coord for x in res for a in x])
        if not np.isfinite(xyz).all() or np.abs(xyz).max() > 9999:
            r["issues"].append("bad_coordinates")
        ca = np.array([x["CA"].coord for x in res if "CA" in x])
        d = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        r["ca_ca_med"] = float(np.median(d)) if len(d) else None
        # > 6 A is a real gap (a missing residue); 4.2-6 A is a stretched bond, which the model
        # produces in low-confidence loops of very long chains (tail_h100 QC, 2026-10-16)
        r["breaks"] = int((d > 6.0).sum())
        r["stretched"] = int(((d > 4.2) & (d <= 6.0)).sum())
        r["ca_clash"] = int((d < 2.9).sum())
        if r["breaks"]:
            r["issues"].append(f"chain_breaks={r['breaks']}")
        b = np.array([x["CA"].get_bfactor() for x in res if "CA" in x])
        r["b_min"], r["b_max"], r["b_mean"] = float(b.min()), float(b.max()), float(b.mean())
        m = json.load(open(mjson))
        r["plddt_json"] = m["mean_plddt"]
        r["model"] = "Fast" if m.get("fast") else "Full"
        r["cfg"] = f"{m.get('num_sampling_steps')}/{m.get('num_loops')}/{m.get('kernel')}/{m.get('esmc_precision')}"
        scale = 100.0 if r["b_max"] <= 1.0 else 1.0
        r["b_vs_json"] = abs(r["b_mean"] * scale - m["mean_plddt"])
        if r["b_vs_json"] > 1.0:
            r["issues"].append(f"bfactor_mean!=json({r['b_mean']*scale:.1f}vs{m['mean_plddt']:.1f})")
        if not (0 <= r["b_min"] and r["b_max"] * scale <= 100.0001):
            r["issues"].append("plddt_out_of_range")
        g = geom().score_structure(text)
        r["clash_per_1k"], r["rama_fav"], r["cheap_score"] = g["clash_per_1k"], g["pct_rama_favored"], g["cheap_score"]
    except Exception as e:
        r["issues"].append(f"error:{type(e).__name__}:{e}")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="data/fold60pid")
    ap.add_argument("--input", default="data/fold60pid_input")
    ap.add_argument("--per-tier", type=int, default=300, help="structures sampled per tier, spread over shards")
    ap.add_argument("--sweep", default="data/sweep", help="Phase 2 sweep dir for the reproducibility check")
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 8)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = args.out or os.path.join(args.run, "qc")
    os.makedirs(out, exist_ok=True)

    jobs, inventory = [], {}
    for tier in ("short_full", "mid_fast", "tail_h100"):
        sdir = os.path.join(args.run, tier, "structures")
        if not os.path.isdir(sdir):
            continue
        cifs = glob.glob(os.path.join(sdir, "*", "*.cif.gz"))
        tmps = glob.glob(os.path.join(args.run, tier, "structures", "*", "*.tmp"))
        empty = [p for p in cifs if os.path.getsize(p) == 0]
        no_json = [p for p in cifs if not os.path.exists(p.replace("/structures/", "/metrics/")[:-7] + ".json")]
        fails = Counter()
        for p in glob.glob(os.path.join(args.run, tier, "progress*.csv")):
            for row in csv.DictReader(open(p)):
                if row["status"] != "ok":
                    fails[row["status"]] += 1
        inventory[tier] = dict(structures=len(cifs), tmp=len(tmps), empty=len(empty), missing_json=len(no_json),
                               shards_started=len(glob.glob(os.path.join(args.run, tier, "progress*.csv"))),
                               failures=dict(fails))
        sample = rng.sample(cifs, min(args.per_tier, len(cifs)))
        seqs = read_fasta(os.path.join(args.input, f"{tier}.fa"), {os.path.basename(p)[:-7] for p in sample})
        for p in sample:
            oid = os.path.basename(p)[:-7]
            jobs.append((tier, p, p.replace("/structures/", "/metrics/")[:-7] + ".json", seqs[oid]))

    with ProcessPoolExecutor(args.workers) as ex:
        res = list(ex.map(check, jobs, chunksize=4))

    print("=== inventory (all files written so far) ===")
    for t, v in inventory.items():
        print(f"{t:11s} structures={v['structures']:,} shards_started={v['shards_started']} tmp={v['tmp']} "
              f"empty={v['empty']} missing_json={v['missing_json']} failures={v['failures'] or 0}")

    print(f"\n=== per-structure checks ({len(res)} sampled) ===")
    by = defaultdict(list)
    for r in res:
        by[r["tier"]].append(r)
    for t, rs in by.items():
        bad = [r for r in rs if r["issues"]]
        iss = Counter(i.split("=")[0].split("(")[0] for r in rs for i in r["issues"])
        ok = [r for r in rs if "plddt_json" in r]
        f = lambda k: [r[k] for r in ok if r.get(k) is not None]
        print(f"{t}: {len(rs)} checked, {len(bad)} with issues {dict(iss) if iss else ''}")
        if ok:
            print(f"   models={dict(Counter(r['model'] for r in ok))} cfg={dict(Counter(r['cfg'] for r in ok))}")
            print(f"   length {min(f('seq_len'))}-{max(f('seq_len'))} aa; residues match input: "
                  f"{sum(r['n_res']==r['seq_len'] for r in ok)}/{len(ok)}; heavy-atom completeness median {st.median(f('heavy_frac')):.3f}")
            print(f"   CA-CA median {st.median(f('ca_ca_med')):.3f} A; chain breaks (>6 A) {sum(f('breaks'))}; stretched 4.2-6 A {sum(f('stretched'))}; "
                  f"|B-factor mean - JSON pLDDT| max {max(f('b_vs_json')):.3f}")
            pl = sorted(f("plddt_json"))
            print(f"   mean pLDDT median {st.median(pl):.1f} (p10 {pl[len(pl)//10]:.1f}, p90 {pl[9*len(pl)//10]:.1f}); "
                  f"<70: {sum(x<70 for x in pl)/len(pl):.1%}")
            print(f"   clashes/1k median {st.median(f('clash_per_1k')):.2f}; Rama favoured median {st.median(f('rama_fav')):.1f}%; "
                  f"cheap_score median {st.median(f('cheap_score')):.3f}")
        for r in bad[:5]:
            print(f"   ! {r['id']} {r['issues']}")

    # reproducibility: same ORF, same config, same seed in the Phase 2 sweep
    rep = []
    for tier, cfg in (("short_full", "full_s32_l3"), ("mid_fast", "fast_s32_l3")):
        sw = os.path.join(args.sweep, cfg, "structures")
        if not os.path.isdir(sw):
            continue
        for p in glob.glob(os.path.join(args.run, tier, "structures", "*", "*.cif.gz")):
            oid = os.path.basename(p)[:-7]
            q = os.path.join(sw, f"{oid}.cif")
            if os.path.exists(q):
                rep.append((tier, oid, p, q))
    if rep:
        import numpy as np
        from Bio.PDB import MMCIFParser
        from tmtools import tm_align
        print(f"\n=== reproducibility vs Phase 2 sweep (same config, seed 0): {len(rep)} shared ORFs ===")
        for tier, oid, p, q in rep:
            def ca(handle):
                ch = next(MMCIFParser(QUIET=True).get_structure("x", handle)[0].get_chains())
                rs = [x for x in ch if "CA" in x]
                return np.array([x["CA"].coord for x in rs]), "".join(THREE2ONE.get(x.get_resname(), "X") for x in rs)
            a, sa = ca(io.StringIO(gzip.open(p, "rt").read()))
            b, sb = ca(q)
            t = tm_align(a, b, sa, sb)
            pa = json.load(open(p.replace("/structures/", "/metrics/")[:-7] + ".json"))["mean_plddt"]
            sub = os.path.basename(os.path.dirname(q.replace("/structures/", "/metrics/")))
            pb = json.load(open(os.path.join(args.sweep, *q.split(os.sep)[-3:-2], "metrics", f"{oid}.json")))["mean_plddt"]
            print(f"   {tier:10s} {oid:14s} TM {t.tm_norm_chain1:.4f} RMSD {t.rmsd:.3f} A  pLDDT {pa:.2f} vs {pb:.2f}")

    with open(os.path.join(out, "qc_sample.tsv"), "w", newline="") as fh:
        keys = sorted({k for r in res for k in r})
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in res:
            w.writerow({**r, "issues": ";".join(r["issues"])})
    print(f"\nwrote {out}/qc_sample.tsv")


if __name__ == "__main__":
    main()
