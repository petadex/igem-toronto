#!/usr/bin/env python
"""Score sweep configs against experimental PDB structures (data/gt, built by gt_build_set.py).

For each config dir under --sweep-dir and each ORF in gt_map.tsv:
  * TM-score normalised by the EXPERIMENTAL chain (how well the prediction reproduces the real
    fold; TM-align is sequence-independent, so the full-ORF prediction vs a partial PDB chain is fine)
  * CA-lDDT over the TM-align residue correspondence (local accuracy; does not saturate the way
    global TM did in issue 110), inclusion radius 15 A, thresholds 0.5/1/2/4 A
  * RMSD over aligned residues
Summaries are per role (clean = released after the ESMFold2 training cutoff, control = before it)
and cluster-weighted (each 90% cluster counts once), since many clean entries are mutants of one
enzyme.

    python score_gt.py --sweep-dir data/sweep_gt --gt data/gt --out data/sweep_gt_scores
"""
import argparse
import csv
import glob
import os
import statistics as st
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


def chain_ca(path, chain=None):
    """CA coords + sequence for one chain (auth chain id); first protein chain if not given/found."""
    import numpy as np
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import three_to_index, index_to_one, is_aa
    model = MMCIFParser(QUIET=True).get_structure("x", path)[0]
    chains = {c.id: c for c in model}
    ch = chains.get(chain) if chain else None
    if ch is None:
        ch = next(c for c in model if any(is_aa(r, standard=True) for r in c))
    xyz, seq = [], []
    for res in ch:
        if res.id[0] != " " or "CA" not in res:
            continue
        try:
            seq.append(index_to_one(three_to_index(res.get_resname())))
        except Exception:
            continue
        xyz.append(res["CA"].coord)
    return np.asarray(xyz, float), "".join(seq)


def lddt_ca(pred, ref, pairs, r0=15.0):
    """CA-lDDT of pred vs ref over aligned residue index pairs [(i_pred, j_ref), ...]."""
    import numpy as np
    if len(pairs) < 3:
        return None
    P = pred[[p for p, _ in pairs]]
    R = ref[[r for _, r in pairs]]
    dR = np.linalg.norm(R[:, None] - R[None], axis=-1)
    dP = np.linalg.norm(P[:, None] - P[None], axis=-1)
    mask = (dR < r0) & ~np.eye(len(pairs), dtype=bool)
    if not mask.any():
        return None
    diff = np.abs(dR - dP)[mask]
    return float(np.mean([(diff < t).mean() for t in (0.5, 1, 2, 4)]))


def align_pairs(seqxA, seqyA):
    """Residue index pairs from TM-align's gapped alignment strings."""
    i = j = 0
    out = []
    for a, b in zip(seqxA, seqyA):
        if a != "-" and b != "-":
            out.append((i, j))
        i += a != "-"
        j += b != "-"
    return out


def score(job):
    from tmtools import tm_align
    cfg, oid, pred_path, pdb_path, chain = job
    try:
        P, sp = chain_ca(pred_path)
        R, sr = chain_ca(pdb_path, chain)
        res = tm_align(P, R, sp, sr)
        pairs = align_pairs(res.seqxA, res.seqyA)
        # lDDT over pairs that TM-align actually superposed (distance cut 5 A after alignment)
        import numpy as np
        Pt = P @ res.u.T + res.t
        close = [(i, j) for i, j in pairs if np.linalg.norm(Pt[i] - R[j]) < 5.0]
        return cfg, oid, round(res.tm_norm_chain2, 4), round(res.rmsd, 3), \
            (round(lddt_ca(P, R, pairs), 4) if pairs else None), len(close) / max(1, len(R)), len(R)
    except Exception as e:
        return cfg, oid, None, None, None, None, f"err:{type(e).__name__}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-dir", default="data/sweep_gt")
    ap.add_argument("--gt", default="data/gt")
    ap.add_argument("--out", default="data/sweep_gt_scores")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 8)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    gt = {f"orf{r['orf_id']}": r for r in csv.DictReader(open(os.path.join(args.gt, "gt_map.tsv")), delimiter="\t")}
    cfgs = sorted(os.path.basename(d.rstrip("/")) for d in glob.glob(os.path.join(args.sweep_dir, "*/"))
                  if os.path.isdir(os.path.join(d, "structures")))
    jobs = []
    for c in cfgs:
        for p in glob.glob(os.path.join(args.sweep_dir, c, "structures", "*.cif")):
            oid = os.path.basename(p)[:-4]
            if oid in gt:
                g = gt[oid]
                jobs.append((c, oid, p, os.path.join(args.gt, "pdb", f"{g['pdb_id']}.cif"), g["chain"] or None))
    with ProcessPoolExecutor(args.workers) as ex:
        res = list(ex.map(score, jobs, chunksize=8))

    with open(os.path.join(args.out, "per_orf.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["config", "id", "role", "pdb_id", "c90_id", "length", "tm_exp", "rmsd", "lddt_ca",
                    "frac_exp_within5A", "exp_len"])
        for c, oid, tm, rmsd, lddt, cov, n in res:
            g = gt[oid]
            w.writerow([c, oid, g["role"], g["pdb_id"], g["c90_id"], g["length"], tm, rmsd, lddt,
                        None if cov is None else round(cov, 3), n])

    by = defaultdict(list)
    for c, oid, tm, rmsd, lddt, cov, n in res:
        if tm is not None:
            by[(c, gt[oid]["role"])].append((gt[oid]["c90_id"], tm, lddt))
    print("config\trole\tn\tn_clusters\ttm_med\ttm_p10\tlddt_med\tlddt_p10\ttm_med_clusterw\tlddt_med_clusterw")
    rows = []
    for (c, role), xs in sorted(by.items()):
        tm = sorted(x[1] for x in xs)
        ld = sorted(x[2] for x in xs if x[2] is not None)
        cl = defaultdict(list)
        for cid, t, l in xs:
            cl[cid].append((t, l))
        tm_cw = st.median(st.mean(t for t, _ in v) for v in cl.values())
        ld_cw = st.median(st.mean(l for _, l in v if l is not None) for v in cl.values() if any(l is not None for _, l in v))
        row = [c, role, len(xs), len(cl), st.median(tm), tm[len(tm) // 10], st.median(ld), ld[len(ld) // 10], tm_cw, ld_cw]
        rows.append(row)
        print("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in row))
    with open(os.path.join(args.out, "summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["config", "role", "n", "n_clusters", "tm_med", "tm_p10", "lddt_med", "lddt_p10",
                    "tm_med_clusterw", "lddt_med_clusterw"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
