"""
domain_from_pae.py -- Stage 0, step 3:  folded structures  ->  per-ORF domain boundaries.

Issue #84.  Self-contained; imports nothing from this repo.

    fetch_clusters.py -> esmfold2_local_predictor.py -> THIS -> [mafft] -> assemble_cores.py

WHAT IT DOES
------------
A thin adapter around AFragmenter.  AFragmenter clusters a PAE matrix into pseudo-rigid
domains; it does not know our metrics-JSON layout, it does not know WHICH of the domains it
returns is the enzyme, and it does not emit the format assemble_cores.py reads.  Those
three things are this script.

    metrics/<id>.json -> symmetrize PAE -> AFragmenter -> pick the domain containing the
    HMM anchor -> pLDDT-trim the ends -> boundaries.tsv (orf_id, start, end)

PARAMETERS ARE ALREADY CALIBRATED -- DO NOT RETUNE PER CLUSTER
--------------------------------------------------------------
`--threshold 2 --resolution 0.3`, established in DOMAIN_CALIBRATION.md against two
independent populations (74 single-domain ORFs for the floor, 6 genuine fusion-carrying
ORFs for the ceiling).  Both bounds agree across 0.1-0.5 and fail only at 1.0.  These are
GLOBAL -- one setting for every cluster.  Retuning them per cluster would reintroduce
exactly the hand-tuned-threshold problem this route exists to remove.

WHY AN ANCHOR IS NEEDED
-----------------------
The PETadex HMM hit region (`core_range`) is NOT a domain boundary -- it under-covers by a
median 92 aa, which is why Stage 0 exists at all.  It is used only to identify WHICH
cluster is the enzyme.  The obvious alternative, "take the largest cluster", fails exactly
where it matters: an MBP fusion is ~370 aa against a ~290 aa PETase, so the tag would win.

CONTIGUITY
----------
AFragmenter returns residue SETS and supports discontiguous domains (an inserted domain
splits its host into two sequence segments).  But we are ordering DNA, so the boundary we
emit has to be one contiguous stretch.  Ragged edges are healed by merging runs separated
by <= --max-gap residues; if a genuinely separate block remains, we keep the run holding
the most anchor residues and flag the ORF as `discontiguous` in the report so it can be
inspected rather than silently spanned.

USAGE
-----
    python domain_from_pae.py \\
        --metrics-dir out/clusters/1/fold/metrics \\
        --anchors cluster_orfids/anchors.tsv \\
        -o out/clusters/1/boundaries.tsv
"""

import argparse
import json
import os
import sys

import numpy as np


# --------------------------------------------------------------------------- #
# input
# --------------------------------------------------------------------------- #

def read_anchors(path):
    """orfid -> (start, end), 1-based inclusive, from a `core_range` column."""
    anchors = {}
    with open(path) as fh:
        cols = fh.readline().rstrip("\n\r").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        for need in ("orfid", "core_range"):
            if need not in idx:
                sys.exit("%s lacks a '%s' column (has %s)" % (path, need, cols))
        for line in fh:
            f = line.rstrip("\n\r").split("\t")
            if len(f) <= idx["core_range"]:
                continue
            rng = f[idx["core_range"]].strip()
            a, _, b = rng.partition("-")
            if a.isdigit() and b.isdigit():
                anchors[f[idx["orfid"]].strip()] = (int(a), int(b))
    return anchors


def orf_id_from_name(name):
    """Fold outputs are named from the FASTA header with separators flattened, so one ORF
    can appear as `orf4981589` or `1219585044__SRR6391592_110893_1_709_1`.  Both reduce to
    the ORF id, which is what anchors and members.tsv are keyed by."""
    n = os.path.basename(name)
    if n.endswith(".json"):
        n = n[:-5]
    if n.startswith("orf"):
        n = n[3:]
    return n.split("|")[0].split("_")[0].strip()


def iter_metrics_local(d):
    for fn in sorted(os.listdir(d)):
        if fn.endswith(".json"):
            with open(os.path.join(d, fn)) as fh:
                yield fn[:-5], json.load(fh)


def iter_metrics_s3(uri):
    """Whole-object GETs -- the 4 KB Range-GET trick used elsewhere works precisely
    because it SKIPS pae, which is the field we need."""
    import boto3
    assert uri.startswith("s3://")
    bucket, _, prefix = uri[5:].partition("/")
    if prefix and not prefix.endswith("/"):
        prefix += "/"
    cl = boto3.client("s3")
    for page in cl.get_paginator("list_objects_v2").paginate(
            Bucket=bucket, Prefix=prefix + "metrics/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith(".json"):
                body = cl.get_object(Bucket=bucket, Key=key)["Body"].read()
                yield os.path.basename(key)[:-5], json.loads(body)


def extract(rec):
    """-> (symmetrized PAE, per-residue pLDDT).  Ours nests these under `confidence`;
    AF2 puts PAE at the top level as `predicted_aligned_error`.  Both accepted."""
    conf = rec.get("confidence", rec)
    pae = conf.get("pae", rec.get("predicted_aligned_error"))
    plddt = conf.get("per_residue_plddt", rec.get("plddt"))
    if pae is None or plddt is None:
        return None, None
    pae = np.asarray(pae, dtype=float)
    if pae.ndim == 3 and pae.shape[0] == 1:
        pae = pae[0]
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        return None, None
    return (pae + pae.T) / 2.0, np.asarray(plddt, float)   # ESMFold2 PAE is not symmetric


# --------------------------------------------------------------------------- #
# clustering
# --------------------------------------------------------------------------- #

def smooth(x, w):
    """Moving average that does NOT decay at the edges.  np.convolve(mode='same')
    zero-pads, which drags the first/last w//2 residues toward zero and would mask away
    exactly the boundary residues we are trying to call."""
    if w <= 1:
        return x
    k = np.ones(w)
    return np.convolve(x, k, mode="same") / np.convolve(np.ones_like(x), k, mode="same")


def cluster_afragmenter(pae, threshold, resolution, n_iterations):
    """See DOMAIN_CALIBRATION.md for the API notes: `threshold` is a CONSTRUCTOR argument,
    `min_size` reshapes the partition rather than filtering it (so pass 1 and filter
    ourselves), and n_iterations=-1 is unbounded and can hang."""
    from afragmenter import AFragmenter
    res = AFragmenter(pae, threshold=threshold).cluster(
        resolution=resolution, min_size=1, n_iterations=n_iterations)
    out = []
    for _, intervals in res.cluster_intervals.items():
        s = set()
        for a, b in intervals:
            s.update(range(int(a), int(b) + 1))
        if s:
            out.append(s)
    return out


def cluster_leiden(pae, threshold, resolution, n_iterations):
    """Croll's pae_to_domains recipe directly, for when the AFragmenter API moves."""
    import igraph
    import leidenalg
    n = len(pae)
    i, j = np.triu_indices(n, k=1)
    m = pae[i, j] < threshold
    if not m.any():
        return [set(range(n))]
    w = 1.0 / np.maximum(pae[i, j][m], 0.25)
    g = igraph.Graph(n=n, edges=list(zip(i[m].tolist(), j[m].tolist())))
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition, weights=w.tolist(),
        resolution_parameter=resolution, n_iterations=n_iterations)
    return [set(c) for c in part]


ENGINES = {"afragmenter": cluster_afragmenter, "leiden": cluster_leiden}


# --------------------------------------------------------------------------- #
# domain selection
# --------------------------------------------------------------------------- #

def runs_of(residues, max_gap):
    """Sorted residue set -> [(start, end)] inclusive, merging runs separated by
    <= max_gap missing residues.  Clusters on real PAE are slightly ragged; without this
    a couple of dropped residues would shatter one domain into several."""
    if not residues:
        return []
    rs = sorted(residues)
    out = [[rs[0], rs[0]]]
    for r in rs[1:]:
        if r - out[-1][1] - 1 <= max_gap:
            out[-1][1] = r
        else:
            out.append([r, r])
    return [(a, b) for a, b in out]


def pick_domain(domains, anchor, n_res, args):
    """-> (start, end, info) with start/end 1-based inclusive, or (None, None, info)."""
    a, b = anchor
    anchor_res = set(range(a - 1, min(b, n_res)))
    info = {"n_domains": len(domains), "anchor_covered": 0.0,
            "frac_of_orf": 0.0, "discontiguous": 0, "note": ""}
    if not anchor_res:
        info["note"] = "anchor outside sequence"
        return None, None, info

    doms = [d for d in domains if len(d) >= args.min_domain]
    info["n_domains"] = len(doms)
    if not doms:
        info["note"] = "no domain >= --min-domain"
        return None, None, info

    best_n, best = max(((len(d & anchor_res), d) for d in doms), key=lambda t: t[0])
    if best_n == 0:
        info["note"] = "no domain overlaps the anchor"
        return None, None, info
    info["anchor_covered"] = round(best_n / len(anchor_res), 3)

    segs = runs_of(best, args.max_gap)
    if len(segs) > 1:
        info["discontiguous"] = 1
        # keep the segment holding the most anchor residues rather than spanning the gap,
        # which would swallow whatever sits between two genuinely separate blocks.
        segs = [max(segs, key=lambda s: len(set(range(s[0], s[1] + 1)) & anchor_res))]
    s, e = segs[0]
    info["frac_of_orf"] = round((e - s + 1) / n_res, 3)
    if info["anchor_covered"] < args.anchor_covered:
        info["note"] = "anchor split across clusters"
    return s + 1, e + 1, info


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--metrics-dir", help="local dir of <orfid>.json")
    src.add_argument("--s3", help="s3://bucket/prefix/  (reads <prefix>/metrics/)")
    ap.add_argument("--anchors", required=True, help="TSV with orfid + core_range")
    ap.add_argument("-o", "--out", required=True, help="boundaries.tsv for assemble_cores")
    ap.add_argument("--report", help="per-ORF diagnostics TSV (default: <out>.report.tsv)")
    ap.add_argument("--engine", default="afragmenter", choices=sorted(ENGINES))
    ap.add_argument("--threshold", type=float, default=2.0,
                    help="AFragmenter contrast threshold. CALIBRATED -- see "
                         "DOMAIN_CALIBRATION.md. Do not retune per cluster.")
    ap.add_argument("--resolution", type=float, default=0.3,
                    help="Leiden resolution. CALIBRATED -- see DOMAIN_CALIBRATION.md.")
    ap.add_argument("--n-iterations", type=int, default=2,
                    help="Leiden iterations (2). -1 = until convergence, unbounded.")
    ap.add_argument("--min-plddt", type=float, default=70.0,
                    help="trim residues below this pLDDT (0 disables)")
    ap.add_argument("--plddt-window", type=int, default=15)
    ap.add_argument("--min-domain", type=int, default=30)
    ap.add_argument("--max-gap", type=int, default=10,
                    help="merge domain segments separated by <= this many residues")
    ap.add_argument("--anchor-covered", type=float, default=0.9,
                    help="warn when less than this fraction of the anchor is in one cluster")
    args = ap.parse_args()

    anchors = read_anchors(args.anchors)
    engine = ENGINES[args.engine]
    it = iter_metrics_local(args.metrics_dir) if args.metrics_dir else iter_metrics_s3(args.s3)

    rows, bounds = [], []
    n_seen = no_anchor = no_pae = failed = 0
    for name, rec in it:
        n_seen += 1
        oid = orf_id_from_name(name)
        if oid not in anchors:
            no_anchor += 1
            rows.append([oid, "", "", 0, 0.0, 0.0, 0, "no anchor"])
            continue
        pae, plddt = extract(rec)
        if pae is None:
            no_pae += 1
            rows.append([oid, "", "", 0, 0.0, 0.0, 0, "no PAE"])
            continue

        doms = engine(pae, args.threshold, args.resolution, args.n_iterations)
        if args.min_plddt > 0:
            good = smooth(plddt, args.plddt_window) >= args.min_plddt
            doms = [{r for r in d if good[r]} for d in doms]

        s, e, info = pick_domain(doms, anchors[oid], len(pae), args)
        rows.append([oid, s if s else "", e if e else "", info["n_domains"],
                     info["anchor_covered"], info["frac_of_orf"],
                     info["discontiguous"], info["note"]])
        if s is None:
            failed += 1
        else:
            bounds.append((oid, s, e))

        if n_seen % 25 == 0:
            print("  %d processed..." % n_seen, file=sys.stderr, flush=True)

    if not bounds:
        sys.exit("no boundaries produced from %d metrics files "
                 "(no_anchor=%d no_pae=%d failed=%d)" % (n_seen, no_anchor, no_pae, failed))

    with open(args.out, "w") as fh:
        fh.write("orf_id\tstart\tend\n")
        for oid, s, e in sorted(bounds, key=lambda t: int(t[0]) if t[0].isdigit() else 0):
            fh.write("%s\t%d\t%d\n" % (oid, s, e))

    rpath = args.report or (args.out + ".report.tsv")
    with open(rpath, "w") as fh:
        fh.write("orf_id\tstart\tend\tn_domains\tanchor_covered\tfrac_of_orf\t"
                 "discontiguous\tnote\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    lens = [e - s + 1 for _, s, e in bounds]
    lens.sort()
    disc = sum(1 for r in rows if r[6] == 1)
    split = sum(1 for r in rows if r[7] == "anchor split across clusters")
    whole = sum(1 for r in rows if r[5] and r[5] >= 0.9)
    print("\n%d metrics files -> %d boundaries" % (n_seen, len(bounds)))
    print("  dropped:      no_anchor=%d  no_pae=%d  no_domain=%d"
          % (no_anchor, no_pae, failed))
    print("  domain length: min %d  median %d  max %d" % (lens[0], lens[len(lens) // 2], lens[-1]))
    print("  discontiguous: %d     anchor-split warnings: %d" % (disc, split))
    print("  spans >=90%% of the ORF (i.e. nothing trimmed): %d" % whole)
    print("\nwrote %s\nwrote %s" % (args.out, rpath))


if __name__ == "__main__":
    main()
