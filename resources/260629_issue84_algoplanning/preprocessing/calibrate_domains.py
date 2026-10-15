"""
calibrate_domains.py -- pick pae_cutoff and resolution ONCE, for the whole corpus.

Issue #84, Stage 0.  Self-contained; imports nothing from this repo.

WHY A SWEEP AND NOT A DEFAULT
-----------------------------
The usual `pae_cutoff = 5.0, resolution = 1.0` is calibrated on AF2 predictions
of general proteins.  Our PAE comes from ESMFold2 and our proteins are a narrow
~250-400 aa band of one enzyme family, so neither number transfers by assumption.
These are GLOBAL parameters -- one setting for every cluster.  If they had to be
retuned per cluster we would be back in the hand-tuned-threshold trap that this
whole route exists to escape.

THE SCORING TRICK: WE ALREADY HAVE A LABEL
------------------------------------------
Every ORF carries the PETadex HMM hit region (`core_range`, e.g. 65-261).  It is
not a domain boundary -- the HMM under-covers, which is the entire reason Stage 0
exists -- but it is a reliable ANCHOR: the true enzyme domain must contain it.
That makes two failure modes detectable with no ground truth at all:

    split   the anchor's residues land in >1 cluster        -> over-segmented
    whole   the anchor's cluster covers ~the entire ORF     -> under-segmented
    ok      neither                                          -> plausible

Pick the widest (cutoff, resolution) region where `ok` dominates, then confirm the
choice sits on a PLATEAU rather than a spike (--report-stability).  ORFs with a
PDB in `paramfold_orfs_metadata.tsv` additionally allow a real residue-error
metric later; this script only does the label-free part.

Note that most ORFs are single-domain, where any setting returns one cluster and
the boundary comes from pLDDT trimming rather than PAE clustering.  The parameters
only bite on the ~10-25% carrying fusion partners, so `--min-excess` restricts
scoring to ORFs long enough for the choice to matter.

ENGINES
-------
    afragmenter  preferred; pip install afragmenter        (Leiden on a PAE network)
    leiden       inline igraph + leidenalg, Croll's recipe (same idea, no wrapper)
    components   numpy only, connected components at the cutoff.  DEBUG ONLY --
                 it is not Leiden and must not be used to choose production values.

USAGE (on the VM -- the PAE matrices live in S3 and are too big to be worth pulling)
-----
    pip install afragmenter
    python calibrate_domains.py \\
        --s3 s3://petadex-protein-structures/esmfold2_paramsweep/s68_l10/ \\
        --anchors ../../260706_issue6_foldpredictions/paramfold_orfs_metadata.tsv \\
        --pae-cutoff 3,5,8,10 --resolution 0.3,0.5,1.0,2.0 \\
        -o sweep.tsv --report-stability
"""

import argparse
import json
import os
import sys

import numpy as np


# --------------------------------------------------------------------------- #
# inputs
# --------------------------------------------------------------------------- #

def read_anchors(path):
    """paramfold_orfs_metadata.tsv -> {orf_id: (start, end, orf_len)}, 1-based
    inclusive.  `core_range` is the PETadex HMM hit region."""
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
            if "-" not in rng:
                continue
            a, _, b = rng.partition("-")
            if not (a.isdigit() and b.isdigit()):
                continue
            olen = 0
            if "orf_len" in idx and f[idx["orf_len"]].isdigit():
                olen = int(f[idx["orf_len"]])
            anchors[f[idx["orfid"]].strip()] = (int(a), int(b), olen)
    return anchors


def orf_id_from_name(name):
    """Fold outputs are named from the FASTA header with separators flattened, so the
    same ORF can appear as any of:
        orf4981589                              (bare-id input, e.g. paramfold)
        1219585044__SRR6391592_110893_1_709_1   (full PETadex header, pipes -> underscores)
    All reduce to the ORF id, which is what anchors are keyed by."""
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
    """Whole-object GETs -- unavoidable here.  The 4 KB Range-GET trick used by
    conf_metrics.py works precisely because it SKIPS pae, which is what we need."""
    import boto3
    assert uri.startswith("s3://")
    bucket, _, prefix = uri[5:].partition("/")
    if prefix and not prefix.endswith("/"):
        prefix += "/"
    cl = boto3.client("s3")
    pages = cl.get_paginator("list_objects_v2").paginate(
        Bucket=bucket, Prefix=prefix + "metrics/")
    for page in pages:
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not key.endswith(".json"):
                continue
            body = cl.get_object(Bucket=bucket, Key=key)["Body"].read()
            yield os.path.basename(key)[:-5], json.loads(body)


def extract(rec):
    """ESMFold2 metrics JSON -> (symmetrized PAE, per-residue pLDDT).
    Ours nests these under `confidence`; AF2 puts PAE at the top level as
    `predicted_aligned_error`.  Both accepted."""
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
    pae = (pae + pae.T) / 2.0            # ESMFold2 PAE is NOT symmetric:
    return pae, np.asarray(plddt, float)  # max |P - P.T| ~16 A on real output


# --------------------------------------------------------------------------- #
# clustering engines
# --------------------------------------------------------------------------- #

# AFragmenter's default is n_iterations=-1, meaning "iterate until convergence".
# That has no bound, so a single awkward graph can stall a whole cell with no output.
# 2 is the usual practical setting; --n-iterations restores -1 if you want it.
N_ITERATIONS = 2


def smooth(x, w):
    """Moving average that does NOT decay at the edges.  A plain
    np.convolve(mode='same') zero-pads, which drags the first and last w//2
    residues toward zero -- on a real ORF that turned genuine pLDDT ~90 termini
    into ~49 and would have masked away exactly the boundary residues we are
    trying to call.  Dividing by the count of real contributions fixes it."""
    if w <= 1:
        return x
    k = np.ones(w)
    return np.convolve(x, k, mode="same") / np.convolve(np.ones_like(x), k, mode="same")


def cluster_leiden(pae, cutoff, resolution, min_domain):
    import igraph
    import leidenalg
    n = len(pae)
    i, j = np.triu_indices(n, k=1)
    m = pae[i, j] < cutoff
    if not m.any():
        return [set(range(n))]
    w = 1.0 / np.maximum(pae[i, j][m], 0.25)      # 0.25 = the PAE bin floor
    g = igraph.Graph(n=n, edges=list(zip(i[m].tolist(), j[m].tolist())))
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        weights=w.tolist(), resolution_parameter=resolution)
    return [set(c) for c in part]


def cluster_afragmenter(pae, cutoff, resolution, min_domain):
    """Real API (verified against the installed package):
        AFragmenter(pae_matrix, threshold=2.0, sequence_file=None)
        .cluster(resolution=None, objective_function='modularity', n_iterations=-1,
                 min_size=10, attempt_merge=True, min_avg_pae=None, **kwargs)
        -> ClusteringResult.cluster_intervals == {cid: [(start, end), ...]}

    NOTE `threshold` is a CONSTRUCTOR argument, not a cluster() one -- cluster()
    forwards unknown kwargs straight to igraph's community_leiden, which rejects
    them with a bare "unexpected keyword argument".
    Intervals are 0-based INCLUSIVE, and a cluster may hold several of them
    (discontiguous domains), which is why we return residue sets.

    NOTE ALSO `min_size` is not a filter -- it feeds `attempt_merge` and RESHAPES
    the partition, non-monotonically.  Measured on a real ORF: raw clusters
    [48, 28, 19] became EMPTY at min_size=20 (all three are >= 20), and [58, 37]
    became [31] -- a size present in neither.  So we take the raw partition and
    apply our own --min-domain filter in score(), which keeps the size rule in
    one predictable place."""
    from afragmenter import AFragmenter
    res = AFragmenter(pae, threshold=cutoff).cluster(
        resolution=resolution, min_size=1, n_iterations=N_ITERATIONS)
    out = []
    for _, intervals in res.cluster_intervals.items():
        s = set()
        for a, b in intervals:
            s.update(range(int(a), int(b) + 1))
        if s:
            out.append(s)
    return out


def cluster_components(pae, cutoff, resolution, min_domain):
    """DEBUG baseline: connected components of the PAE<cutoff graph.  Ignores
    resolution entirely.  Present so the pipeline is exercisable without igraph;
    NOT a substitute for Leiden."""
    n = len(pae)
    adj = pae < cutoff
    seen, out = np.zeros(n, bool), []
    for s in range(n):
        if seen[s]:
            continue
        stack, comp = [s], set()
        seen[s] = True
        while stack:
            u = stack.pop()
            comp.add(u)
            for v in np.nonzero(adj[u] & ~seen)[0]:
                seen[v] = True
                stack.append(v)
        out.append(comp)
    return out


ENGINES = {"afragmenter": cluster_afragmenter,
           "leiden": cluster_leiden,
           "components": cluster_components}


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #

def score(domains, anchor, n_res, args):
    """-> (verdict, n_domains, anchor_dom_size, frac_of_orf, discontiguous)"""
    a, b = anchor
    anchor_res = set(range(a - 1, min(b, n_res)))
    if not anchor_res:
        return "no_anchor", len(domains), 0, 0.0, 0

    doms = [d for d in domains if len(d) >= args.min_domain]
    if not doms:
        return "no_domain", len(domains), 0, 0.0, 0

    overlaps = [(len(d & anchor_res), d) for d in doms]
    best_n, best = max(overlaps, key=lambda t: t[0])
    if best_n == 0:
        return "no_domain", len(doms), 0, 0.0, 0

    covered = best_n / len(anchor_res)
    frac = len(best) / n_res
    disc = int((max(best) - min(best) + 1) > len(best) * 1.1)

    if covered < args.anchor_covered:
        v = "split"
    elif frac >= args.whole_frac:
        v = "whole"
    else:
        v = "ok"
    return v, len(doms), len(best), frac, disc


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--metrics-dir", help="local dir of <orfid>.json")
    src.add_argument("--s3", help="s3://bucket/prefix/  (reads <prefix>/metrics/)")
    ap.add_argument("--anchors", required=True, help="paramfold_orfs_metadata.tsv")
    ap.add_argument("--pae-cutoff", default="3,5,8,10")
    ap.add_argument("--resolution", default="0.3,0.5,1.0,2.0")
    ap.add_argument("--engine", default="afragmenter", choices=sorted(ENGINES))
    ap.add_argument("--min-plddt", type=float, default=70.0)
    ap.add_argument("--plddt-window", type=int, default=15)
    ap.add_argument("--min-domain", type=int, default=30)
    ap.add_argument("--anchor-covered", type=float, default=0.9,
                    help="anchor fraction that must land in ONE cluster (0.9)")
    ap.add_argument("--whole-frac", type=float, default=0.9,
                    help="anchor cluster covering this much of the ORF = no split (0.9)")
    ap.add_argument("--min-excess", type=int, default=0,
                    help="only score ORFs at least this many aa longer than their "
                         "anchor -- the parameters only matter for those")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--max-len", type=int, default=0,
                    help="skip ORFs longer than this. Our design targets are ~250-400 aa "
                         "and resolution does not transfer across very different lengths, "
                         "so long controls are both slow and off-target for calibration.")
    ap.add_argument("--n-iterations", type=int, default=2,
                    help="Leiden iterations (2). -1 = until convergence, which is "
                         "AFragmenter's default and has no time bound.")
    ap.add_argument("-o", "--out", default="sweep.tsv")
    ap.add_argument("--report-stability", action="store_true")
    args = ap.parse_args()

    global N_ITERATIONS
    N_ITERATIONS = args.n_iterations

    cutoffs = [float(x) for x in args.pae_cutoff.split(",")]
    resols = [float(x) for x in args.resolution.split(",")]
    engine = ENGINES[args.engine]
    if args.engine == "components":
        print("WARNING: --engine components is a DEBUG baseline, not Leiden. "
              "Do not choose production parameters from it.", file=sys.stderr)

    anchors = read_anchors(args.anchors)
    it = iter_metrics_local(args.metrics_dir) if args.metrics_dir else iter_metrics_s3(args.s3)

    proteins = []
    no_anchor = no_pae = too_long = too_short_excess = 0
    for oid, rec in it:
        key = orf_id_from_name(oid)
        if key not in anchors:
            no_anchor += 1
            continue
        pae, plddt = extract(rec)
        if pae is None:
            no_pae += 1
            continue
        a, b, _ = anchors[key]
        if args.max_len and len(pae) > args.max_len:
            too_long += 1
            continue
        if args.min_excess and len(pae) - (b - a + 1) < args.min_excess:
            too_short_excess += 1
            continue
        proteins.append((key, pae, plddt, (a, b)))
        if len(proteins) % 25 == 0:
            print("  fetched %d..." % len(proteins), file=sys.stderr, flush=True)
        if args.limit and len(proteins) >= args.limit:
            break

    tally = ("kept %d  |  dropped: no_anchor=%d no_pae=%d too_long(--max-len)=%d "
             "excess_below(--min-excess)=%d"
             % (len(proteins), no_anchor, no_pae, too_long, too_short_excess))
    if not proteins:
        sys.exit("no usable proteins.\n  " + tally +
                 "\nIf the drops are all from --min-excess, this fold set has no ORFs "
                 "with that much material beyond their anchor -- it cannot calibrate a "
                 "ceiling, and you need folds of longer ORFs.")
    print(tally + "\n")

    rows, grid = [], {}
    ncells = len(cutoffs) * len(resols)
    cell = 0
    for cut in cutoffs:
        for res in resols:
            cell += 1
            print("  cell %d/%d  cutoff=%s resolution=%s" % (cell, ncells, cut, res),
                  file=sys.stderr, flush=True)
            counts = {}
            for pi, (key, pae, plddt, anchor) in enumerate(proteins, 1):
                if pi % 25 == 0:
                    print("      %d/%d proteins" % (pi, len(proteins)),
                          file=sys.stderr, flush=True)
                doms = engine(pae, cut, res, args.min_domain)
                if args.min_plddt > 0:
                    good = smooth(plddt, args.plddt_window) >= args.min_plddt
                    doms = [{r for r in d if good[r]} for d in doms]
                v, nd, sz, frac, disc = score(doms, anchor, len(pae), args)
                counts[v] = counts.get(v, 0) + 1
                rows.append([cut, res, key, nd, sz, len(pae), round(frac, 3), disc, v])
            grid[(cut, res)] = counts

    with open(args.out, "w") as fh:
        fh.write("pae_cutoff\tresolution\torf_id\tn_domains\tanchor_dom_size\t"
                 "orf_len\tfrac_of_orf\tdiscontiguous\tverdict\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    n = len(proteins)
    print("%% 'ok' by (pae_cutoff x resolution) -- pick a PLATEAU, not a peak\n")
    print("cutoff\\res  " + "".join("%8s" % r for r in resols))
    for cut in cutoffs:
        cells = ["%7.0f%%" % (100.0 * grid[(cut, r)].get("ok", 0) / n) for r in resols]
        print("%-11s " % cut + "".join("%8s" % c for c in cells))
    print("\nfailure modes at each cell (split = over-segmented, whole = under):")
    for cut in cutoffs:
        for res in resols:
            c = grid[(cut, res)]
            print("  cut %-5s res %-5s  " % (cut, res) +
                  "  ".join("%s=%d" % (k, c[k]) for k in sorted(c)))

    if args.report_stability:
        print("\nstability -- ORFs whose verdict changes between adjacent resolutions:")
        for cut in cutoffs:
            for r1, r2 in zip(resols, resols[1:]):
                d = {}
                for r in rows:
                    if r[0] == cut and r[1] in (r1, r2):
                        d.setdefault(r[2], {})[r[1]] = r[8]
                ch = sum(1 for v in d.values() if len(v) == 2 and v[r1] != v[r2])
                print("  cut %-5s %.2f -> %.2f : %d / %d changed" % (cut, r1, r2, ch, n))

    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
