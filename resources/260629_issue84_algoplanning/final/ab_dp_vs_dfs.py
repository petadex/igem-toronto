"""A/B the whole designer with the stock DFS cut search vs the DP/CSP one.

cutsearch_design.py is NOT edited.  `evaluate_K` looks `place_cuts` up as a
module global at call time, so rebinding that name in this driver swaps the cut
search for one run and restores it for the other.  Everything else -- the greedy,
the caps, the seed, the random arm, the selection rule -- is identical between
arms, so any difference in the frontier is attributable to the cut search alone.

  python ab_dp_vs_dfs.py <cores.aln.fasta> [--max-library 2500] [--k-max 6]
"""

import argparse
import sys
import time

from _paths import inputs_from_argv          # noqa: F401  (puts HERE on sys.path)
import cutsearch_design as U
import dp_cutsearch as D


def dp_shim(aligned, L, K, min_block, const, chemistry, arm_codons, reserved,
            node_budget=2_000_000, n_keep=1, pool_factor=10, max_layer_cols=None):
    """place_cuts' exact call signature, backed by the DP/CSP search."""
    return D.dp_cut_search(aligned, L, K, min_block, const, chemistry,
                           arm_codons, reserved, n_keep=n_keep,
                           pool_factor=pool_factor,
                           max_layer_cols=max_layer_cols)


def run_arm(name, aligned, weights, L, const, args, reserved, overhead):
    stock = U.place_cuts
    if name == "dp":
        U.place_cuts = dp_shim
    try:
        results = []
        t0 = time.time()
        for K in range(1, args.k_max + 1):
            if K > 1 and K * args.min_block_cols > L:
                break
            r = U.evaluate_K(aligned, weights, K, args.min_block_cols, const,
                             "gg", 6, reserved, L, 1.0, 3,
                             n_candidates=args.cut_candidates,
                             max_layer_cols=None,
                             node_budget=args.cut_node_budget,
                             max_library=args.max_library,
                             max_nt=args.max_nt,
                             oligo_overhead=overhead,
                             proxy_candidates=args.proxy_candidates,
                             seed=args.seed,
                             exhaustive_max=args.exhaustive_max)
            if r:
                results.append(r)
            print("    [%s] K=%d done (%.1fs elapsed)" % (name, K, time.time() - t0))
            sys.stdout.flush()
        return results, U.recommend(results) if results else None, time.time() - t0
    finally:
        U.place_cuts = stock


def frontier(name, results, rec, secs, total_w, n):
    print()
    print("  %s -- %.1f s" % (name, secs))
    print("     K  cores      natseq      library  junk%%  oligos   nt ord  "
          "longest  seqs/oligo  stopped by         winner     oligos per layer")
    for r in results:
        mark = "  <==" if rec is not None and r is rec else ""
        per_layer = [len(u) for u in r["layers"]]
        print("    %2d  %2d/%-3d  %4d/%-4d  %11s  %5.1f  %5d  %7s  %5dnt  "
              "%9.2f  %-17s  %-11s %s%s"
              % (r["K"], r["n_cores_encoded"], n, r["encoded_weight"], total_w,
                 format(r["library"], ","), r["junk_pct"], r["oligos"],
                 format(r.get("nt_ordered", r["nt"]), ","),
                 r["longest_oligo_nt"], r["seqs_per_oligo"],
                 str(r.get("stopped_by"))[:17],
                 "%s#%s" % (r.get("winner_arm"), r.get("winner_arm_rank")),
                 per_layer, mark))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aln_fasta", nargs="*")
    ap.add_argument("--max-library", type=int, default=2500)
    ap.add_argument("--max-nt", type=int, default=45000)
    ap.add_argument("--k-max", type=int, default=6)
    ap.add_argument("--min-block-cols", type=int, default=20)
    ap.add_argument("--cut-candidates", type=int, default=50)
    ap.add_argument("--proxy-candidates", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cut-node-budget", type=int, default=2_000_000)
    ap.add_argument("--exhaustive-max", type=int, default=150,
                    help="K=2 enumerates ALL segmentations when the site count is "
                         "below this, which is 331 greedy runs on cluster 1 and "
                         "dominates the wall time; lower it to keep K=2 in pool "
                         "mode and comparable with the other rows")
    ap.add_argument("--arms", default="both", choices=["both", "dfs", "dp"],
                    help="'dp' alone is the only way to reach high K -- the DFS "
                         "cannot enumerate K>=7 in reasonable time")
    args = ap.parse_args()

    label, path = inputs_from_argv([sys.argv[0]] + args.aln_fasta)[0]
    seqs = U.read_aligned_cores(path)
    aligned = [s for s, _ in seqs]
    weights = [w for _, w in seqs]
    L = len(aligned[0])
    const = U.constant_columns(aligned)
    reserved = U.BACKBONE_OVERHANGS
    overhead = 24
    total_w, n = sum(weights), len(aligned)

    banner = {"both": "A/B", "dp": "SINGLE ARM (dp only)",
              "dfs": "SINGLE ARM (dfs only)"}[args.arms]
    print("=" * 104)
    print("%s  %s -- %d cores (%d natural), %d cols"
          % (banner, label, n, total_w, L))
    print("caps: library <= %s, nt <= %s | K 1..%d | %d candidates "
          "(%d proxy + %d random) | seed %d"
          % (format(args.max_library, ","), format(args.max_nt, ","), args.k_max,
             args.cut_candidates, args.proxy_candidates,
             args.cut_candidates - args.proxy_candidates, args.seed))
    print("=" * 104)

    arms = ("dfs", "dp") if args.arms == "both" else (args.arms,)
    out = {}
    for name in arms:
        print()
        print("  running %s arm ..." % name.upper())
        out[name] = run_arm(name, aligned, weights, L, const, args, reserved, overhead)

    for name, title in (("dfs", "STOCK  (place_cuts, branch-and-bound DFS)"),
                        ("dp", "NEW    (dp_cut_search, k-best DP + CSP)")):
        if name not in out:
            continue
        res, rec, secs = out[name]
        frontier(title, res, rec, secs, total_w, n)

    if len(out) < 2:                       # single arm: no comparison to make
        name = arms[0]
        res, rec, secs = out[name]
        if rec:
            print()
            print("  %s recommends K=%d: %d/%d natural, library %s, %d oligos, %s nt"
                  % (name, rec["K"], rec["encoded_weight"], total_w,
                     format(rec["library"], ","), rec["oligos"],
                     format(rec.get("nt_ordered", rec["nt"]), ",")))
        return

    print()
    print("=" * 104)
    print("VERDICT")
    print("=" * 104)
    a_res, a_rec, a_t = out["dfs"]
    b_res, b_rec, b_t = out["dp"]
    by_k = {r["K"]: r for r in a_res}
    print("     K   stock natseq / nt ord      dp natseq / nt ord     change")
    for r in b_res:
        a = by_k.get(r["K"])
        if a is None:
            print("    %2d   %-24s  %d / %s   (K absent from stock)"
                  % (r["K"], "-", r["encoded_weight"],
                     format(r.get("nt_ordered", r["nt"]), ",")))
            continue
        aw, an = a["encoded_weight"], a.get("nt_ordered", a["nt"])
        bw, bn = r["encoded_weight"], r.get("nt_ordered", r["nt"])
        if bw > aw:
            verdict = "DP +%d natural sequences" % (bw - aw)
        elif bw < aw:
            verdict = "*** DP WORSE by %d ***" % (aw - bw)
        elif bn < an:
            verdict = "DP same coverage, -%s nt" % format(an - bn, ",")
        elif bn > an:
            verdict = "DP same coverage, +%s nt" % format(bn - an, ",")
        else:
            verdict = "identical"
        print("    %2d   %4d / %-15s   %4d / %-15s  %s"
              % (r["K"], aw, format(an, ","), bw, format(bn, ","), verdict))

    print()
    for nm, rec in (("stock", a_rec), ("dp", b_rec)):
        if rec:
            print("  %-6s recommends K=%d: %d/%d natural, library %s, %d oligos, %s nt, cuts %s"
                  % (nm, rec["K"], rec["encoded_weight"], total_w,
                     format(rec["library"], ","), rec["oligos"],
                     format(rec.get("nt_ordered", rec["nt"]), ","), rec["cuts"]))
    print()
    print("  wall time: stock %.1f s, dp %.1f s" % (a_t, b_t))


if __name__ == "__main__":
    main()
