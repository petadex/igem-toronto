"""Exact cut search: k-best DP over the RASPP graph, then a junction CSP.

WHY THIS EXISTS
---------------
`place_cuts()` in cutsearch_design.py ranks segmentations by a cost that is
ADDITIVE over blocks -- log(distinct_count(a, b)) depends only on that block's
two endpoints.  That is RASPP's structure exactly (Endelman et al. 2004) and
would be a shortest path, except that Golden Gate junctions must be mutually
orthogonal: a constraint over the whole token SET, with no edge-local
expression.  Bellman's principle fails, so that module enumerates paths by
branch-and-bound DFS instead.

The enumeration is exponential.  On cluster 1 (101 legal sites) the number of
complete placements is ~5e3 at K=3, ~1.7e5 at K=4 and ~4.1e6 at K=5 against a
2,000,000 node budget -- so K>=5 truncates, returns "best found" rather than
best, and the truncation is BIASED: the DFS walks sites left to right, so
segmentations with a late first cut are never reached at all.

This module separates the two concerns that were tangled together:

    1. WHERE to cut  -- exact k-best DP over the same graph, tokens ignored.
                        Tokens are what break optimal substructure; drop them
                        and RASPP's DP applies unchanged.
    2. WHICH overhangs -- a tiny backtracking CSP over the chosen cuts only.

Walking segmentations in increasing DP cost and taking the first feasible ones
gives the PROVABLE GLOBAL OPTIMUM, not a heuristic.  Cluster 1 at K=6 is 507
nodes and ~17k edges, under 1 MB -- about 120x less work than the truncated
search it replaces.

cutsearch_design.py is imported read-only and never modified; this is the same
pattern runtime.py, poolcost.py and annotate_degenerate.py already use.

  python dp_cutsearch.py                       # whichever default inputs exist
  python dp_cutsearch.py cores.aln.fasta       # explicit
  python dp_cutsearch.py cores.aln.fasta --k-max 4 --no-compare
"""

import argparse
import heapq
import math
import sys
import time

from _paths import inputs_from_argv          # noqa: F401  (also puts HERE on path)
import cutsearch_design as U


# =========================================================================== #
# STAGE C -- assigning overhangs to a FIXED set of cuts.
#
# This is the part that cannot live in the graph.  With the cut positions fixed
# it is a tiny constraint-satisfaction problem: K-1 variables, each with the
# domain of legal overhangs at its site (2-30 entries), and pairwise
# non-conflict constraints.  K-1 <= 11, so backtracking settles it immediately.
# =========================================================================== #

def assign_tokens(domains, chemistry, stats=None):
    """One mutually-orthogonal token per cut, or None if no assignment exists.

    Returns the tokens in CUT ORDER (indexed as `domains` is), not in the order
    the search happened to bind them.  Variables are tried smallest-domain-first
    (static MRV): a site with two options is far likelier to be the one that
    fails, so testing it early prunes the tree instead of discovering the
    conflict at the last level."""
    n = len(domains)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: len(domains[i]))
    chosen = [None] * n

    def bt(k):
        if stats is not None:
            stats["csp_nodes"] += 1
        if k == n:
            return True
        i = order[k]
        for t in domains[i]:
            clash = False
            for j in range(n):
                cj = chosen[j]
                if cj is not None and U.tokens_conflict(t, cj, chemistry):
                    clash = True
                    break
            if clash:
                continue
            chosen[i] = t
            if bt(k + 1):
                return True
            chosen[i] = None
        return False

    return chosen if bt(0) else None


class PairFilter:
    """Is there ANY non-conflicting token pair between two sites?

    Full backtracking is wasted on the common failure, which is a single pair of
    sites that cannot coexist at all.  This is arc consistency over pairs,
    computed lazily and cached, and it is what makes correlated failure cheap:
    near-duplicate segmentations reuse the same sites, so one dead pair kills a
    whole block of top-ranked paths and we want to find that out in O(1)."""

    def __init__(self, domains, chemistry):
        self.domains = domains
        self.chemistry = chemistry
        self._memo = {}

    def ok(self, i, j):
        key = (i, j) if i < j else (j, i)
        hit = self._memo.get(key)
        if hit is None:
            hit = any(not U.tokens_conflict(t, u, self.chemistry)
                      for t in self.domains[key[0]]
                      for u in self.domains[key[1]])
            self._memo[key] = hit
        return hit

    def any_dead_pair(self, idxs):
        for a in range(len(idxs)):
            for b in range(a + 1, len(idxs)):
                if not self.ok(idxs[a], idxs[b]):
                    return True
        return False


# =========================================================================== #
# STAGE A -- k-best DP over the graph, tokens ignored.
#
# Nodes are (site index, cuts placed so far); edges are blocks; edge weight is
# log(distinct_count).  Keeping the top `beam` partial paths at each node rather
# than a single best yields the true top-`beam` complete paths:
#
#   if a globally top-`beam` path had a prefix outside its node's top-`beam`,
#   each of the `beam` better prefixes could take that path's suffix, producing
#   `beam` strictly better complete paths -- contradiction.
#
# That argument needs suffix cost to be independent of how the node was reached,
# which is true ONLY because tokens are not part of the state here.  That is the
# whole reason the DP/CSP split works.
# =========================================================================== #

def _kbest_paths(aligned, L, K, min_block, sites, beam, max_layer_cols):
    """All complete segmentations the beam can see, cheapest first.

    Returns [(cost, (p1, ..., p_{K-1}), (i1, ..., i_{K-1})), ...] with the second
    tuple the cut columns and the third their indices into `sites`."""
    pos = [p for p, _ in sites]
    S = len(pos)

    def width_ok(w):
        return w >= min_block and (max_layer_cols is None or w <= max_layer_cols)

    # layer[i] = up to `beam` (cost, cuts_tuple, idx_tuple) for reaching site i
    # having placed the current number of cuts.  Rebuilt in place per level.
    layer = [[] for _ in range(S)]
    for i, p in enumerate(pos):
        if not width_ok(p):
            continue
        # forward room: K-1 blocks still to place after this cut
        if L - p < (K - 1) * min_block:
            continue
        layer[i] = [(math.log(U.distinct_count(aligned, 0, p)), (p,), (i,))]

    for j in range(2, K):                       # placing the j-th cut
        nxt = [[] for _ in range(S)]
        remaining = K - j                       # blocks left after this cut
        for i, p in enumerate(pos):
            if L - p < remaining * min_block:
                continue
            cand = []
            for h in range(i):
                if not layer[h]:
                    continue
                w = p - pos[h]
                if w < min_block:
                    continue
                if max_layer_cols is not None and w > max_layer_cols:
                    continue                    # pos ascends, but earlier h are WIDER
                step = math.log(U.distinct_count(aligned, pos[h], p))
                for c, cuts, idxs in layer[h]:
                    cand.append((c + step, cuts + (p,), idxs + (i,)))
            if cand:
                nxt[i] = heapq.nsmallest(beam, cand, key=lambda r: r[0])
        layer = nxt

    out = []
    for i, p in enumerate(pos):
        if not layer[i]:
            continue
        w = L - p
        if not width_ok(w):
            continue
        close = math.log(U.distinct_count(aligned, p, L))
        for c, cuts, idxs in layer[i]:
            out.append((c + close, cuts, idxs))
    out.sort(key=lambda r: r[0])
    return out


# =========================================================================== #
# STAGE B -- walk in cost order, CSP each, double the beam if the pool is short.
# =========================================================================== #

def dp_cut_search(aligned, L, K, min_block, const, chemistry, arm_codons,
                  reserved, node_budget=None, n_keep=1, pool_factor=10,
                  max_layer_cols=None, beam=None, beam_max=8192, stats=None):
    """Drop-in replacement for `place_cuts` with the same return contract.

    Returns (candidates, truncated) with candidates = [(cost, cuts, tokens), ...]
    cheapest first.  `node_budget` is accepted so this is signature-compatible
    and IGNORED -- there is no enumeration to bound.

    `truncated` keeps its old meaning ("the pool is not provably the best") but
    fires for a different reason: only when `beam_max` was reached without
    filling the pool.  Infeasible paths ABOVE the winner are harmless -- with
    beam m the DP yields the true top-m, so a feasible path at rank r <= m is the
    global optimum however many infeasible ones preceded it."""
    if stats is None:
        stats = new_stats()

    if K == 1:
        if max_layer_cols is not None and L > max_layer_cols:
            return [], False
        return [(0.0, [], [])], False

    sites = U.legal_sites(L, min_block, const, chemistry, arm_codons, reserved)
    if len(sites) < K - 1:
        return [], False

    domains = [toks for _p, toks in sites]
    pf = PairFilter(domains, chemistry)
    pool_size = max(1, n_keep * pool_factor)
    if beam is None:
        beam = max(200, 4 * pool_size)

    kept, truncated, prev_paths = [], False, -1
    while True:
        stats["dp_passes"] += 1
        t0 = time.time()
        paths = _kbest_paths(aligned, L, K, min_block, sites, beam, max_layer_cols)
        stats["dp_seconds"] += time.time() - t0
        stats["paths_seen"] = len(paths)

        kept, skipped = [], 0
        for cost, cuts, idxs in paths:
            if pf.any_dead_pair(idxs):
                stats["prefiltered"] += 1
                skipped += 1
                continue
            toks = assign_tokens([domains[i] for i in idxs], chemistry, stats)
            stats["csp_calls"] += 1
            if toks is None:
                stats["csp_failed"] += 1
                skipped += 1
                continue
            if not kept:
                stats["first_feasible_rank"] = skipped        # 0 == rank 1 won
            kept.append((cost, list(cuts), toks))
            if len(kept) >= pool_size:
                break

        if len(kept) >= pool_size or not paths:
            break
        # A short pool is NOT evidence of truncation on its own: the caller may
        # simply have asked for more candidates than exist.  K=2 in exhaustive
        # mode asks for 1010 when only ~101 segmentations exist, and an earlier
        # version doubled the beam to the ceiling and then reported truncation
        # for a search that had in fact enumerated everything.  Growing the beam
        # without finding new paths is the real signal that the space is
        # exhausted, and that case is complete, not truncated.
        if len(paths) == prev_paths:
            break
        prev_paths = len(paths)
        if beam >= beam_max:
            truncated = True
            break
        beam *= 2
        stats["beam_doublings"] += 1

    stats["beam_final"] = beam
    stats["pool"] = len(kept)
    return kept, truncated


def new_stats():
    return {"dp_passes": 0, "dp_seconds": 0.0, "paths_seen": 0, "prefiltered": 0,
            "csp_calls": 0, "csp_failed": 0, "csp_nodes": 0, "beam_doublings": 0,
            "beam_final": 0, "pool": 0, "first_feasible_rank": None}


# =========================================================================== #
# Verification and comparison.
# =========================================================================== #

def check_candidates(cands, L, K, min_block, sites_by_pos, chemistry):
    """Every returned candidate must be structurally legal.  Catches a CSP that
    returns a malformed assignment, which no downstream check would notice until
    the wet lab did."""
    for cost, cuts, toks in cands:
        assert len(cuts) == K - 1, "wrong number of cuts"
        assert len(toks) == K - 1, "tokens do not match cuts"
        assert list(cuts) == sorted(cuts), "cuts not ascending"
        b = [0] + list(cuts) + [L]
        for i in range(len(b) - 1):
            assert b[i + 1] - b[i] >= min_block, "block under min_block"
        for p, t in zip(cuts, toks):
            assert t in sites_by_pos[p], "token not legal at its site"
        for a in range(len(toks)):
            for c in range(a + 1, len(toks)):
                assert not U.tokens_conflict(toks[a], toks[c], chemistry), \
                    "returned tokens conflict"
    return True


def brute_force(aligned, L, K, min_block, sites, chemistry):
    """Every legal segmentation with a feasible token assignment, cheapest first.
    Only tractable on toys -- this is the ground truth the DP is checked against."""
    from itertools import combinations
    pos = [p for p, _ in sites]
    domains = [toks for _p, toks in sites]
    out = []
    for idxs in combinations(range(len(pos)), K - 1):
        cuts = [pos[i] for i in idxs]
        b = [0] + cuts + [L]
        if any(b[i + 1] - b[i] < min_block for i in range(len(b) - 1)):
            continue
        if assign_tokens([domains[i] for i in idxs], chemistry) is None:
            continue
        cost = sum(math.log(U.distinct_count(aligned, b[i], b[i + 1]))
                   for i in range(len(b) - 1))
        out.append((cost, cuts))
    out.sort(key=lambda r: r[0])
    return out


TOY = ["MNPTDGSAKLPVART", "MNPTDGSAKLPVSRT", "MNPTDGSAKLPVTRT",
       "MNPTEGSARLPVART", "MNPTEGSARLPVSRT", "MNPTEGSARLPVTRT"]


def self_test():
    """The 6-core / 15-column example from figures/cutsearch_example.md, where
    the answer is known by hand: (3,10) and (3,11) cost 6, everything else 12."""
    print("=" * 72)
    print("SELF TEST -- toy from figures/cutsearch_example.md")
    print("=" * 72)
    al, L, K, mb = TOY, len(TOY[0]), 3, 3
    const = U.constant_columns(al)
    sites = U.legal_sites(L, mb, const, "gg", 6, U.BACKBONE_OVERHANGS)
    sites_by_pos = {p: toks for p, toks in sites}
    print("  sites: %s" % [p for p, _ in sites])

    bf = brute_force(al, L, K, mb, sites, "gg")
    print("  brute force, all %d feasible segmentations:" % len(bf))
    for cost, cuts in bf:
        print("     %-10s product %5.1f" % (str(cuts), math.exp(cost)))

    cands, trunc = dp_cut_search(al, L, K, mb, const, "gg", 6,
                                 U.BACKBONE_OVERHANGS, n_keep=len(bf),
                                 pool_factor=1)
    check_candidates(cands, L, K, mb, sites_by_pos, "gg")
    # Sort both with the SAME tie-break before comparing.  Costs tie constantly
    # here (six of the eight segmentations score 12), and the two enumerations
    # visit tied paths in different orders -- which is not a disagreement about
    # anything, so comparing raw order would fail for no reason.
    dp = sorted((round(c, 9), list(cu)) for c, cu, _t in cands)
    gt = sorted((round(c, 9), list(cu)) for c, cu in bf)
    ok = dp == gt
    print("  DP returned %d, truncated=%s" % (len(cands), trunc))
    print("  ordering identical to brute force: %s" % ("YES" if ok else "NO"))
    if not ok:
        print("     dp: %s" % dp[:4])
        print("     bf: %s" % gt[:4])
    cheapest = math.exp(cands[0][0])
    print("  cheapest product = %.1f (expected 6)" % cheapest)

    # The retry path must actually be exercised, or it could be silently broken.
    print()
    print("  forced-failure test (every domain crushed to one shared token):")
    dom = [t for _p, t in sites]
    shared = dom[0][0]
    fake = [[shared] for _ in sites]
    pf = PairFilter(fake, "gg")
    dead = pf.any_dead_pair((0, 1))
    none_assign = assign_tokens([fake[0], fake[1]], "gg")
    print("     identical tokens at two sites -> dead pair: %s, CSP: %s"
          % (dead, "None" if none_assign is None else "assigned"))
    assert dead and none_assign is None, "pre-filter/CSP failed to reject"
    print("     both correctly rejected")
    return ok


def compare(label, path, k_max, compare_dfs, min_block, node_budget):
    seqs = U.read_aligned_cores(path)
    al = [s for s, _ in seqs]
    L = len(al[0])
    const = U.constant_columns(al)
    sites = U.legal_sites(L, min_block, const, "gg", 6, U.BACKBONE_OVERHANGS)
    sites_by_pos = {p: toks for p, toks in sites}

    print()
    print("=" * 108)
    print("%s -- %d cores, %d cols, %d legal sites (min_block %d)"
          % (label, len(al), L, len(sites), min_block))
    print("=" * 108)
    hdr = ("   K |      DP cost   pool  trunc   rank  pre  csp!  x2   "
           "DP sec |    DFS cost  pool  trunc   DFS sec |  verdict")
    print(hdr)
    print("   " + "-" * (len(hdr) - 3))

    for K in range(2, k_max + 1):
        st = new_stats()
        t0 = time.time()
        cands, trunc = dp_cut_search(al, L, K, min_block, const, "gg", 6,
                                     U.BACKBONE_OVERHANGS, n_keep=5,
                                     pool_factor=10, stats=st)
        dt = time.time() - t0
        if cands:
            check_candidates(cands, L, K, min_block, sites_by_pos, "gg")
        dp_cost = cands[0][0] if cands else float("nan")
        rank = st["first_feasible_rank"]

        d_cost, d_pool, d_tr, d_dt, verdict = float("nan"), 0, "-", float("nan"), ""
        if compare_dfs:
            t0 = time.time()
            dc, d_tr_b = U.place_cuts(al, L, K, min_block, const, "gg", 6,
                                      U.BACKBONE_OVERHANGS,
                                      node_budget=node_budget, n_keep=5)
            d_dt = time.time() - t0
            d_cost = dc[0][0] if dc else float("nan")
            d_pool, d_tr = len(dc), "yes" if d_tr_b else "no"
            if dc and cands:
                if dp_cost < d_cost - 1e-9:
                    verdict = "DP BETTER by %.4f" % (d_cost - dp_cost)
                elif dp_cost > d_cost + 1e-9:
                    verdict = "*** DP WORSE -- BUG ***"
                else:
                    verdict = "equal"

        print("   %d | %11.4f %6d  %5s %6s %4d %5d %3d %8.2f | %11.4f %5d  %5s "
              "%9.2f | %s"
              % (K, dp_cost, len(cands), "yes" if trunc else "no",
                 "-" if rank is None else str(rank + 1), st["prefiltered"],
                 st["csp_failed"], st["beam_doublings"], dt,
                 d_cost, d_pool, d_tr, d_dt, verdict))

    print()
    print("   rank = position of the first FEASIBLE path in DP cost order "
          "(1 = the optimum was directly assignable)")
    print("   pre  = paths killed by the pairwise pre-filter;  csp! = paths that "
          "reached the CSP and failed;  x2 = beam doublings")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aln_fasta", nargs="*", help="ALIGNED core FASTA(s)")
    ap.add_argument("--k-max", type=int, default=6)
    ap.add_argument("--min-block-cols", type=int, default=20)
    ap.add_argument("--no-compare", action="store_true",
                    help="skip the place_cuts baseline (it is slow at K>=5)")
    ap.add_argument("--no-self-test", action="store_true")
    ap.add_argument("--cut-node-budget", type=int, default=2_000_000,
                    help="node budget for the place_cuts baseline only")
    args = ap.parse_args()

    ok = True
    if not args.no_self_test:
        ok = self_test()

    pairs = inputs_from_argv([sys.argv[0]] + args.aln_fasta)
    for label, path in pairs:
        compare(label, path, args.k_max, not args.no_compare,
                args.min_block_cols, args.cut_node_budget)

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
