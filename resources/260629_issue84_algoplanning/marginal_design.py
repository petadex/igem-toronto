#!/usr/bin/env python3
"""
Marginal-gain fragment/oligo design -- the general (K-fragment) algorithm.

This is the successor to `fragment_split.py`. Two things changed, both driven by
the team's design philosophy (see below):

  1. The number of fragments K is NO LONGER FIXED at 3. We search every
     K = 1 .. K_max (default 12, the realistic wet-lab range) and report the
     whole frontier so we can SEE where junk starts to explode. K-1 cuts are
     placed by a RASPP/SwiftLib-style dynamic program (the general case of
     fragment_split's brute-forced 2 cuts).

  2. The hard junk *budget* is replaced by a JUNK-FRACTION cutoff. The cutoff is
     a single, easy-to-read percentage: the maximum share of the produced library
     that may be junk (phantom recombinants). We add whole cores in order of
     *cheapest junk per newly-encoded natural sequence*, but only while the
     resulting design's junk fraction stays at or below `--max-junk-pct`, and we
     STOP when no core can be added without breaching it. This encodes "smaller
     gains for less junk": we keep spending diversity only while the library
     stays at least (100 - cutoff)% real.

----------------------------------------------------------------------------- #
DESIGN PHILOSOPHY (from the team, 2026-07-01)
----------------------------------------------------------------------------- #
  * We would rather encode a FEW extra full-length sequences per order than
    blow up the junk by cramming in diversity -- whether via too many fragments
    or via over-wide degenerate codons. Coverage is NOT the objective; marginal
    efficiency is.
  * So the headline number this script reports is concrete: how many full-length
    natural sequences (of the whole cluster) the chosen encoding actually
    produces, and at what junk / oligo cost.

----------------------------------------------------------------------------- #
THE MODEL (unchanged from fragment_split.py, generalized to K fragments)
----------------------------------------------------------------------------- #
Cores are aligned. K-1 cuts fall between alignment columns, splitting every core
into K contiguous pieces. Golden Gate / HR assembles one piece per layer, so the
producible library is the CARTESIAN PRODUCT across layers:

    library size = product over fragments f of (# distinct pieces kept in f)
    a natural core is ENCODED  iff  all K of its pieces are in the library
    junk (phantom recombinants) = library size - (# distinct natural cores encoded)

Junk is therefore the inter-fragment cartesian-product junk from the GGAssembler
discussion. Cut placement keeps co-varying columns inside one fragment (so real
haplotypes are preserved and phantom recombinants are minimized); marginal
selection then decides WHICH cores are worth their junk.

Degenerate codons are the finer, intra-fragment lever. They are included here but
JUNK-GATED: a fragment's equal-length pieces are degenerate-compressed only when
the degenerate codon does not enlarge that fragment's library (pure oligo saving,
no extra junk), consistent with the philosophy. Default is pure discrete so the
junk accounting stays exact.

Run:
  # cores must be ALIGNED first (gaps '-'); e.g. in WSL:
  #   mafft --auto cluster1.core.fasta > cluster1.core.aln.fasta
  py marginal_design.py ninetypidorfs/cluster1.core.aln.fasta
  py marginal_design.py ninetypidorfs/cluster2.core.aln.fasta --max-junk-pct 80
  py marginal_design.py ... --k-max 12 --min-block-cols 20 --chemistry agnostic
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import defaultdict
from datetime import datetime
from functools import lru_cache

# Reuse the single-oligo prototype for the optional (junk-gated) degenerate pass.
try:
    from greedy_oligo import (
        design_oligo, expand_degenerate_codon, AA_BY_CODON,
    )
    _HAVE_GREEDY = True
except Exception:  # pragma: no cover
    design_oligo = None
    _HAVE_GREEDY = False


# --------------------------------------------------------------------------- #
# 1. Input: an ALIGNED core FASTA (gaps '-'), headers like >coreN_n<k>
# --------------------------------------------------------------------------- #

def read_aligned_cores(path):
    """Return [(aligned_seq, weight), ...]; all aligned_seq must be equal length.
    The `_n<k>` header suffix is how many natural sequences collapsed onto this
    unique core (its weight)."""
    seqs, header, buf = [], None, []

    def flush():
        if header is not None and buf:
            weight = 1
            if "_n" in header:
                tail = header.rsplit("_n", 1)[1]
                if tail.isdigit():
                    weight = int(tail)
            seqs.append(("".join(buf), weight))

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header, buf = line[1:], []
            else:
                buf.append(line.upper())
    flush()

    lengths = {len(s) for s, _ in seqs}
    if len(lengths) != 1:
        sys.exit(
            f"input is not aligned: found {len(lengths)} distinct lengths {sorted(lengths)}.\n"
            "Align the cores first, e.g.:  mafft --auto <cores.fasta> > <cores.aln.fasta>"
        )
    return seqs


# --------------------------------------------------------------------------- #
# 2. Piece extraction + distinct-count objective (cached over column ranges)
# --------------------------------------------------------------------------- #
# Module-level handles so the lru_cache'd distinct-count can see the alignment
# without rebuilding it. Set once per run by prepare().

_ALIGNED: list[str] = []      # ungapped-later aligned strings, one per unique core
_WEIGHTS: list[int] = []


def prepare(seqs):
    global _ALIGNED, _WEIGHTS
    _ALIGNED = [s for s, _ in seqs]
    _WEIGHTS = [w for _, w in seqs]
    _distinct.cache_clear()


def piece(s, a, b):
    """The real (ungapped) sub-sequence a core contributes over columns [a, b)."""
    return s[a:b].replace("-", "")


@lru_cache(maxsize=None)
def _distinct(a, b):
    """Number of distinct ungapped pieces across all cores over columns [a, b).
    This is the per-fragment library contribution at full coverage; the product
    over fragments is the junk objective the cut placement minimizes."""
    return len({piece(s, a, b) for s in _ALIGNED})


# --------------------------------------------------------------------------- #
# 3. Chemistry: SELF-CONTAINED validity of a cut for Golden Gate / Gibson-HR.
#
# A segmentation is only ever proposed if it is physically buildable, so no
# downstream tool is needed to VALIDATE it. There are two levels of validity, and
# the second is exactly why plain per-edge legality is not enough:
#
#   Level 1 (per site): can a junction physically live at this boundary?
#     GG -- a 4-nt BsaI overhang must sit in immutable (constant) sequence and be
#           individually high-fidelity (not self-complementary, balanced GC). The
#           overhang straddles the cut: last 2 nt of the left residue's codon + the
#           first 2 nt of the right residue's codon. Because we synthesize the
#           constant fragments, we may pick any synonymous codons, so a site offers
#           a SET of achievable overhangs.
#     HR -- a homology arm of >= arm_codons constant residues must straddle the cut.
#
#   Level 2 (per set): the chosen junctions must not cross-react with EACH OTHER.
#     GG -- two overhangs mis-ligate if they present the same sticky end or one
#           anneals to the other's reverse complement within 1 mismatch. This is a
#           property of the whole SET, not one edge -- exactly GGAssembler's
#           "rainbow" constraint -- so we enforce it WHILE searching the cuts
#           (Section 4). It is therefore impossible for this tool to emit a set
#           that mis-assembles.
#     HR -- two homology arms must not be near-identical (would mis-recombine).
#
# Fidelity here is a transparent, literature-grounded MODEL: Potapov 2018 / Pryor
# 2020 show <=1-mismatch ligation dominates infidelity, palindromes self-ligate,
# and extreme-GC overhangs ligate poorly. The exact empirical 256x256 ligation
# table can be dropped into `_overhang_ok` / `_gg_conflict` verbatim without
# touching the algorithm -- its structure is already correct.
# --------------------------------------------------------------------------- #

# Standard genetic code (built exactly like greedy_oligo's, then inverted) so the
# chemistry model is self-contained even if the optional greedy import failed.
_BASES = "TCAG"
_AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_AA_BY_CODON = dict(zip([a + b + c for a in _BASES for b in _BASES for c in _BASES],
                        _AA_STRING))
_CODON_BY_AA = defaultdict(list)
for _cod, _aa in _AA_BY_CODON.items():
    if _aa != "*":
        _CODON_BY_AA[_aa].append(_cod)

_COMP = str.maketrans("ACGT", "TGCA")
GG_OH_LEN = 4                      # BsaI leaves a 4-nt overhang


def _revcomp(s):
    return s.translate(_COMP)[::-1]


def _hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def _const_aa(j):
    """The single residue at column j if it is constant & gap-free, else None."""
    col = {s[j] for s in _ALIGNED}
    if len(col) == 1 and "-" not in col:
        return next(iter(col))
    return None


def _overhang_ok(oh):
    """Level-1 fidelity of a single overhang: reject self-complementary (palindrome
    self-ligation), monotone, and extreme-GC (0 or all) overhangs."""
    if oh == _revcomp(oh):
        return False
    gc = sum(1 for b in oh if b in "GC")
    if gc == 0 or gc == GG_OH_LEN:
        return False
    if len(set(oh)) == 1:
        return False
    return True


def _gg_conflict(a, b):
    """Level-2: two overhangs cross-react if they present the same sticky end or
    one anneals to the other's reverse complement within a single mismatch."""
    return a == b or _hamming(a, _revcomp(b)) <= 1


def _gg_overhangs_at(p, L):
    """Achievable, individually high-fidelity 4-nt overhangs straddling boundary p
    (synonymous codon choice on the two constant flanking residues). Empty set =>
    Level-1 illegal for Golden Gate here."""
    if p <= 0 or p >= L:
        return frozenset()
    laa, raa = _const_aa(p - 1), _const_aa(p)
    if laa is None or raa is None:
        return frozenset()
    ohs = set()
    for lc in _CODON_BY_AA[laa]:
        for rc in _CODON_BY_AA[raa]:
            oh = lc[1:] + rc[:2]        # 2 nt from each side = 4-nt overhang
            if _overhang_ok(oh):
                ohs.add(oh)
    return frozenset(ohs)


def _hr_arm_at(p, L, arm_codons):
    """Constant residue window straddling boundary p (arm_codons//2 each side).
    Returns the arm residue string, or None if the window is not fully constant."""
    half = arm_codons // 2
    if p - half < 0 or p + half > L:
        return None
    residues = []
    for j in range(p - half, p + half):
        aa = _const_aa(j)
        if aa is None:
            return None
        residues.append(aa)
    return "".join(residues)


def _hr_conflict(a, b, max_ident=0.8):
    """Two homology arms risk mis-recombination if near-identical."""
    if len(a) != len(b):
        return False
    same = sum(1 for x, y in zip(a, b) if x == y)
    return same / len(a) >= max_ident


def boundary_tokens(p, L, chemistry, arm_codons):
    """Candidate junction 'tokens' at boundary p: overhang 4-mers (GG, possibly
    several synonymous options), a single homology-arm string (HR), or the trivial
    None token (agnostic). Empty list => Level-1 illegal here for this chemistry."""
    if chemistry == "gg":
        return list(_gg_overhangs_at(p, L))
    if chemistry == "hr":
        arm = _hr_arm_at(p, L, arm_codons)
        return [arm] if arm is not None else []
    return [None]                      # agnostic: no physical junction constraint


def tokens_conflict(a, b, chemistry):
    if a is None or b is None:
        return False
    if chemistry == "gg":
        return _gg_conflict(a, b)
    if chemistry == "hr":
        return _hr_conflict(a, b)
    return False


# --------------------------------------------------------------------------- #
# 4. Cut placement: minimize the product of per-fragment distinct counts for each
#    K, subject to chemistry validity. 'agnostic' is a plain RASPP/SwiftLib DP.
#    'gg'/'hr' need the SET-level orthogonality constraint (Level 2), which is not
#    a per-edge property, so they use an exact DFS that carries the chosen junction
#    tokens along the path (the rainbow-shortest-path idea, done exactly since K
#    and the number of valid sites are small). Returns (cuts, tokens) or None.
# --------------------------------------------------------------------------- #

def place_cuts(L, K, min_block, chemistry, arm_codons):
    if K == 1:
        return [], []                          # one fragment: no junctions at all
    if chemistry == "agnostic":
        cuts = _place_cuts_dp(L, K, min_block)
        return (cuts, [None] * (K - 1)) if cuts is not None else None
    return _place_cuts_orthogonal(L, K, min_block, chemistry, arm_codons)


def _place_cuts_dp(L, K, min_block):
    """RASPP/SwiftLib shortest path: K-1 cuts minimizing sum(log distinct). Every
    boundary is legal (agnostic chemistry makes no physical assumption)."""
    dp = [dict() for _ in range(K + 1)]
    dp[1] = {j: (math.log(_distinct(0, j)), 0)
             for j in range(min_block, L - (K - 1) * min_block + 1)}
    for k in range(2, K + 1):
        cur, prev = dp[k], dp[k - 1]
        for j in range(k * min_block, L - (K - k) * min_block + 1):
            best = None
            for i, (cost_i, _) in prev.items():
                if i > j - min_block:
                    continue
                cost = cost_i + math.log(_distinct(i, j))
                if best is None or cost < best[0]:
                    best = (cost, i)
            if best is not None:
                cur[j] = best
    if L not in dp[K]:
        return None
    cuts, j = [], L
    for k in range(K, 1, -1):
        _, i = dp[k][j]
        cuts.append(i)
        j = i
    cuts.reverse()
    return cuts


def _place_cuts_orthogonal(L, K, min_block, chemistry, arm_codons, node_budget=3_000_000):
    """Exact DFS for gg/hr: place K-1 cuts left-to-right, assigning each a junction
    token (overhang / arm) that is orthogonal to all tokens already chosen, and
    minimize sum(log distinct). Enforcing Level-2 during the search guarantees the
    returned set cannot mis-assemble. Bounded by node_budget as a safety valve."""
    cand = []                                  # Level-1-legal sites and their tokens
    for p in range(min_block, L - min_block + 1):
        toks = boundary_tokens(p, L, chemistry, arm_codons)
        if toks:
            cand.append((p, toks))
    best = {"cost": math.inf, "cuts": None, "tokens": None}
    nodes = [0]

    def dfs(idx_start, last_pos, cuts, tokens, cost):
        nodes[0] += 1
        if nodes[0] > node_budget:
            return
        if len(cuts) == K - 1:                 # close the final fragment [last_pos, L]
            if L - last_pos < min_block:
                return
            total = cost + math.log(_distinct(last_pos, L))
            if total < best["cost"]:
                best.update(cost=total, cuts=list(cuts), tokens=list(tokens))
            return
        # cost only grows (log distinct >= 0), so prune whole subtrees past best
        if cost >= best["cost"]:
            return
        for ci in range(idx_start, len(cand)):
            p, toks = cand[ci]
            if p - last_pos < min_block:
                continue
            # leave room for the remaining fragments after this cut
            frags_left = K - (len(cuts) + 1)
            if L - p < frags_left * min_block:
                continue
            newcost = cost + math.log(_distinct(last_pos, p))
            if newcost >= best["cost"]:
                continue
            for t in toks:
                if all(not tokens_conflict(t, u, chemistry) for u in tokens):
                    cuts.append(p)
                    tokens.append(t)
                    dfs(ci + 1, p, cuts, tokens, newcost)
                    cuts.pop()
                    tokens.pop()

    dfs(0, 0, [], [], 0.0)
    if best["cuts"] is None:
        return None
    return best["cuts"], best["tokens"]


def core_pieces(cuts, L):
    """For the chosen cuts, each core's tuple of K ungapped pieces."""
    bounds = [0] + cuts + [L]
    out = []
    for s in _ALIGNED:
        out.append(tuple(piece(s, bounds[f], bounds[f + 1]) for f in range(len(bounds) - 1)))
    return out


# --------------------------------------------------------------------------- #
# 5. Marginal-gain core selection under a fixed fragmentation.
# --------------------------------------------------------------------------- #

def _coverage(present, pieces, weights):
    """All cores producible from the current piece sets (includes free
    recombinants), and their unique count + natural weight."""
    K = len(present)
    covered = [i for i in range(len(pieces))
               if all(pieces[i][f] in present[f] for f in range(K))]
    return covered, len(covered), sum(weights[i] for i in covered)


def marginal_select(pieces, weights, max_junk_frac, max_junk):
    """Grow a set of cores under a JUNK-FRACTION cutoff. At each step we add the
    core with the lowest phantom-junk cost per newly-encoded natural sequence,
    but ONLY among cores whose addition keeps the design's junk fraction
    (junk / produced library) at or below `max_junk_frac`. We stop when no core
    can be added without breaching the cutoff (or the hard library cap `max_junk`).

    Free recombinants (all pieces already present) are absorbed automatically by
    _coverage, so they never cost anything and only lower the junk fraction.
    Multi-start from the heaviest cores (seed choice barely matters once free
    cores are absorbed) and keep the best run by (covered weight, then least
    junk, then fewest oligos)."""
    n = len(pieces)
    K = len(pieces[0])
    order = sorted(range(n), key=lambda i: -weights[i])
    seeds = order[:min(n, 12)]

    def run_from(seed):
        present = [set() for _ in range(K)]
        for f in range(K):
            present[f].add(pieces[seed][f])
        covered, U, W = _coverage(present, pieces, weights)
        traj = [_snapshot(present, U, W, pieces)]
        cov_set = set(covered)
        while True:
            base_lib = _libsize(present)
            best = None  # (score, marg_seqs, c, add)
            for c in range(n):
                if c in cov_set:
                    continue
                # library if we add c's pieces
                new_lib = 1
                for f in range(K):
                    new_lib *= len(present[f]) + (0 if pieces[c][f] in present[f] else 1)
                if new_lib > max_junk:          # hard safety cap on library size
                    continue
                # what becomes covered if we add c
                add = [(f, pieces[c][f]) for f in range(K) if pieces[c][f] not in present[f]]
                for f, pc in add:
                    present[f].add(pc)
                _, U2, W2 = _coverage(present, pieces, weights)
                for f, pc in add:
                    present[f].discard(pc)
                # feasibility: resulting design must stay under the junk-fraction cutoff
                junk2 = new_lib - U2
                frac2 = junk2 / new_lib if new_lib > 0 else 0.0
                if frac2 > max_junk_frac:
                    continue
                marg_seqs = W2 - W
                if marg_seqs <= 0:
                    continue
                marg_junk = junk2 - (base_lib - U)
                score = marg_junk / marg_seqs   # phantom recombinants per real seq
                key = (score, -marg_seqs, c)
                if best is None or key < best[0]:
                    best = (key, c, add)
            if best is None:
                break
            _, c, add = best
            for f, pc in add:
                present[f].add(pc)
            covered, U, W = _coverage(present, pieces, weights)
            cov_set = set(covered)
            traj.append(_snapshot(present, U, W, pieces))
        return present, U, W, traj

    best = None
    for seed in seeds:
        present, U, W, traj = run_from(seed)
        key = (W, -(_libsize(present) - U), -_n_oligos(present))
        if best is None or key > best[0]:
            best = (key, present, U, W, traj)
    return best[1], best[2], best[3], best[4]


def _libsize(present):
    lib = 1
    for f in present:
        lib *= len(f)
    return lib


def _n_oligos(present):
    return sum(len(f) for f in present)


def _snapshot(present, U, W, pieces):
    lib = _libsize(present)
    return {
        "oligos": _n_oligos(present),
        "library": lib,
        "junk": lib - U,
        "junk_pct": 100.0 * (lib - U) / lib if lib > 0 else 0.0,
        "covered_cores": U,
        "covered_weight": W,
    }


# --------------------------------------------------------------------------- #
# 6. Fragment encoding: discrete (default, exact junk) + junk-gated degenerate.
# --------------------------------------------------------------------------- #

def encode_fragment(present_pieces, degenerate):
    """Turn one fragment's kept distinct pieces into orderable oligos.

    Discrete: one oligo per distinct piece. Junk-gated degenerate: for each group
    of equal-length pieces, try a degenerate codon oligo (greedy_oligo) and accept
    it ONLY if it does not enlarge the fragment's library (library_size <= number
    of pieces in the group) -- pure oligo saving, zero extra junk. Otherwise keep
    the pieces discrete. This honors 'degenerate codons must not blow up junk'."""
    pieces = sorted(present_pieces)
    lengths = sorted({len(p) for p in pieces})
    if not degenerate or not _HAVE_GREEDY:
        return {"n_oligos": len(pieces), "nt": sum(len(p) * 3 for p in pieces),
                "deg_bases": 0, "lengths": lengths, "n_pieces": len(pieces),
                "lib": len(pieces)}

    by_len = defaultdict(list)
    for p in pieces:
        by_len[len(p)].append(p)
    n_oligos, nt, deg_bases, lib = 0, 0, 0, 0
    for width, rows in sorted(by_len.items()):
        if len(rows) == 1:
            n_oligos += 1
            nt += width * 3
            lib += 1
            continue
        res = design_oligo([(r, 1) for r in rows], max_degenerate=3 * width)
        # accept degenerate only if it adds no junk vs. keeping these rows discrete
        if res["library_size"] <= len(rows) and res["n_cores_covered"] == len(rows):
            n_oligos += 1
            nt += len(res["oligo"])
            deg_bases += res["degenerate_bases"]
            lib += res["library_size"]
        else:
            n_oligos += len(rows)
            nt += sum(len(r) * 3 for r in rows)
            lib += len(rows)
    return {"n_oligos": n_oligos, "nt": nt, "deg_bases": deg_bases,
            "lengths": lengths, "n_pieces": len(pieces), "lib": lib}


# --------------------------------------------------------------------------- #
# 7. Evaluate one K end-to-end.
# --------------------------------------------------------------------------- #

def evaluate_K(seqs, K, min_block, chemistry, arm_codons, L, max_junk_frac,
               max_junk, degenerate):
    placed = place_cuts(L, K, min_block, chemistry, arm_codons)
    if placed is None:
        return None                            # no chemistry-valid segmentation for this K
    cuts, tokens = placed
    pieces = core_pieces(cuts, L)
    weights = _WEIGHTS
    present, U, W, traj = marginal_select(pieces, weights, max_junk_frac, max_junk)

    frags = []
    for f in range(K):
        enc = encode_fragment(present[f], degenerate)
        bounds = [0] + cuts + [L]
        enc["a"], enc["b"] = bounds[f], bounds[f + 1]
        frags.append(enc)

    library = _libsize(present)
    n_cores = len(seqs)
    total_w = sum(weights)
    return {
        "K": K,
        "cuts": cuts,
        "tokens": tokens,               # per-junction overhangs (gg) / arms (hr) / None
        "frags": frags,
        "library_size": library,
        "junk": library - U,
        "n_cores_total": n_cores,
        "n_cores_encoded": U,
        "total_weight": total_w,
        "encoded_weight": W,
        "coverage_pct": 100.0 * W / total_w,
        "total_oligos": sum(fr["n_oligos"] for fr in frags),
        "total_nt": sum(fr["nt"] for fr in frags),
        "seqs_per_oligo": W / max(1, sum(fr["n_oligos"] for fr in frags)),
        "junk_pct": 100.0 * (library - U) / library if library > 0 else 0.0,
        "trajectory": traj,
    }


def recommend(results):
    """Pick the K that best fits the philosophy. Every design already respects the
    junk-fraction cutoff, so the remaining trade-off is coverage vs. order size.
    The philosophy is explicit: encode a FEW EXTRA full-length sequences PER ORDER
    rather than pay for more oligos -- i.e. maximize natural sequences encoded per
    oligo ordered. Ties break toward more coverage, then fewer fragments (simpler
    in the wet lab). Note K=1 (synthesize every gene) is the 0%-junk baseline; a
    design only earns the recommendation by beating its sequences-per-oligo."""
    return max(results, key=lambda r: (round(r["seqs_per_oligo"], 3),
                                       r["encoded_weight"], -r["K"]))


# --------------------------------------------------------------------------- #
# 8. Reporting + run persistence
# --------------------------------------------------------------------------- #

def build_report(args, results, rec, L, n_cores, total_w):
    lines = []
    lines.append(f"input: {args.aln_fasta}")
    lines.append(f"loaded {n_cores} unique aligned cores ({total_w} natural sequences), "
                 f"alignment width {L} columns")
    lines.append(f"chemistry: {args.chemistry}   min block cols: {args.min_block_cols}   "
                 f"max junk: {args.max_junk_pct:.0f}% of library   K range: 1..{args.k_max}   "
                 f"degenerate: {args.degenerate}")
    lines.append("")
    lines.append("FRONTIER over number of fragments K "
                 f"(each design keeps junk <= {args.max_junk_pct:.0f}% of its library):")
    lines.append(f"  {'K':>2}  {'encoded cores':>13}  {'nat seqs':>12}  "
                 f"{'library':>10}  {'junk%':>6}  {'oligos':>7}  {'seq/oligo':>9}")
    for r in results:
        star = "  <== recommended" if r["K"] == rec["K"] else ""
        lines.append(f"  {r['K']:>2}  {r['n_cores_encoded']:>4}/{n_cores:<8}  "
                     f"{r['encoded_weight']:>5}/{total_w:<6}  {r['library_size']:>10,}  "
                     f"{r['junk_pct']:>5.1f}%  {r['total_oligos']:>7}  "
                     f"{r['seqs_per_oligo']:>9.2f}{star}")
    lines.append("")
    lines.append("=" * 70)
    lines.append(f"RECOMMENDED DESIGN:  K = {rec['K']} fragments")
    lines.append("=" * 70)
    bounds = [0] + rec["cuts"] + [L]
    seg = "  |  ".join(f"[{bounds[i]},{bounds[i+1]})" for i in range(len(bounds) - 1))
    lines.append(f"cuts / segments: {seg}")
    lines.append("")
    for i, f in enumerate(rec["frags"], 1):
        lvar = (f"{f['lengths'][0]}-{f['lengths'][-1]} aa  <-- length variation (indel)"
                if len(f["lengths"]) > 1 else f"{f['lengths'][0]} aa")
        lines.append(f"  fragment {i}: cols [{f['a']},{f['b']})  "
                     f"{f['n_pieces']} distinct pieces, length {lvar}")
        lines.append(f"      -> {f['n_oligos']} oligos, {f['deg_bases']} degenerate nt, ~{f['nt']} nt")
    lines.append("")
    # Chemistry-validated junctions (empty for K=1 / agnostic).
    if args.chemistry != "agnostic" and rec["cuts"]:
        if args.chemistry == "gg":
            lines.append("JUNCTIONS (Golden Gate) -- validated overhangs, "
                         "Level-1 high-fidelity + Level-2 mutually orthogonal:")
            for cut, oh in zip(rec["cuts"], rec["tokens"]):
                lines.append(f"  col {cut:>4}:  overhang 5'-{oh}-3'  (rc {_revcomp(oh)})")
        else:  # hr
            lines.append("JUNCTIONS (Gibson/HR) -- validated homology arms, "
                         "Level-1 constant window + Level-2 mutually distinct:")
            for cut, arm in zip(rec["cuts"], rec["tokens"]):
                lines.append(f"  col {cut:>4}:  arm {arm}  ({len(arm)} constant residues)")
        lines.append("")
    lines.append("FULL-LENGTH SEQUENCES ENCODED (the headline number):")
    lines.append(f"  {rec['n_cores_encoded']}/{n_cores} unique cores  "
                 f"= {rec['encoded_weight']}/{total_w} natural sequences "
                 f"({rec['coverage_pct']:.0f}%)")
    lines.append(f"  dropped: {n_cores - rec['n_cores_encoded']} cores "
                 f"(adding any would push junk over the {args.max_junk_pct:.0f}% cutoff)")
    lines.append("")
    lines.append(f"producible library size: {rec['library_size']:,}")
    lines.append(f"phantom recombinants (junk): {rec['junk']:,}  "
                 f"({rec['junk_pct']:.1f}% of the produced library)")
    lines.append(f"total oligos to order: {rec['total_oligos']}  "
                 f"({rec['seqs_per_oligo']:.2f} natural seqs per oligo)")
    lines.append(f"total synthesis (GGAssembler-style): ~{rec['total_nt']} nt")
    lines.append("")
    lines.append("marginal trajectory of the recommended design "
                 "(each step = cheapest next core added):")
    lines.append(f"      {'oligos':>7}  {'encoded':>7}  {'nat seqs':>8}  "
                 f"{'library':>10}  {'junk%':>6}")
    last = None
    for s in rec["trajectory"]:
        if last is None or s["covered_weight"] != last:
            lines.append(f"      {s['oligos']:>7}  {s['covered_cores']:>7}  "
                         f"{s['covered_weight']:>8}  {s['library']:>10,}  "
                         f"{s['junk_pct']:>5.1f}%")
            last = s["covered_weight"]
    lines.append("")
    if args.chemistry == "agnostic":
        lines.append("note: chemistry is 'agnostic' -- cuts are NOT validated for any "
                     "assembly method. Re-run with --chemistry gg (or hr) for a design "
                     "that is buildable as-is.")
    else:
        lines.append(f"SELF-CONTAINED: every junction above carries a validated "
                     f"{'Golden Gate overhang' if args.chemistry=='gg' else 'homology arm'} "
                     f"(Level-1) and the whole set is mutually {'orthogonal' if args.chemistry=='gg' else 'distinct'} "
                     f"(Level-2), so this design assembles as specified with NO GGAssembler "
                     f"validation step. GGAssembler would only be an optional extra to shave "
                     f"codon cost further; it is not required for correctness.")
    return "\n".join(lines)


def save_run(out_root, stem, args, results, rec, L, seqs):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(out_root, f"{ts}_{stem}_{args.chemistry}_marginal")
    os.makedirs(run_dir, exist_ok=True)
    report = build_report(args, results, rec, L, len(seqs), sum(w for _, w in seqs))
    with open(os.path.join(run_dir, "report.txt"), "w") as fh:
        fh.write(report + "\n")
    summary = {
        "input": os.path.abspath(args.aln_fasta),
        "alignment_width": L,
        "n_cores_total": len(seqs),
        "total_weight": sum(w for _, w in seqs),
        "chemistry": args.chemistry,
        "min_block_cols": args.min_block_cols,
        "max_junk_pct": args.max_junk_pct,
        "k_max": args.k_max,
        "degenerate": args.degenerate,
        "recommended_K": rec["K"],
        "frontier": [{k: r[k] for k in (
            "K", "cuts", "n_cores_encoded", "encoded_weight", "junk",
            "total_oligos", "total_nt", "library_size", "seqs_per_oligo",
            "junk_pct", "coverage_pct")} for r in results],
        "recommended": {k: rec[k] for k in (
            "K", "cuts", "tokens", "frags", "n_cores_encoded", "encoded_weight",
            "total_weight", "junk", "junk_pct", "library_size", "total_oligos",
            "total_nt", "coverage_pct")},
    }
    with open(os.path.join(run_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    # one FASTA per fragment of the recommended design (what you'd order)
    bounds = [0] + rec["cuts"] + [L]
    aligned = [s for s, _ in seqs]
    for i in range(len(bounds) - 1):
        a, b = bounds[i], bounds[i + 1]
        # recover this fragment's kept pieces from the encoded cores
        kept = sorted({piece(s, a, b) for s in aligned})  # all distinct; report set
        with open(os.path.join(run_dir, f"fragment{i+1}.fasta"), "w") as fh:
            for k, frag in enumerate(kept, 1):
                fh.write(f">frag{i+1}_p{k}_len{len(frag)}\n{frag}\n")
    return run_dir, report


# --------------------------------------------------------------------------- #
# 9. CLI
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aln_fasta", help="ALIGNED core FASTA (gaps '-'), headers >coreN_n<k>")
    ap.add_argument("--k-max", type=int, default=12,
                    help="max number of fragments to try (default 12; K sweeps 1..k-max)")
    ap.add_argument("--min-block-cols", type=int, default=20,
                    help="minimum fragment width in alignment columns (default 20)")
    ap.add_argument("--chemistry", choices=["agnostic", "gg", "hr"], default="agnostic",
                    help="cut-site validity model: 'gg' Golden Gate (validated 4-nt "
                         "overhangs, orthogonal), 'hr' Gibson/HR (validated homology "
                         "arms), 'agnostic' no assembly assumption (default)")
    ap.add_argument("--arm-codons", type=int, default=6,
                    help="HR only: constant residues required across a homology arm "
                         "(default 6 ~= 18 bp)")
    ap.add_argument("--max-junk-pct", type=float, default=80.0,
                    help="JUNK CUTOFF: max percent of the produced library allowed to "
                         "be junk (phantom recombinants), 0-100 (default 80). Lower = "
                         "less junk, fewer sequences; higher = more coverage, more junk.")
    ap.add_argument("--max-junk", type=int, default=1_000_000,
                    help="hard safety cap on producible library size (default 1e6)")
    ap.add_argument("--degenerate", action="store_true",
                    help="enable junk-gated degenerate-codon oligo compression "
                         "(only applied when it adds no junk)")
    ap.add_argument("--out-dir", default="algoruns",
                    help="parent folder for per-run output subfolders")
    args = ap.parse_args()

    seqs = read_aligned_cores(args.aln_fasta)
    prepare(seqs)
    L = len(seqs[0][0])

    max_junk_frac = args.max_junk_pct / 100.0
    results = []
    for K in range(1, args.k_max + 1):
        if K > 1 and K * args.min_block_cols > L:
            break                       # can't fit this many fragments
        r = evaluate_K(seqs, K, args.min_block_cols, args.chemistry, args.arm_codons,
                       L, max_junk_frac, args.max_junk, args.degenerate)
        if r is not None:
            results.append(r)
    if not results:
        sys.exit("no chemistry-valid segmentation for any K (try smaller "
                 "--min-block-cols, a different --chemistry, or --arm-codons).")

    rec = recommend(results)
    stem = os.path.splitext(os.path.basename(args.aln_fasta))[0]
    run_dir, report = save_run(args.out_dir, stem, args, results, rec, L, seqs)
    print(report)
    print(f"\nsaved run to {run_dir}/  (report.txt, summary.json, fragment*.fasta)")


if __name__ == "__main__":
    main()
