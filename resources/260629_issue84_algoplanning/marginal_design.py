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

Degenerate codons are the finer, intra-fragment lever (`--degenerate`). They are
JUNK-BUDGETED: starting from the discrete design, column-homologous piece groups
are compressed into degenerate oligos -- biggest oligo saving first -- while the
whole design's library stays under the SAME --max-junk-pct cutoff, spending only
the junk head-room the cutoff leaves. This shrinks the oligo order by trading a
little phantom junk for fewer oligos (it does not buy new coverage). Default is
pure discrete so the junk accounting stays exact.

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
        design_oligo, expand_degenerate_codon, AA_BY_CODON, IUPAC, best_codon_for,
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


# --------------------------------------------------------------------------- #
# Type IIS domestication + reserved backbone overhangs (Golden Gate).
#
#   * FORBIDDEN_SITES: the assembly enzyme's own recognition site must not occur
#     ANYWHERE inside a produced full-length sequence, or the enzyme would cut it
#     internally. Default BsmBI/Esp3I = CGTCTC (+ its reverse complement GAGACG);
#     overridden from --gg-enzyme in main(). Enforced by actually BACK-TRANSLATING
#     the design to concrete codons (Section 6b): synonymous codons are chosen so
#     no expansion of any assembled full-length sequence contains the site.
#   * BACKBONE_OVERHANGS: the destination vector opens with these two fusion sites
#     to receive the insert's outer ends, so INTERNAL junctions must not use them
#     (nor cross-react with them) -- unless we deliberately share them in a
#     shared-overhang / minimal-plasmid design (then the first/last fragment, which
#     mates the backbone, is excluded from the ban).
# --------------------------------------------------------------------------- #
FORBIDDEN_SITES = frozenset({"CGTCTC", "GAGACG"})
BACKBONE_OVERHANGS = frozenset({"CGGA", "GGTG"})


def _has_forbidden_site(nt):
    """True if a concrete DNA string contains any forbidden Type IIS site."""
    return any(site in nt for site in FORBIDDEN_SITES)


def _degenerate_has_forbidden(oligo_iupac):
    """True if SOME expansion of an IUPAC-coded oligo could contain a forbidden
    site -- i.e. a window whose per-position ambiguity sets can all match it. Used
    to reject a degenerate-codon oligo (which carries IUPAC ambiguity codes, not
    plain ACGT) unless every expansion is site-free."""
    for site in FORBIDDEN_SITES:
        for i in range(len(oligo_iupac) - len(site) + 1):
            if all(site[k] in IUPAC[oligo_iupac[i + k]] for k in range(len(site))):
                return True
    return False


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


def _gg_overhangs_at(p, L, reserved=frozenset()):
    """Achievable, individually high-fidelity 4-nt overhangs straddling boundary p
    (synonymous codon choice on the two constant flanking residues). An overhang is
    kept only if SOME synonymous flanking-codon pair realizes it without spelling a
    forbidden Type IIS site in the junction (domesticable), and it does not cross-
    react with any reserved backbone overhang. Empty set => illegal for GG here."""
    if p <= 0 or p >= L:
        return frozenset()
    laa, raa = _const_aa(p - 1), _const_aa(p)
    if laa is None or raa is None:
        return frozenset()
    ohs = set()
    for lc in _CODON_BY_AA[laa]:
        for rc in _CODON_BY_AA[raa]:
            oh = lc[1:] + rc[:2]        # 2 nt from each side = 4-nt overhang
            if not _overhang_ok(oh):
                continue
            if _has_forbidden_site(lc + rc):
                continue                # this synonymous context spells the site
            if any(_gg_conflict(oh, r) for r in reserved):
                continue                # collides with a reserved backbone overhang
            ohs.add(oh)                 # kept: >=1 domesticated context exists
    return frozenset(ohs)


def _gg_pin_codons(p, oh):
    """Recover a concrete, domesticated synonymous codon pair (lc, rc) for the two
    constant residues flanking boundary p that realizes overhang `oh` with no
    forbidden site in lc+rc. These become the pinned junction codons (last codon of
    the left fragment, first codon of the right fragment). Returns (lc, rc)."""
    laa, raa = _const_aa(p - 1), _const_aa(p)
    for lc in _CODON_BY_AA[laa]:
        for rc in _CODON_BY_AA[raa]:
            if lc[1:] + rc[:2] == oh and not _has_forbidden_site(lc + rc):
                return lc, rc
    raise ValueError(f"no domesticated codon pair for overhang {oh} at col {p}")


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


def boundary_tokens(p, L, chemistry, arm_codons, reserved=frozenset()):
    """Candidate junction 'tokens' at boundary p: overhang 4-mers (GG, possibly
    several synonymous options), a single homology-arm string (HR), or the trivial
    None token (agnostic). Empty list => Level-1 illegal here for this chemistry."""
    if chemistry == "gg":
        return list(_gg_overhangs_at(p, L, reserved))
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

def place_cuts(L, K, min_block, chemistry, arm_codons, reserved=frozenset()):
    if K == 1:
        return [], []                          # one fragment: no junctions at all
    if chemistry == "agnostic":
        cuts = _place_cuts_dp(L, K, min_block)
        return (cuts, [None] * (K - 1)) if cuts is not None else None
    return _place_cuts_orthogonal(L, K, min_block, chemistry, arm_codons, reserved)


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


def _place_cuts_orthogonal(L, K, min_block, chemistry, arm_codons, reserved=frozenset(),
                           node_budget=3_000_000):
    """Exact DFS for gg/hr: place K-1 cuts left-to-right, assigning each a junction
    token (overhang / arm) that is orthogonal to all tokens already chosen, and
    minimize sum(log distinct). Enforcing Level-2 during the search guarantees the
    returned set cannot mis-assemble. Bounded by node_budget as a safety valve."""
    cand = []                                  # Level-1-legal sites and their tokens
    for p in range(min_block, L - min_block + 1):
        toks = boundary_tokens(p, L, chemistry, arm_codons, reserved)
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
# 5b. Densify: degeneracy as a COVERAGE move (DeCoDe-lite), linkage-guided.
#
# #1 spends degeneracy only to shrink the oligo order. #2 lets it BUY coverage: a
# degenerate codon over an INDEPENDENTLY-varying column encodes the whole cartesian
# of natural haplotypes at that column, so it can cover cores that discrete
# selection dropped -- with little junk BECAUSE the column is independent (the
# cartesian is densely natural). Columns that co-vary (LINKED) are kept as discrete
# haplotypes instead, so we never manufacture off-haplotype phantoms.
#
# Mechanics: within each fragment, cores split into CELLS keyed by (gap pattern,
# linked-column haplotype). A cell can be encoded discretely (its kept pieces, one
# oligo each) or as ONE degenerate oligo that degenerates the cell's INDEPENDENT
# columns -- covering every natural piece of the cell. We greedily adopt degenerate
# cells that add the most newly-covered natural sequence (then fewest oligos, least
# junk) while the whole design's junk fraction stays <= the cutoff. Junk stays
# exact (variant count = codon cartesian). Linkage is scored with APC-corrected
# mutual information; the junk budget is the ultimate guard, so MI errors are safe.
# --------------------------------------------------------------------------- #

def _pair_mi(a_res, b_res, weights):
    """Weighted mutual information (bits) between two aligned columns."""
    total = sum(weights)
    if total <= 0:
        return 0.0
    pa, pb, pab = defaultdict(float), defaultdict(float), defaultdict(float)
    for a, b, w in zip(a_res, b_res, weights):
        pa[a] += w
        pb[b] += w
        pab[(a, b)] += w
    mi = 0.0
    for (a, b), c in pab.items():
        pxy = c / total
        mi += pxy * math.log2(pxy / ((pa[a] / total) * (pb[b] / total)))
    return max(0.0, mi)


def _independent_cols(variable, colres, weights, tau):
    """Variable columns whose max APC-corrected MI to any other variable column is
    <= tau (near-independent -> safe to degenerate). `colres[r]` = residues (one per
    core, aligned to `weights`) at relative column r."""
    if len(variable) <= 1:
        return set(variable)
    mi = {}
    for x in range(len(variable)):
        for y in range(x + 1, len(variable)):
            i, j = variable[x], variable[y]
            mi[(i, j)] = mi[(j, i)] = _pair_mi(colres[i], colres[j], weights)
    mean_i = {i: sum(mi[(i, j)] for j in variable if j != i) / (len(variable) - 1)
              for i in variable}
    overall = sum(mean_i.values()) / len(variable)
    indep = set()
    for i in variable:
        worst = 0.0
        for j in variable:
            if j == i:
                continue
            apc = mi[(i, j)] - (mean_i[i] * mean_i[j] / overall if overall > 0 else 0.0)
            worst = max(worst, apc)
        if worst <= tau:
            indep.add(i)
    return indep


def _build_group_oligo(block_rows, deg_cols, lead, tail, lctx, rctx):
    """One degenerate IUPAC oligo over `block_rows` (equal-length, column-homologous):
    pinned junction codons; degenerate codons on `deg_cols` (covering the block's AA
    set there); fixed codons elsewhere (must be constant across the block). Returns
    (oligo, variants, deg_bases, nt) or None if un-coverable or an expansion (with
    junction context) could contain a forbidden site."""
    width = len(block_rows[0])
    codons, variants, deg = [], 1, 0
    for r in range(width):
        if r < len(lead):
            codons.append(lead[r])
            continue
        if r >= width - len(tail):
            codons.append(tail[r - (width - len(tail))])
            continue
        aas = frozenset(row[r] for row in block_rows) if r in deg_cols \
            else frozenset(block_rows[0][r])
        bc = best_codon_for(aas)
        if bc is None:
            return None
        codon, covered, ndeg, _ = bc
        codons.append(codon)
        variants *= len(covered)
        deg += ndeg
    oligo = "".join(codons)
    if _degenerate_has_forbidden(lctx + oligo + rctx):
        return None
    return oligo, variants, deg, len(oligo)


def densify_select(present, pieces, weights, cuts, tokens, L, chemistry,
                   arm_codons, max_junk_frac, max_junk, link_tau):
    """DeCoDe-lite: choose, per fragment gap-group, one of three encodings under the
    junk-fraction cutoff -- DISCRETE (kept pieces, one oligo each), HAP-SPLIT (one
    degenerate oligo per linked-column haplotype, degenerating only INDEPENDENT
    columns), or WHOLE-GROUP (a single degenerate oligo over all variable columns).
    The greedy adopts the change that adds the most newly-covered natural sequence
    (then fewest oligos, least junk). HAP-SPLIT/WHOLE-GROUP present every natural
    piece of the group, so they can cover cores discrete selection dropped; keeping
    linked columns discrete (HAP-SPLIT) avoids off-haplotype phantoms. Subsumes #1
    (WHOLE-GROUP is exactly its compression). Returns (frags, frag_units, library,
    U, W); U/W are the possibly-larger covered cores / natural weight."""
    K = len(present)
    n = len(pieces)
    bounds = [0] + cuts + [L]
    pins = (_pins_for_design(chemistry, cuts, tokens, arm_codons)
            if chemistry in ("gg", "hr") and cuts else [((), ())] * len(cuts))

    frag_groups = []                      # per f: list of group dicts
    group_of = [[None] * K for _ in range(n)]     # core -> (per f) group index
    for f in range(K):
        a, b = bounds[f], bounds[f + 1]
        lead, tail, lctx, rctx = _frag_pins(pins, f, K)
        pat_cores = defaultdict(list)
        for c in range(n):
            sl = _ALIGNED[c][a:b]
            pat_cores[tuple(i for i, ch in enumerate(sl) if ch != "-")].append(c)
        groups = []
        for pat, corelist in pat_cores.items():
            width = len(pat)
            piecemap = {c: pieces[c][f] for c in corelist}
            colres = [[piecemap[c][r] for c in corelist] for r in range(width)]
            cw = [weights[c] for c in corelist]
            pinned = set(range(len(lead))) | set(range(width - len(tail), width))
            variable = [r for r in range(width)
                        if len(set(colres[r])) > 1 and r not in pinned]
            indep = set(_independent_cols(variable, colres, cw, link_tau)) & set(variable)
            linked = [r for r in variable if r not in indep]
            rows = sorted(set(piecemap.values()))
            g = {"rows": rows, "kept": set(r for r in rows if r in present[f]),
                 "indep": sorted(indep), "variable": variable}
            # option DISCRETE
            g["opts"] = {0: {"oligos": len(g["kept"]), "variants": len(g["kept"]),
                             "units": [("discrete", r) for r in sorted(g["kept"])]}}
            # option WHOLE-GROUP (== #1 compression): degenerate every variable column
            if variable:
                o2 = _build_group_oligo(rows, set(variable), lead, tail, lctx, rctx)
                if o2:
                    g["opts"][2] = {"oligos": 1, "variants": o2[1],
                                    "units": [("degenerate", tuple(rows), o2[0], o2[2], o2[3])]}
            # option HAP-SPLIT: one oligo per linked-haplotype, degenerate indep cols
            blocks = defaultdict(list)
            for r in rows:
                blocks[tuple(r[c] for c in linked)].append(r)
            units, var, ok = [], 0, True
            for block in blocks.values():
                bo = _build_group_oligo(block, indep, lead, tail, lctx, rctx)
                if bo is None:
                    ok = False
                    break
                units.append(("degenerate", tuple(block), bo[0], bo[2], bo[3]))
                var += bo[1]
            if ok and units and (2 not in g["opts"] or len(units) < g["opts"][2]["oligos"]
                                 or var < g["opts"][2]["variants"]):
                g["opts"][1] = {"oligos": len(units), "variants": var, "units": units}
            groups.append(g)
            for c in corelist:
                group_of[c][f] = len(groups) - 1
        frag_groups.append(groups)

    def library_of(mode):
        lib = 1
        for f in range(K):
            lib *= max(1, sum(frag_groups[f][gi]["opts"][mode[f][gi]]["variants"]
                              for gi in range(len(frag_groups[f]))))
        return lib

    def oligos_of(mode):
        return sum(frag_groups[f][gi]["opts"][mode[f][gi]]["oligos"]
                   for f in range(K) for gi in range(len(frag_groups[f])))

    def covered_of(mode):
        cnt = wsum = 0
        for c in range(n):
            ok = True
            for f in range(K):
                gi = group_of[c][f]
                g = frag_groups[f][gi]
                if mode[f][gi] != 0:            # degenerate options present all rows
                    continue
                if pieces[c][f] not in g["kept"]:
                    ok = False
                    break
            if ok:
                cnt += 1
                wsum += weights[c]
        return cnt, wsum

    mode = [[0] * len(frag_groups[f]) for f in range(K)]
    cur_cnt, cur_W = covered_of(mode)
    cur_olig = oligos_of(mode)
    while True:
        best = None
        for f in range(K):
            for gi, g in enumerate(frag_groups[f]):
                for m in g["opts"]:
                    if m == mode[f][gi]:
                        continue
                    prev = mode[f][gi]
                    mode[f][gi] = m
                    lib = library_of(mode)
                    cnt, W = covered_of(mode)
                    olig = oligos_of(mode)
                    mode[f][gi] = prev
                    if lib > max_junk:
                        continue
                    frac = (lib - cnt) / lib if lib > 0 else 0.0
                    if frac > max_junk_frac:
                        continue
                    dW, dO = W - cur_W, olig - cur_olig
                    if dW < 0 or (dW == 0 and dO >= 0):
                        continue           # must add coverage or shrink the order
                    key = (-dW, dO, lib)
                    if best is None or key < best[0]:
                        best = (key, f, gi, m)
        if best is None:
            break
        _, f, gi, m = best
        mode[f][gi] = m
        cur_cnt, cur_W = covered_of(mode)
        cur_olig = oligos_of(mode)

    frag_units, variants = [[] for _ in range(K)], [0] * K
    for f in range(K):
        for gi, g in enumerate(frag_groups[f]):
            opt = g["opts"][mode[f][gi]]
            frag_units[f].extend(opt["units"])
            variants[f] += opt["variants"]

    frags, library = [], 1
    for f in range(K):
        library *= max(1, variants[f])
        nt = sum((len(u[1]) * 3 if u[0] == "discrete" else u[4]) for u in frag_units[f])
        frags.append({
            "a": bounds[f], "b": bounds[f + 1],
            "n_oligos": len(frag_units[f]),
            "natural_pieces": len(present[f]),
            "variants": variants[f],
            "n_pieces": len(present[f]),
            "deg_bases": sum(u[3] for u in frag_units[f] if u[0] == "degenerate"),
            "nt": nt,
            "lengths": sorted({len(p) for p in present[f]}),
            "degenerate_cells": sum(1 for u in frag_units[f] if u[0] == "degenerate"),
            "indep_cols": sorted({r for gi, g in enumerate(frag_groups[f])
                                  if mode[f][gi] == 1 for r in g["indep"]}),
        })
    U, W = covered_of(mode)
    return frags, frag_units, library, U, W


# --------------------------------------------------------------------------- #
# 6. Fragment encoding: discrete pieces + JUNK-BUDGETED degenerate compression.
#
# Degeneracy is a second lever alongside fragmentation. A degenerate codon at a
# variable column lets ONE oligo stand in for several natural pieces, but it emits
# the full cartesian product of its degenerate positions, so it can add phantom
# (junk) variants. We treat it as a budgeted trade: starting from the discrete
# design, greedily accept degenerate group-compressions -- biggest oligo saving
# first -- while the WHOLE design's library stays under the SAME --max-junk-pct
# cutoff (spending only the junk head-room the cutoff leaves). This does not buy
# new coverage (a later step); it shrinks the oligo order by spending head-room.
# Off unless --degenerate.
#
# Correctness: a degenerate oligo may only merge pieces that are column-homologous,
# i.e. share the SAME gap pattern within the fragment (identical length + identical
# alignment columns). Pieces with different indel patterns stay discrete.
# --------------------------------------------------------------------------- #

def _lib_cap(U, max_junk_frac, max_junk):
    """Largest library whose junk fraction (library - U)/library stays <= cutoff."""
    if max_junk_frac >= 1.0:
        return max_junk
    return min(max_junk, int(U / (1.0 - max_junk_frac)))


def _fragment_groups(a, b, kept_set):
    """Kept pieces of fragment [a,b) grouped by gap pattern, so each group is
    column-homologous and equal-length (safely degenerable as one oligo). Each
    distinct kept piece is assigned to exactly one group."""
    groups = defaultdict(list)
    seen = set()
    for s in _ALIGNED:
        sl = s[a:b]
        pc = sl.replace("-", "")
        if pc in kept_set and pc not in seen:
            seen.add(pc)
            pattern = tuple(i for i, ch in enumerate(sl) if ch != "-")
            groups[pattern].append(pc)
    return list(groups.values())


def _frag_pins(pins, f, K):
    """(lead, tail, left_ctx, right_ctx) codons for fragment f: its own pinned
    first/last junction codons, plus the neighbouring fragments' pinned codons as
    fixed context (so cross-junction windows are covered)."""
    lead = pins[f - 1][1] if f > 0 else ()
    tail = pins[f][0] if f < K - 1 else ()
    left_ctx = "".join(pins[f - 1][0]) if f > 0 else ""
    right_ctx = "".join(pins[f][1]) if f < K - 1 else ""
    return lead, tail, left_ctx, right_ctx


def _degenerate_group(rows, lead=(), tail=(), left_ctx="", right_ctx=""):
    """One degenerate IUPAC oligo covering all `rows` (equal-length, column-
    homologous). The first len(lead) / last len(tail) codons are PINNED to the
    junction codons; the rest use the cheapest codon covering that column's
    residues. Returns (oligo, n_aa_variants, deg_bases, nt) or None if a column is
    un-coverable or some expansion (with junction context) could contain a
    forbidden Type IIS site. n_aa_variants = distinct AA pieces the oligo yields
    (>= len(rows); the excess is phantom junk)."""
    width = len(rows[0])
    if width < len(lead) + len(tail):
        return None
    codons, variants, deg = [], 1, 0
    for j in range(width):
        if j < len(lead):
            codons.append(lead[j])
            continue
        if j >= width - len(tail):
            codons.append(tail[j - (width - len(tail))])
            continue
        req = frozenset(r[j] for r in rows)
        bc = best_codon_for(req)
        if bc is None:
            return None
        deg_codon, aas, ndeg, _ = bc
        codons.append(deg_codon)
        variants *= len(aas)
        deg += ndeg
    oligo = "".join(codons)
    if _degenerate_has_forbidden(left_ctx + oligo + right_ctx):
        return None
    return oligo, variants, deg, len(oligo)


def encode_design(present, U, cuts, tokens, L, chemistry, arm_codons,
                  max_junk_frac, max_junk, degenerate):
    """Turn the selected pieces into orderable oligos. Discrete = one oligo per
    kept piece. With `degenerate`, greedily compress column-homologous groups into
    degenerate oligos while the whole design's library stays under the junk cutoff.

    Returns (frags, frag_units, library):
      frags[f]      = summary dict (a,b,n_oligos,natural_pieces,variants,deg_bases,nt,lengths)
      frag_units[f] = list of ('discrete', aa) | ('degenerate', rows, oligo, deg, nt)
      library       = product over fragments of the presented variant count."""
    K = len(present)
    bounds = [0] + cuts + [L]
    pins = (_pins_for_design(chemistry, cuts, tokens, arm_codons)
            if chemistry in ("gg", "hr") and cuts else [((), ())] * len(cuts))

    frag_units = [[("discrete", p) for p in sorted(present[f])] for f in range(K)]
    variants = [len(present[f]) for f in range(K)]

    if degenerate and _HAVE_GREEDY:
        cap = _lib_cap(U, max_junk_frac, max_junk)
        library = 1
        for v in variants:
            library *= v
        cands = []
        for f in range(K):
            lead, tail, lctx, rctx = _frag_pins(pins, f, K)
            for rows in _fragment_groups(bounds[f], bounds[f + 1], set(present[f])):
                if len(rows) < 2:
                    continue
                opt = _degenerate_group(rows, lead, tail, lctx, rctx)
                if opt is None:
                    continue
                oligo, var, deg, nt = opt
                cands.append({"f": f, "rows": rows, "oligo": oligo, "deg": deg,
                              "nt": nt, "saved": len(rows) - 1,
                              "added": max(0, var - len(rows))})
        # spend head-room on the biggest oligo savings first (least junk breaks ties)
        cands.sort(key=lambda c: (-c["saved"], c["added"]))
        for c in cands:
            f = c["f"]
            new_lib = library // variants[f] * (variants[f] + c["added"])
            if new_lib > cap:
                continue
            rowset = set(c["rows"])
            frag_units[f] = [u for u in frag_units[f]
                             if not (u[0] == "discrete" and u[1] in rowset)]
            frag_units[f].append(("degenerate", tuple(c["rows"]), c["oligo"],
                                  c["deg"], c["nt"]))
            variants[f] += c["added"]
            library = new_lib

    frags, library = [], 1
    for f in range(K):
        library *= variants[f]
        nt = sum((len(u[1]) * 3 if u[0] == "discrete" else u[4]) for u in frag_units[f])
        frags.append({
            "a": bounds[f], "b": bounds[f + 1],
            "n_oligos": len(frag_units[f]),
            "natural_pieces": len(present[f]),
            "variants": variants[f],
            "n_pieces": len(present[f]),
            "deg_bases": sum(u[3] for u in frag_units[f] if u[0] == "degenerate"),
            "nt": nt,
            "lengths": sorted({len(p) for p in present[f]}),
        })
    return frags, frag_units, library


# --------------------------------------------------------------------------- #
# 6b. Back-translation / domestication.
#
# Turn the recommended AA design into concrete DNA and GUARANTEE that no forbidden
# Type IIS site (FORBIDDEN_SITES) occurs in ANY producible full-length sequence.
#
# For gg/hr the junction residues are constant (gg: the two overhang flanks; hr:
# the homology arm), so their codons are PINNED once and shared by every piece in
# the adjacent layers. Each remaining residue of each kept piece is back-translated
# by a DP that picks synonymous codons avoiding the site, padded on both sides by
# the neighbouring pinned codons. Because those boundary codons are identical for
# all pieces in a layer, EVERY cartesian-product assembly of domesticated pieces is
# itself site-free -- so we verify the whole (possibly huge) library by checking
# each piece and each pairwise junction context once. Agnostic makes no assembly
# promise, so we merely prove each kept full-length core is site-free as one ORF.
# --------------------------------------------------------------------------- #

def back_translate(aa, lead=(), tail=(), left_ctx="", right_ctx=""):
    """Concrete codon string for residues `aa`, with the first len(lead) and last
    len(tail) codons PINNED to the supplied codons, choosing synonymous codons for
    the rest so that left_ctx + <codons> + right_ctx has no forbidden site.
    Returns the DNA, or None if no site-free assignment exists.

    DP over codon choices; state = trailing 5 nt (enough to catch a 6-nt site that
    straddles the next codon). One predecessor kept per state for reconstruction."""
    n = len(aa)
    lead, tail = list(lead), list(tail)
    layers = [{left_ctx[-5:]: (None, None)}]     # layer 0 = fixed left-context suffix
    for idx, res in enumerate(aa):
        if idx < len(lead):
            cands = [lead[idx]]
        elif idx >= n - len(tail):
            cands = [tail[idx - (n - len(tail))]]
        else:
            cands = _CODON_BY_AA.get(res, [])
        cur = {}
        for suf in layers[-1]:
            for cod in cands:
                s = suf + cod
                if _has_forbidden_site(s):
                    continue
                ns = s[-5:]
                if ns not in cur:
                    cur[ns] = (suf, cod)
        if not cur:
            return None
        layers.append(cur)
    end = next((suf for suf in layers[-1] if not _has_forbidden_site(suf + right_ctx)), None)
    if end is None:
        return None
    codons, suf = [], end
    for layer in range(len(layers) - 1, 0, -1):
        prev, cod = layers[layer][suf]
        codons.append(cod)
        suf = prev
    codons.reverse()
    return "".join(codons)


def _pins_for_design(chemistry, cuts, tokens, arm_codons):
    """Per internal junction, the concrete codons each adjacent fragment owns:
    (left_codons, right_codons) = the last codons of the left fragment and the
    first codons of the right fragment, shared by all pieces at that junction."""
    pins = []
    for cut, tok in zip(cuts, tokens):
        if chemistry == "gg":
            lc, rc = _gg_pin_codons(cut, tok)
            pins.append(((lc,), (rc,)))
        elif chemistry == "hr":
            half = arm_codons // 2
            arm_aa = "".join(_const_aa(j) for j in range(cut - half, cut + half))
            arm_dna = back_translate(arm_aa)          # domesticate the arm once
            if arm_dna is None:
                raise ValueError(f"cannot domesticate homology arm at col {cut}")
            cods = [arm_dna[i:i + 3] for i in range(0, len(arm_dna), 3)]
            pins.append((tuple(cods[:half]), tuple(cods[half:])))
        else:
            pins.append(((), ()))
    return pins


def domesticate(rec, L, chemistry, arm_codons):
    """Attach concrete DNA/oligos to the recommended design and verify site-freedom.
    Discrete units are back-translated to concrete DNA; degenerate units keep their
    IUPAC oligo but are screened so NO expansion contains the site. Sets
    rec['backtranslation_ok'], rec['frag_oligos'] (gg/hr), and rec['bt_*']. Returns
    bool. Guarantee holds over the whole library because pinned junction codons make
    each piece + each pairwise junction context sufficient to check."""
    if chemistry not in ("gg", "hr"):
        return _domesticate_agnostic(rec, L)

    cuts, tokens = rec["cuts"], rec["tokens"]
    units = rec["frag_units"]
    K = len(units)
    pins = _pins_for_design(chemistry, cuts, tokens, arm_codons) if cuts else []

    frag_oligos = []
    for f in range(K):
        lead, tail, lctx, rctx = _frag_pins(pins, f, K)
        outs = []
        for u in units[f]:
            if u[0] == "discrete":
                dna = back_translate(u[1], lead=lead, tail=tail,
                                     left_ctx=lctx, right_ctx=rctx)
                if dna is None or _has_forbidden_site(dna):
                    rec["backtranslation_ok"] = False
                    rec["backtranslation_fail"] = f"fragment {f + 1} discrete piece"
                    return False
                outs.append({"type": "discrete", "aa": u[1], "seq": dna})
            else:
                rows, oligo = u[1], u[2]
                if _degenerate_has_forbidden(lctx + oligo + rctx):
                    rec["backtranslation_ok"] = False
                    rec["backtranslation_fail"] = f"fragment {f + 1} degenerate oligo"
                    return False
                outs.append({"type": "degenerate", "rows": list(rows), "seq": oligo})
        frag_oligos.append(outs)

    # Verify every cross-junction pair (IUPAC-aware, so it covers discrete+degenerate
    # mixes). Redundant with the pinning proof, but an explicit check.
    pairs = 0
    for f in range(K - 1):
        for lo in frag_oligos[f]:
            for ro in frag_oligos[f + 1]:
                pairs += 1
                if _degenerate_has_forbidden(lo["seq"][-5:] + ro["seq"][:5]):
                    rec["backtranslation_ok"] = False
                    rec["backtranslation_fail"] = f"site across junction {f + 1}/{f + 2}"
                    return False

    rec["frag_oligos"] = frag_oligos
    rec["backtranslation_ok"] = True
    rec["bt_n_oligos"] = sum(len(o) for o in frag_oligos)
    rec["bt_junction_pairs"] = pairs
    rec["bt_examples"] = _assemble_examples(rec, L, chemistry, arm_codons)
    return True


def _assemble_examples(rec, L, chemistry, arm_codons, k=3):
    """A few encoded cores rendered as assembled full-length DNA (for the report),
    by back-translating each core's ACTUAL pieces with the junction pins -- concrete
    and site-free even where the design uses a degenerate oligo for that layer."""
    cuts, tokens, units = rec["cuts"], rec["tokens"], rec["frag_units"]
    K = len(units)
    bounds = [0] + cuts + [L]
    pins = _pins_for_design(chemistry, cuts, tokens, arm_codons) if cuts else []
    # a core is encoded iff each of its pieces is presented by that fragment's units
    presented = []
    for f in range(K):
        s = set()
        for u in units[f]:
            s.add(u[1]) if u[0] == "discrete" else s.update(u[1])
        presented.append(s)
    out, seen = [], set()
    for s in _ALIGNED:
        parts = [piece(s, bounds[i], bounds[i + 1]) for i in range(K)]
        if not all(parts[i] in presented[i] for i in range(K)):
            continue
        dnas, ok = [], True
        for f in range(K):
            lead, tail, lctx, rctx = _frag_pins(pins, f, K) if pins else ((), (), "", "")
            d = back_translate(parts[f], lead=lead, tail=tail, left_ctx=lctx, right_ctx=rctx)
            if d is None:
                ok = False
                break
            dnas.append(d)
        if ok and "".join(dnas) not in seen:
            seen.add("".join(dnas))
            out.append("".join(dnas))
        if len(out) >= k:
            break
    return out


def _domesticate_agnostic(rec, L):
    """Agnostic makes no assembly promise; just prove each kept full-length core is
    site-free as one continuous ORF (single-sequence back-translation)."""
    cuts, kept = rec["cuts"], rec["kept_pieces"]
    K = len(kept)
    bounds = [0] + cuts + [L]
    cores = set()
    for s in _ALIGNED:
        parts = tuple(piece(s, bounds[i], bounds[i + 1]) for i in range(K))
        if all(parts[i] in kept[i] for i in range(K)):
            cores.add("".join(parts))
    examples, ok = [], True
    for aa in sorted(cores):
        dna = back_translate(aa)
        if dna is None:
            ok = False
            break
        if len(examples) < 3:
            examples.append(dna)
    rec["backtranslation_ok"] = ok
    rec["bt_examples"] = examples
    rec["bt_n_oligos"] = sum(len(k_) for k_ in kept)
    rec["bt_junction_pairs"] = 0
    return ok


# --------------------------------------------------------------------------- #
# 7. Evaluate one K end-to-end.
# --------------------------------------------------------------------------- #

def evaluate_K(seqs, K, min_block, chemistry, arm_codons, L, max_junk_frac,
               max_junk, degenerate, reserved=frozenset(), densify=False, link_tau=0.2):
    placed = place_cuts(L, K, min_block, chemistry, arm_codons, reserved)
    if placed is None:
        return None                            # no chemistry-valid segmentation for this K
    cuts, tokens = placed
    pieces = core_pieces(cuts, L)
    weights = _WEIGHTS
    present, U, W, traj = marginal_select(pieces, weights, max_junk_frac, max_junk)
    disc = {"covered_cores": U, "covered_weight": W,
            "library": _libsize(present), "oligos": _n_oligos(present)}

    if densify:
        frags, frag_units, library, U, W = densify_select(
            present, pieces, weights, cuts, tokens, L, chemistry, arm_codons,
            max_junk_frac, max_junk, link_tau)
    else:
        frags, frag_units, library = encode_design(
            present, U, cuts, tokens, L, chemistry, arm_codons,
            max_junk_frac, max_junk, degenerate)

    n_cores = len(seqs)
    total_w = sum(weights)
    return {
        "K": K,
        "cuts": cuts,
        "tokens": tokens,               # per-junction overhangs (gg) / arms (hr) / None
        "kept_pieces": [sorted(present[f]) for f in range(K)],  # discrete kept pieces
        "frag_units": frag_units,       # per-fragment discrete/degenerate oligo plan
        "frags": frags,
        "discrete": disc,               # baseline before densify (for report delta)
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
        pcs = (f"{f['natural_pieces']} natural pieces -> {f['variants']} library variants"
               if f["variants"] != f["natural_pieces"] else f"{f['natural_pieces']} distinct pieces")
        lines.append(f"  fragment {i}: cols [{f['a']},{f['b']})  {pcs}, length {lvar}")
        lines.append(f"      -> {f['n_oligos']} oligos, {f['deg_bases']} degenerate nt, ~{f['nt']} nt")
    if args.densify:
        d = rec["discrete"]
        deg_cells = sum(f.get("degenerate_cells", 0) for f in rec["frags"])
        lines.append("")
        lines.append(f"  DENSIFY (linkage-guided degenerate coverage, "
                     f"--link-tau {args.link_tau}):")
        lines.append(f"    discrete baseline : {d['covered_cores']}/{n_cores} cores, "
                     f"{d['covered_weight']}/{total_w} seqs, {d['library']:,} library, "
                     f"{d['oligos']} oligos")
        lines.append(f"    after densify     : {rec['n_cores_encoded']}/{n_cores} cores, "
                     f"{rec['encoded_weight']}/{total_w} seqs, {rec['library_size']:,} library "
                     f"({rec['junk_pct']:.1f}% junk), {rec['total_oligos']} oligos")
        lines.append(f"    => +{rec['encoded_weight'] - d['covered_weight']} seqs covered, "
                     f"{d['oligos'] - rec['total_oligos']} fewer oligos, via {deg_cells} "
                     f"degenerate cell(s) (independent columns only; linked columns kept discrete)")
    elif args.degenerate:
        disc_oligos = sum(f["natural_pieces"] for f in rec["frags"])
        disc_lib = 1
        for f in rec["frags"]:
            disc_lib *= f["natural_pieces"]
        deg_bases = sum(f["deg_bases"] for f in rec["frags"])
        lines.append("")
        lines.append(f"  DEGENERACY (junk-budgeted): {disc_oligos} discrete oligos -> "
                     f"{rec['total_oligos']} oligos ({disc_oligos - rec['total_oligos']} saved, "
                     f"{deg_bases} degenerate nt); library {disc_lib:,} -> {rec['library_size']:,} "
                     f"(spent junk head-room under the {args.max_junk_pct:.0f}% cutoff)")
    lines.append("")
    # Chemistry-validated junctions (empty for K=1 / agnostic).
    if args.chemistry != "agnostic" and rec["cuts"]:
        if args.chemistry == "gg":
            lines.append("JUNCTIONS (Golden Gate) -- validated overhangs, "
                         "Level-1 high-fidelity + Level-2 mutually orthogonal:")
            for cut, oh in zip(rec["cuts"], rec["tokens"]):
                lines.append(f"  col {cut:>4}:  overhang 5'-{oh}-3'  (rc {_revcomp(oh)})")
            bb = "  |  ".join(sorted(BACKBONE_OVERHANGS))
            if args.shared_backbone_overhangs:
                lines.append(f"  [backbone overhangs {bb} SHARED with these junctions "
                             f"(first/last fragment mates the backbone and is excluded)]")
            else:
                lines.append(f"  [backbone overhangs {bb} reserved -- excluded from all "
                             f"internal junctions above]")
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
    # Domestication / back-translation: concrete DNA + Type IIS site guarantee.
    sites = "  ".join(sorted(FORBIDDEN_SITES))
    lines.append(f"DOMESTICATION ({args.gg_enzyme.upper()} site {sites} forbidden in "
                 f"every full-length sequence):")
    if rec.get("backtranslation_ok"):
        if args.chemistry in ("gg", "hr"):
            lines.append(f"  PASS -- rendered {rec['bt_n_oligos']} oligos (discrete pieces "
                         f"back-translated, degenerate oligos screened over all expansions) "
                         f"with pinned junction codons; verified {rec['bt_n_oligos']} oligos and "
                         f"{rec['bt_junction_pairs']} junction contexts site-free => ALL "
                         f"{rec['library_size']:,} producible sequences are site-free by "
                         f"construction.")
        else:
            lines.append(f"  PASS -- every kept full-length core back-translates to a "
                         f"site-free ORF ({rec['bt_n_oligos']} pieces). NOTE: agnostic "
                         f"gives no per-fragment reuse guarantee across junctions.")
        for i, dna in enumerate(rec.get("bt_examples", []), 1):
            shown = f"{dna[:60]}...{dna[-30:]}" if len(dna) > 90 else dna
            lines.append(f"    example full-length DNA {i} ({len(dna)} nt): {shown}")
    else:
        reason = rec.get("backtranslation_fail", "no site-free codon assignment")
        lines.append(f"  FAILED -- {reason}. Try --min-block-cols / different cuts "
                     f"or chemistry.")
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
        "domestication": {
            "forbidden_sites": sorted(FORBIDDEN_SITES),
            "backbone_overhangs_reserved": (not args.shared_backbone_overhangs
                                            and args.chemistry == "gg"),
            "backtranslation_ok": rec.get("backtranslation_ok"),
            "verified_oligos": rec.get("bt_n_oligos"),
            "verified_junction_pairs": rec.get("bt_junction_pairs"),
            "example_full_length_dna": rec.get("bt_examples", []),
        },
    }
    with open(os.path.join(run_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    # one FASTA per fragment: ONLY the oligos the design actually orders. Discrete
    # oligos are concrete domesticated DNA; degenerate oligos are IUPAC (a mix that
    # expands to several pieces). When domestication ran (gg/hr) we use frag_oligos;
    # otherwise fall back to the kept amino-acid pieces.
    frag_oligos = rec.get("frag_oligos")
    if frag_oligos is not None:
        for i, outs in enumerate(frag_oligos):
            with open(os.path.join(run_dir, f"fragment{i+1}.fasta"), "w") as fh:
                for k, o in enumerate(outs, 1):
                    if o["type"] == "discrete":
                        fh.write(f">frag{i+1}_p{k}_len{len(o['aa'])}aa_{len(o['seq'])}nt\n"
                                 f"{o['seq']}\n")
                    else:
                        fh.write(f">frag{i+1}_deg{k}_covers{len(o['rows'])}pieces_"
                                 f"{len(o['seq'])}nt_IUPAC\n{o['seq']}\n")
    else:
        for i, kept in enumerate(rec["kept_pieces"]):
            with open(os.path.join(run_dir, f"fragment{i+1}.fasta"), "w") as fh:
                for k, aa in enumerate(sorted(kept), 1):
                    fh.write(f">frag{i+1}_p{k}_len{len(aa)}aa\n{aa}\n")
    # assembled example full-length sequences as DNA (domestication-verified)
    if rec.get("bt_examples"):
        with open(os.path.join(run_dir, "examples_full_length_dna.fasta"), "w") as fh:
            for k, dna in enumerate(rec["bt_examples"], 1):
                fh.write(f">example{k}_len{len(dna)}nt\n{dna}\n")
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
                    help="enable junk-BUDGETED degenerate-codon oligo compression "
                         "(spends leftover junk head-room under --max-junk-pct to "
                         "cut the oligo count; does not add coverage)")
    ap.add_argument("--densify", action="store_true",
                    help="enable linkage-guided degenerate COVERAGE moves (DeCoDe-lite): "
                         "degenerate independent columns to cover cores discrete "
                         "selection dropped, keeping linked columns discrete. Supersedes "
                         "--degenerate. Needs greedy_oligo importable.")
    ap.add_argument("--link-tau", type=float, default=0.2,
                    help="APC-MI threshold (bits) below which a column counts as "
                         "independent and may be degenerated by --densify (default 0.2)")
    ap.add_argument("--gg-enzyme", choices=["bsmbi", "esp3i", "bsai"], default="bsmbi",
                    help="Type IIS enzyme whose recognition site is domesticated out "
                         "of every full-length sequence (default bsmbi/esp3i = CGTCTC; "
                         "bsai = GGTCTC). Its reverse complement is forbidden too.")
    ap.add_argument("--shared-backbone-overhangs", action="store_true",
                    help="GG only: permit the reserved backbone overhangs (CGGA/GGTG) "
                         "at internal junctions (shared-overhang / minimal-plasmid "
                         "designs, first/last fragment excluded). Default forbids them.")
    ap.add_argument("--out-dir", default="algoruns",
                    help="parent folder for per-run output subfolders")
    args = ap.parse_args()

    # Set the forbidden Type IIS recognition site (+ reverse complement) from the
    # chosen enzyme, and the reserved backbone overhangs for Golden Gate.
    global FORBIDDEN_SITES
    site = {"bsmbi": "CGTCTC", "esp3i": "CGTCTC", "bsai": "GGTCTC"}[args.gg_enzyme]
    FORBIDDEN_SITES = frozenset({site, _revcomp(site)})
    reserved = (BACKBONE_OVERHANGS if args.chemistry == "gg"
                and not args.shared_backbone_overhangs else frozenset())
    if args.densify and not _HAVE_GREEDY:
        sys.exit("--densify needs greedy_oligo.py importable (run from the script's dir).")

    seqs = read_aligned_cores(args.aln_fasta)
    prepare(seqs)
    L = len(seqs[0][0])

    max_junk_frac = args.max_junk_pct / 100.0
    results = []
    for K in range(1, args.k_max + 1):
        if K > 1 and K * args.min_block_cols > L:
            break                       # can't fit this many fragments
        r = evaluate_K(seqs, K, args.min_block_cols, args.chemistry, args.arm_codons,
                       L, max_junk_frac, args.max_junk, args.degenerate, reserved,
                       args.densify, args.link_tau)
        if r is not None:
            results.append(r)
    if not results:
        sys.exit("no chemistry-valid segmentation for any K (try smaller "
                 "--min-block-cols, a different --chemistry, or --arm-codons).")

    rec = recommend(results)
    # Back-translate the recommended design to concrete DNA and prove no forbidden
    # Type IIS site occurs in any producible full-length sequence.
    domesticate(rec, L, args.chemistry, args.arm_codons)
    stem = os.path.splitext(os.path.basename(args.aln_fasta))[0]
    run_dir, report = save_run(args.out_dir, stem, args, results, rec, L, seqs)
    print(report)
    print(f"\nsaved run to {run_dir}/  (report.txt, summary.json, fragment*.fasta)")


if __name__ == "__main__":
    main()
