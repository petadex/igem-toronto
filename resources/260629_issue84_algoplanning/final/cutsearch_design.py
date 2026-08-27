#!/usr/bin/env python3
"""
Unified degenerate-codon / fragment oligo design.  Issue #84, rewrite.

Successor to marginal_design.py.  The change is architectural, not cosmetic: the
old pipeline chose discrete pieces first and only then tried to compress them
with degenerate codons, so the two levers were never priced against each other.
Here they are the SAME lever.

  KEY SIMPLIFICATION -- there is no separate "discrete" code path.
  A layer's contents are UNITS.  A unit is a list of per-column amino-acid sets.
    * a literal natural piece      = a unit whose sets are all singletons
    * a degenerate oligo           = a unit with some sets of size > 1
  So the two moves the search can make are:
    * ADD  a new unit   (buy another oligo)      -> +1 variant, +3*len nt
    * WIDEN a unit's sets (widen a codon)        -> +variants,   +0 nt
  and they are scored in the same two currencies, junk and nucleotides ordered.
  "Prefer the degenerate encoding when the junk is equal" is therefore not a rule
  we assert -- it falls out of the tie-break.

----------------------------------------------------------------------------- #
THE THREE STAGES
----------------------------------------------------------------------------- #
STAGE 1  CONTEXT       Precompute everything needed to PRICE a move in O(1):
                       the degenerate-codon table (all 15^3 IUPAC triplets ->
                       amino-acid set, cost, stop-freeness), which alignment
                       columns are constant, and where a chemistry-valid
                       junction could physically go.  Decides nothing.

STAGE 2  SEARCH        For each candidate segmentation, ONE greedy: seed from the
                       heaviest real core, then repeatedly add the natural core
                       with the cheapest junk per newly-encoded natural sequence
                       -- taking, per layer, whichever of {reuse, widen, add} is
                       cheapest.  Stop when no core can be added without breaching
                       the junk cap.  Works on amino-acid sets; DNA enters only as
                       a sliding-window check that a widened codon cannot spell a
                       forbidden Type IIS site.

STAGE 3  MATERIALIZE   Turn the winning symbolic design into concrete DNA: pin the
                       junction codons, emit the IUPAC oligos, assemble example
                       full-length sequences and prove no forbidden site occurs,
                       then report coverage / junk / nucleotides ordered.

----------------------------------------------------------------------------- #
DEFINITIONS (as written up in the notebook)
----------------------------------------------------------------------------- #
  L_T  target library of unique cores, weighted by w(s) = how many natural ORFs
       collapsed onto that core (the `_n<k>` FASTA suffix).
  L_O  protein library the ordered DNA produces = the cartesian product across
       fragments of each layer's producible pieces.
  covered   core s is covered iff every one of its K pieces is producible.
  coverage  weighted:   sum w(s) over covered / sum w(s) over all.
  counting junk  UNweighted:  1 - |L_T & L_O| / |L_O|.   This is the cap.
  (|L_O| is computed as the product of per-layer piece counts, an upper bound --
   two piece-tuples can concatenate to the same protein when an indel shifts
   material across a cut -- so reported junk is conservative.)

Run:
  # cores must be ALIGNED first (gaps '-'), e.g. in WSL:
  #   mafft --auto cluster2.core.fasta > cluster2.core.aln.fasta
  python cutsearch_design.py ../stage0_cluster1/out/clusters/1/c1.core.aln.fasta
  python cutsearch_design.py <cores.aln.fasta> --chemistry gg
  python cutsearch_design.py <cores.aln.fasta> --chemistry gg --k-max 5 \n      --cut-candidates 50 --proxy-candidates 5 --seed 7 \n      --max-library 300 --max-nt 25000
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import os
import random
import sys
from collections import defaultdict
from datetime import datetime

# =========================================================================== #
# STAGE 1a -- the degenerate-codon table.
#
# Built once.  Everything the search needs to price a codon move is a lookup in
# here: given the amino acids a column must produce, what is the cheapest
# stop-free IUPAC triplet that covers them, and what does it ACTUALLY produce
# (which may be a superset -- e.g. {D,K} is only reachable via RAN = {D,E,N,K},
# and those two extra amino acids are junk we are forced to buy).
# =========================================================================== #

_BASES = "TCAG"
_AA_STRING = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
AA_BY_CODON = dict(zip([a + b + c for a in _BASES for b in _BASES for c in _BASES],
                       _AA_STRING))

CODONS_BY_AA = defaultdict(list)
for _c, _a in AA_BY_CODON.items():
    if _a != "*":
        CODONS_BY_AA[_a].append(_c)

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}



# Relative synonymous codon usage, E. coli K-12 (Kazusa).  Used ONLY as a
# tie-break: the primary keys (amino acids produced, then degenerate bases) are
# untouched, so coverage, junk and library size are unaffected -- this changes
# only WHICH synonymous codon spells a choice that was already made.
#
# Without it the sort tied for every synonymous concrete codon and fell through
# to alphabetical order, which selected GGA (Gly), ATA (Ile), CTA (Leu) and AGA
# (Arg) -- the four rarest E. coli codons, the ones Rosetta/RIL strains exist to
# supplement.  That put 24-29% of positions on rare codons against a 5-10%
# baseline for native genes, which is an expression problem, not a design one.
ECOLI_USAGE = {
    "GCT":.16,"GCC":.27,"GCA":.21,"GCG":.36,
    "CGT":.38,"CGC":.40,"CGA":.06,"CGG":.10,"AGA":.04,"AGG":.02,
    "AAT":.45,"AAC":.55, "GAT":.63,"GAC":.37, "TGT":.45,"TGC":.55,
    "CAA":.35,"CAG":.65, "GAA":.69,"GAG":.31,
    "GGT":.34,"GGC":.40,"GGA":.11,"GGG":.15, "CAT":.57,"CAC":.43,
    "ATT":.51,"ATC":.42,"ATA":.07,
    "TTA":.13,"TTG":.13,"CTT":.10,"CTC":.10,"CTA":.04,"CTG":.50,
    "AAA":.77,"AAG":.23, "ATG":1.0, "TTT":.57,"TTC":.43,
    "CCT":.16,"CCC":.12,"CCA":.19,"CCG":.53,
    "TCT":.15,"TCC":.15,"TCA":.12,"TCG":.15,"AGT":.15,"AGC":.28,
    "ACT":.17,"ACC":.44,"ACA":.13,"ACG":.27, "TGG":1.0,
    "TAT":.57,"TAC":.43, "GTT":.26,"GTC":.22,"GTA":.15,"GTG":.37,
}
USE_CODON_USAGE = True
USAGE_OK = 0.16


def _build_degenerate_table():
    """All 15^3 IUPAC triplets -> (triplet, amino-acid set, #degenerate bases,
    #concrete triplets), stop-free only, cheapest first."""
    out = []
    for b1 in IUPAC:
        for b2 in IUPAC:
            for b3 in IUPAC:
                aas, stop = set(), False
                for n1 in IUPAC[b1]:
                    for n2 in IUPAC[b2]:
                        for n3 in IUPAC[b3]:
                            aa = AA_BY_CODON[n1 + n2 + n3]
                            if aa == "*":
                                stop = True
                            else:
                                aas.add(aa)
                if stop:
                    continue                     # never buy a truncation
                trip = len(IUPAC[b1]) * len(IUPAC[b2]) * len(IUPAC[b3])
                ndeg = sum(1 for b in (b1, b2, b3) if len(IUPAC[b]) > 1)
                # WORST usage over the expansions: a degenerate codon forces
                # every triplet it spans, so the weakest one is what limits
                # translation.
                exp3 = [n1 + n2 + n3 for n1 in IUPAC[b1] for n2 in IUPAC[b2]
                        for n3 in IUPAC[b3] if AA_BY_CODON[n1 + n2 + n3] != "*"]
                worst = min((ECOLI_USAGE.get(c, 0.0) for c in exp3), default=0.0)
                # GC deviation from 50%, averaged over the expansions.
                gcdev = (abs(sum(sum(1 for x in c if x in "GC") / 3.0
                                 for c in exp3) / len(exp3) - 0.5)
                         if exp3 else 0.5)
                out.append((b1 + b2 + b3, frozenset(aas), ndeg, trip, worst, gcdev))
    # cheapest = fewest amino acids produced (least junk), then fewest degenerate
    # bases (Twist charges for them), then fewest concrete triplets, then -- only
    # as a tie-break -- the best-translated synonymous option.
    # Usage is CLAMPED at USAGE_OK: the goal is to avoid codons E. coli
    # translates badly, not to maximise usage.  Maximising it picked the
    # GC-richest wobble every time and pushed 50-nt windows to 84% GC, which is
    # a synthesis-failure risk -- one problem swapped for another.  Clamping
    # makes every adequately-used codon equal, and GC balance then decides.
    # gcdev is ROUNDED before comparison: |1/3 - 0.5| and |2/3 - 0.5| are equal
    # in exact arithmetic but differ in the last bits as floats, so an unrounded
    # key let floating-point noise pick between e.g. GAA (used 0.69) and GAG
    # (0.31).  With the tie made genuine, TRUE usage -- unclamped -- is the final
    # arbiter, so among codons that are equally safe and equally GC-balanced we
    # take the one E. coli actually prefers.
    out.sort(key=lambda r: (len(r[1]), r[2], r[3],
                            -min(r[4], USAGE_OK) if USE_CODON_USAGE else 0.0,
                            round(r[5], 6) if USE_CODON_USAGE else 0.0,
                            -r[4] if USE_CODON_USAGE else 0.0))
    return [(t, aas, nd, ntr) for t, aas, nd, ntr, _u, _g in out]


DEG_TABLE = _build_degenerate_table()
_codon_cache: dict[frozenset, list] = {}


def codons_for(aa_set):
    """Every stop-free IUPAC triplet covering `aa_set`, cheapest first.
    Returns [(triplet, produced_aa_set, n_degenerate_bases), ...]."""
    key = frozenset(aa_set)
    hit = _codon_cache.get(key)
    if hit is None:
        hit = [(t, aas, nd) for t, aas, nd, _ in DEG_TABLE if key <= aas]
        _codon_cache[key] = hit
    return hit


def best_codon(aa_set):
    """Cheapest stop-free triplet covering `aa_set`, or None."""
    c = codons_for(aa_set)
    return c[0] if c else None


# =========================================================================== #
# STAGE 1b -- forbidden Type IIS sites.
#
# The recognition site of the assembly enzyme must not occur ANYWHERE in a
# produced sequence or the enzyme cuts the insert internally.  Because our oligos
# carry IUPAC ambiguity, the test is "could SOME expansion spell the site", which
# is a sliding window over the ambiguity sets.  A 6-nt site spans at most three
# codons, so this is a local check -- which is why it can run inside the search.
# =========================================================================== #

FORBIDDEN_SITES = frozenset({"CGTCTC", "GAGACG"})   # BsmBI/Esp3I, set from --gg-enzyme
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def iupac_may_contain_site(seq):
    """True if some concrete expansion of an IUPAC string contains a forbidden
    site.  `seq` may include plain ACGT (which are singleton IUPAC codes)."""
    for site in FORBIDDEN_SITES:
        n = len(site)
        for i in range(len(seq) - n + 1):
            if all(site[k] in IUPAC[seq[i + k]] for k in range(n)):
                return True
    return False


# =========================================================================== #
# 2. Input
# =========================================================================== #

def read_aligned_cores(path):
    """[(aligned_seq, weight)].  `>coreN_n<k>` -- k natural ORFs collapsed here."""
    seqs, header, buf = [], None, []

    def flush():
        if header is not None and buf:
            w = 1
            if "_n" in header:
                tail = header.rsplit("_n", 1)[1]
                if tail.isdigit():
                    w = int(tail)
            seqs.append(("".join(buf), w))

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
    if not seqs:
        sys.exit(f"no sequences read from {path}")
    lengths = {len(s) for s, _ in seqs}
    if len(lengths) != 1:
        sys.exit(f"input is not aligned: {len(lengths)} distinct lengths {sorted(lengths)}")
    return seqs


# =========================================================================== #
# STAGE 1c -- where a junction could physically go.
#
# Level 1 (per site): can a junction live at this boundary at all?
#   gg -- the 4-nt overhang straddles the cut (last 2 nt of the left residue's
#         codon + first 2 nt of the right residue's codon), so BOTH flanking
#         columns must be constant across the cluster; the overhang must be
#         non-palindromic and GC-balanced, and must avoid the backbone's own
#         overhangs.  Synonymous codons give a SET of achievable overhangs.
#   hr -- a window of >= arm_codons constant residues must straddle the cut.
# Level 2 (per SET) is not a property of one boundary -- two junctions must not
# cross-react -- so it is enforced during cut placement, not here.
# =========================================================================== #

BACKBONE_OVERHANGS = frozenset({"CGGA", "GGTG"})
GG_OH_LEN = 4


def constant_columns(aligned):
    """Columns where every core carries the same NON-GAP residue."""
    L = len(aligned[0])
    const = [None] * L
    for j in range(L):
        c = aligned[0][j]
        if c == "-":
            continue
        if all(s[j] == c for s in aligned):
            const[j] = c
    return const


def overhang_ok(oh, reserved):
    """Individually high-fidelity: not self-complementary (palindromes
    self-ligate), GC content 25-75% (extreme GC ligates poorly), and neither
    equal to NOR CROSS-REACTIVE WITH one of the destination vector's reserved
    overhangs.  Potapov 2018 / Pryor 2020.

    The cross-reactivity clause matters: rejecting only EXACT matches let an
    internal overhang sit one mismatch from CGGA, which is the threshold at
    which Golden Gate actually mis-ligates -- so it would have passed every
    check and still ligated into the backbone.  gg_conflict() is applied among
    the chosen internal overhangs during cut placement but never against the
    reserved pair, and this is where that hole is closed.  Because the reserved
    overhangs are FIXED constants, this is a per-site (Level 1) property and can
    be enforced here rather than in the set-level search."""
    if oh in reserved or revcomp(oh) in reserved:
        return False
    if any(gg_conflict(oh, r) for r in reserved):
        return False
    if oh == revcomp(oh):
        return False
    gc = sum(1 for b in oh if b in "GC")
    return 1 <= gc <= 3


def gg_conflict(a, b):
    """Two overhangs mis-ligate if identical, or if one anneals to the other's
    reverse complement within a single mismatch (the dominant GGA failure)."""
    if a == b:
        return True
    rb = revcomp(b)
    return sum(1 for x, y in zip(a, rb) if x != y) <= 1


def hr_conflict(a, b, max_ident=0.8):
    """Homology arms that are near-identical mis-recombine."""
    if len(a) != len(b):
        return False
    same = sum(1 for x, y in zip(a, b) if x == y)
    return same / len(a) > max_ident


def junction_options(p, const, chemistry, arm_codons, reserved):
    """The junction tokens available at boundary `p` (cut between columns p-1
    and p), or [] if no junction can live there.

    gg -> list of (overhang, left_codon, right_codon); the codons are PINNED into
          the flanking oligos so the overhang is actually realised.
    hr -> list of (arm_string,)   agnostic -> [(None,)]"""
    if chemistry == "agnostic":
        return [(None,)]
    if chemistry == "gg":
        la, ra = const[p - 1], const[p]
        if la is None or ra is None:
            return []
        seen, out = set(), []
        for cl in CODONS_BY_AA[la]:
            for cr in CODONS_BY_AA[ra]:
                oh = cl[1:] + cr[:2]
                if oh in seen or not overhang_ok(oh, reserved):
                    continue
                seen.add(oh)
                out.append((oh, cl, cr))
        return out
    if chemistry == "hr":
        lo, hi = p - arm_codons, p + arm_codons
        if lo < 0 or hi > len(const):
            return []
        if any(const[j] is None for j in range(lo, hi)):
            return []
        return [("".join(const[j] for j in range(lo, hi)),)]
    raise ValueError(chemistry)


def tokens_conflict(a, b, chemistry):
    if chemistry == "gg":
        return gg_conflict(a[0], b[0])
    if chemistry == "hr":
        return hr_conflict(a[0], b[0])
    return False


# =========================================================================== #
# STAGE 1d -- candidate segmentations.
#
# Cut placement is now only a CANDIDATE GENERATOR: it proposes segmentations and
# the stage-2 search is the arbiter, so cuts are finally judged by the coverage
# they actually achieve under the cap rather than by a proxy computed before we
# know what we will buy.  The proxy is still how candidates are RANKED: minimise
# the product of per-layer distinct piece counts (== the library at full
# coverage), which is additive in logs, i.e. a shortest path (RASPP/SwiftLib).
# =========================================================================== #

def piece_of(s, a, b):
    """The real (ungapped) subsequence a core contributes over columns [a, b)."""
    return s[a:b].replace("-", "")


_DC_MEMO = {}
_DC_OWNER = [None]


def distinct_count(aligned, a, b):
    """Memoised per alignment.  The cache is keyed on WHICH alignment object we
    are looking at and cleared when that changes -- an earlier version keyed on
    id(aligned) alone, and CPython recycles ids after garbage collection, so a
    driver processing several clusters in one process could have read another
    cluster's piece counts and ranked cuts against the wrong data.  Holding the
    owner reference also stops that id being reused while it is live."""
    if _DC_OWNER[0] is not aligned:
        _DC_MEMO.clear()
        _DC_OWNER[0] = aligned
    key = (a, b)
    hit = _DC_MEMO.get(key)
    if hit is None:
        hit = len({piece_of(s, a, b) for s in aligned})
        _DC_MEMO[key] = hit
    return hit


def layer_widths(cuts, L):
    b = [0] + list(cuts) + [L]
    return [b[i + 1] - b[i] for i in range(len(b) - 1)]


def place_cuts(aligned, L, K, min_block, const, chemistry, arm_codons, reserved,
               node_budget=2_000_000, n_keep=1, pool_factor=10,
               max_layer_cols=None):
    """The cheapest `n_keep`-ish K-1 cut sets for this K.

    Returns (candidates, truncated) with candidates = [(proxy_cost, cuts, tokens),
    ...] cheapest first, or ([], truncated) when no segmentation exists.

    CHANGED vs unified_design.py, which returned only the single cheapest set:
    the proxy ranked here (product of per-layer distinct piece counts) is a
    FULL-COVERAGE quantity that the design never actually builds, so handing the
    greedy one segmentation let a proxy decide something it cannot see.  We now
    keep a pool and let stage 2 arbitrate -- which is what this module always
    claimed to do.

    All chemistries use the same bounded DFS.  For gg/hr that is forced: the
    chosen junctions must be mutually orthogonal, a property of the whole SET
    which does not decompose over edges, so no shortest-path DP applies.  For
    'agnostic' a DP would suffice, but keeping one code path keeps the candidate
    pool and the truncation flag identical across chemistries.

    `truncated` is True when the node budget stopped the search: the pool is then
    whatever was reached, NOT provably the best.  The caller must surface it --
    at the old 400k budget this fired silently at K=4 on cluster 1 and returned
    a worse segmentation (cost 9.663) than an exhaustive search (9.512)."""
    if K == 1:
        if max_layer_cols is not None and L > max_layer_cols:
            return [], False
        return [(0.0, [], [])], False

    sites = []
    for p in range(min_block, L - min_block + 1):
        toks = junction_options(p, const, chemistry, arm_codons, reserved)
        if toks:
            sites.append((p, toks))
    if len(sites) < K - 1:
        return [], False

    pool_size = max(1, n_keep * pool_factor)
    keep = []                    # max-heap by cost: (-cost, serial, cuts, tokens)
    serial, nodes, truncated = [0], [0], [False]

    def bound():
        """Prune against the WORST candidate we are still keeping, not the best;
        with a pool of size P the P-th best cost is the admissible bound."""
        return math.inf if len(keep) < pool_size else -keep[0][0]

    def offer(cost, cuts, tokens):
        serial[0] += 1
        item = (-cost, serial[0], list(cuts), list(tokens))
        if len(keep) < pool_size:
            heapq.heappush(keep, item)
        elif cost < -keep[0][0]:
            heapq.heapreplace(keep, item)

    def dfs(idx, last, cuts, tokens, cost):
        nodes[0] += 1
        if nodes[0] > node_budget:
            truncated[0] = True
            return
        if len(cuts) == K - 1:
            w = L - last
            if w < min_block:
                return
            if max_layer_cols is not None and w > max_layer_cols:
                return
            offer(cost + math.log(distinct_count(aligned, last, L)), cuts, tokens)
            return
        if cost >= bound():
            return                                # log-costs are >= 0, so prune
        for ci in range(idx, len(sites)):
            p, toks = sites[ci]
            if p - last < min_block:
                continue
            if max_layer_cols is not None and p - last > max_layer_cols:
                break            # sites ascend in p, so every later p is wider too
            if L - p < (K - len(cuts) - 1) * min_block:
                continue
            newcost = cost + math.log(distinct_count(aligned, last, p))
            if newcost >= bound():
                continue
            for t in toks:
                if any(tokens_conflict(t, u, chemistry) for u in tokens):
                    continue
                cuts.append(p)
                tokens.append(t)
                dfs(ci + 1, p, cuts, tokens, newcost)
                cuts.pop()
                tokens.pop()
                if chemistry == "agnostic":
                    break                         # tokens are all None; one suffices

    dfs(0, 0, [], [], 0.0)
    cands = sorted(((-c, cu, tk) for c, _s, cu, tk in keep), key=lambda r: r[0])
    return cands, truncated[0]


def legal_sites(L, min_block, const, chemistry, arm_codons, reserved):
    out = []
    for p in range(min_block, L - min_block + 1):
        toks = junction_options(p, const, chemistry, arm_codons, reserved)
        if toks:
            out.append((p, toks))
    return out


def random_segmentations(aligned, L, K, min_block, const, chemistry, arm_codons,
                         reserved, n, rng, max_layer_cols=None, max_tries=50_000):
    """Uniformly sample valid cut sets, under EXACTLY the constraints the DFS
    obeys -- same min_block spacing, same per-site Level-1 filtering, same
    set-level orthogonality.  Random therefore cannot place two cuts closer than
    min_block; the arms differ only in how candidates are proposed, which is what
    makes the comparison between them fair.

    This is not a control that we keep for tidiness.  On cluster 1 the proxy pool
    collapses to 9 distinct designs at K=3 (4 at K=5) out of 120 candidates,
    while random gives ~119 -- so random is the arm that supplies real diversity,
    and past N~20 at K=3 it beats the proxy outright.  See CUT_SEARCH.md."""
    if K < 2:
        return []
    sites = legal_sites(L, min_block, const, chemistry, arm_codons, reserved)
    if len(sites) < K - 1:
        return []
    out, seen, tries = [], set(), 0
    while len(out) < n and tries < max_tries:
        tries += 1
        pick = sorted(rng.sample(range(len(sites)), K - 1))
        ps = tuple(sites[i][0] for i in pick)
        if ps in seen:
            continue
        b = [0] + list(ps) + [L]
        widths = [b[i + 1] - b[i] for i in range(len(b) - 1)]
        if min(widths) < min_block:
            continue
        if max_layer_cols is not None and max(widths) > max_layer_cols:
            continue
        toks, ok = [], True
        for i in pick:
            cand = [t for t in sites[i][1]
                    if not any(tokens_conflict(t, u, chemistry) for u in toks)]
            if not cand:
                ok = False
                break
            toks.append(rng.choice(cand))
        if not ok:
            continue
        seen.add(ps)
        cost = sum(math.log(distinct_count(aligned, b[i], b[i + 1]))
                   for i in range(len(b) - 1))
        out.append((cost, list(ps), toks))
    return out


def diversify(cands, n_keep, L, bucket_cols=10):
    """Pick `n_keep` candidates SPREAD over segmentation shape, not the n_keep
    cheapest -- those are near-duplicates (`[63,83]`, `[63,84]`, `[62,83]` ...)
    that all give the greedy the same design.  Bucket by widest layer, then take
    the cheapest unused candidate from each bucket in turn.

    Widest layer is the right axis because it decides the longest oligo, which is
    what synthesis actually constrains, and because the proxy is systematically
    biased toward one huge layer plus narrow ones."""
    if len(cands) <= n_keep:
        return cands
    buckets = defaultdict(list)
    for c in cands:
        buckets[max(layer_widths(c[1], L)) // bucket_cols].append(c)
    order = sorted(buckets)
    out, i = [], 0
    while len(out) < n_keep and any(buckets[k] for k in order):
        k = order[i % len(order)]
        if buckets[k]:
            out.append(buckets[k].pop(0))
        i += 1
    return sorted(out, key=lambda r: r[0])


# =========================================================================== #
# STAGE 2 -- the design state.
#
# A UNIT is one ordered oligo.  It covers a set of alignment columns (its gap
# pattern) and carries one amino-acid set per column.  Two cores can share a unit
# only if they have the SAME gap pattern in that layer, because a degenerate
# codon can encode "V or I" but cannot encode "residue or nothing".
# =========================================================================== #

class Unit:
    __slots__ = ("cols", "sets", "codons", "exp", "variants", "nt", "pins")

    def __init__(self, cols, sets, codons, exp, variants, nt, pins):
        self.cols = cols          # tuple of alignment column indices, ascending
        self.sets = sets          # list[frozenset] -- what we ASKED for
        self.codons = codons      # list[str]       -- IUPAC triplet per column
        self.exp = exp            # list[frozenset] -- what we actually GET
        self.variants = variants  # product of |exp| = this oligo's library share
        self.nt = nt              # 3 * len(cols)
        self.pins = pins          # {col: codon} forced by a junction overhang

    def oligo(self):
        return "".join(self.codons)

    def covers(self, cols, residues):
        if cols != self.cols:
            return False
        return all(r in e for r, e in zip(residues, self.exp))


def build_unit(cols, sets, pins):
    """Realise a unit: choose a codon per column, honour junction pins, and make
    sure no expansion of the resulting oligo can spell a forbidden Type IIS site.

    The site check is the ONE place DNA enters stage 2.  It is local (a 6-nt site
    spans <= 3 codons), and where it fails we repair it by swapping in the next
    cheapest synonymous/degenerate codon at a column inside the offending window
    -- only rejecting the unit if no repair exists."""
    codons, exp = [], []
    for col, aa_set in zip(cols, sets):
        if col in pins:
            cod = pins[col]
            produced = frozenset(AA_BY_CODON[cod])
            if not aa_set <= produced:
                return None                       # pin cannot express this column
            codons.append(cod)
            exp.append(produced)
            continue
        bc = best_codon(aa_set)
        if bc is None:
            return None                           # no stop-free codon covers it
        codons.append(bc[0])
        exp.append(bc[1])

    if iupac_may_contain_site("".join(codons)):
        if not _repair_sites(cols, sets, pins, codons, exp):
            return None

    variants = 1
    for e in exp:
        variants *= len(e)
    return Unit(tuple(cols), list(sets), codons, exp, variants, 3 * len(cols), pins)


def _repair_sites(cols, sets, pins, codons, exp):
    """Try alternative codons, one column at a time, until no expansion of the
    oligo can spell a forbidden site.  Only unpinned columns may be changed."""
    for _ in range(len(cols)):
        seq = "".join(codons)
        if not iupac_may_contain_site(seq):
            return True
        # find the first offending window and try to break it
        fixed = False
        for i in range(len(seq) - 5):
            window = seq[i:i + 6]
            if not iupac_may_contain_site(window):
                continue
            for col_i in range(i // 3, min(len(cols), (i + 6 + 2) // 3)):
                if cols[col_i] in pins:
                    continue
                for cand, produced, _nd in codons_for(sets[col_i])[1:6]:
                    old_c, old_e = codons[col_i], exp[col_i]
                    codons[col_i], exp[col_i] = cand, produced
                    if not iupac_may_contain_site("".join(codons)[max(0, i - 6):i + 12]):
                        fixed = True
                        break
                    codons[col_i], exp[col_i] = old_c, old_e
                if fixed:
                    break
            if fixed:
                break
        if not fixed:
            return False
    return not iupac_may_contain_site("".join(codons))


def layer_view(core, a, b):
    """(cols, residues) -- the columns this core actually occupies in [a, b) and
    the residues it puts there.  The cols tuple IS the gap pattern."""
    cols, res = [], []
    for j in range(a, b):
        if core[j] != "-":
            cols.append(j)
            res.append(core[j])
    return tuple(cols), tuple(res)


# =========================================================================== #
# STAGE 2 -- the unified greedy.
# =========================================================================== #

def cross_junction_site(layers):
    """True if some pair of adjacent-layer oligos could spell a forbidden site
    ACROSS their boundary.

    build_unit() only ever sees one unit at a time, so a site straddling a
    junction is invisible to it: the window is
    [codon p-2][pinned p-1][pinned p][codon p+1], and while the two pinned
    codons are identical in every oligo, the flanking ones vary.  Before this,
    the only cross-junction test in the program ran inside assemble_examples()
    -- after the design was final, on three examples, with no way to act on the
    result.  Checking it here lets evaluate_K discard the candidate and take
    another, which turns "we checked and it was clean" into "the search cannot
    return one".  A 6-nt site needs at most 5 nt either side of the boundary."""
    for f in range(len(layers) - 1):
        for a in layers[f]:
            oa = a.oligo()[-5:]
            for b in layers[f + 1]:
                if iupac_may_contain_site(oa + b.oligo()[:5]):
                    return True
    return False


def layer_total(units):
    return sum(u.variants for u in units)


def library_size(layers):
    lib = 1
    for units in layers:
        lib *= max(1, layer_total(units))
    return lib


def _covered_in_layer(units, view):
    for u in units:
        if u.covers(*view):
            return True
    return False


def move_options(units, view, pins, widen_candidates):
    """How this layer could accommodate a core, cheapest first.
    Returns [(delta_variants, delta_nt, kind, index, unit_or_None), ...].

      kind 'reuse' -- already producible, free
      kind 'widen' -- extend an existing unit's amino-acid sets (0 extra oligos)
      kind 'add'   -- buy another oligo (+1 variant, +3*len nt)

    Both paid moves are priced in the SAME currencies, which is the whole point:
    'widen' wins ties on nucleotides, so a degenerate encoding is preferred
    whenever it costs no more junk -- without that preference being hard-coded."""
    cols, residues = view
    if _covered_in_layer(units, view):
        return [(0, 0, "reuse", -1, None)]

    opts = []
    fresh = build_unit(cols, [frozenset(r) for r in residues], pins)
    if fresh is not None:
        opts.append((fresh.variants, fresh.nt, "add", -1, fresh))

    # Only try to widen the units that are closest to already covering this core;
    # widening a distant unit is never the cheapest move and costs time to price.
    ranked = []
    for i, u in enumerate(units):
        if u.cols != cols:
            continue                              # different gap pattern
        need = sum(1 for r, e in zip(residues, u.exp) if r not in e)
        ranked.append((need, i))
    ranked.sort()
    for _need, i in ranked[:widen_candidates]:
        u = units[i]
        new_sets = [s | {r} for s, r in zip(u.sets, residues)]
        w = build_unit(cols, new_sets, u.pins)
        if w is not None:
            opts.append((w.variants - u.variants, 0, "widen", i, w))

    opts.sort(key=lambda o: (o[0], o[1]))
    return opts


def greedy(aligned, weights, cuts, tokens, L, chemistry, max_junk_frac,
           widen_candidates, max_library=None, max_nt=None, oligo_overhead=0):
    """Seed from the heaviest real core, then repeatedly add the natural core that
    consumes the LEAST LIBRARY per newly-encoded natural sequence, while every cap
    still holds.

    OBJECTIVE (unchanged): maximise natural sequences encoded.  HEURISTIC (changed
    from unified_design.py): rank by `delta_library / gain` rather than
    `delta_junk / gain`.  Once the cap is on library size rather than junk
    fraction, the library is the resource that runs out, and a greedy has to price
    what it runs out of -- picking the move that leaves the most room for later
    ones.  Ranking by `gain` alone is NOT the same thing and is measurably worse
    (37 vs 52 natural sequences at K=5 under a 1000-library cap): different cores
    cost wildly different amounts of library for the same gain, and greedy-by-
    objective spends the budget on the expensive ones first.  `delta_junk / gain`
    and `delta_library / gain` gave identical designs in all ten configurations
    tested, so this is a relabelling in practice, not a behaviour change.

    NB the seed is a real core, never the column-wise consensus: the consensus is
    a chimera that covers nothing, so seeding there would start the library at
    100% junk and permanently widen every column where the majority residue is
    not one we end up needing."""
    n = len(aligned)
    bounds = [0] + list(cuts) + [L]
    K = len(bounds) - 1
    pins = _pins_by_column(cuts, tokens, chemistry)

    views = [[layer_view(aligned[i], bounds[f], bounds[f + 1]) for f in range(K)]
             for i in range(n)]

    layers = [[] for _ in range(K)]
    # cov[f] = set of cores layer f can already produce.  Only ever grows, and
    # only in the layer we touched, so coverage never needs a full recompute.
    cov = [set() for _ in range(K)]

    seed = max(range(n), key=lambda i: (weights[i], -i))
    for f in range(K):
        u = build_unit(views[seed][f][0], [frozenset(r) for r in views[seed][f][1]],
                       pins[f])
        if u is None:
            return None
        layers[f].append(u)
        _absorb(layers[f], views, f, cov[f], n)

    covered = _intersect(cov, n)
    W = sum(weights[i] for i in covered)
    lib = library_size(layers)
    traj = [_snapshot(layers, covered, W, lib)]
    stopped_by = "no improving move"

    while True:
        base_junk = lib - len(covered)
        cur_oligos = sum(len(u) for u in layers)
        cur_nt = (sum(u.nt for units in layers for u in units)
                  + oligo_overhead * cur_oligos)
        blocked = {"junk": 0, "library": 0, "nt": 0}
        best = None
        for c in range(n):
            if c in covered:
                continue
            plan, ok = [], True
            for f in range(K):
                opts = move_options(layers[f], views[c][f], pins[f], widen_candidates)
                if not opts:
                    ok = False
                    break
                plan.append(opts[0])
            if not ok:
                continue

            new_lib = 1
            for f in range(K):
                new_lib *= max(1, layer_total(layers[f]) + plan[f][0])
            new_cov = _hypothetical_coverage(layers, cov, plan, views, n, K)
            new_W = sum(weights[i] for i in new_cov)
            gain = new_W - W
            if gain <= 0:
                continue
            junk = new_lib - len(new_cov)
            nt = sum(p[1] for p in plan)
            new_oligos = cur_oligos + sum(1 for pf in plan if pf[2] == "add")
            new_nt = cur_nt + nt + oligo_overhead * (new_oligos - cur_oligos)

            # ---- feasibility: every cap is a FILTER, never part of the rank --
            if new_lib > 0 and junk / new_lib > max_junk_frac:
                blocked["junk"] += 1
                continue
            if max_library is not None and new_lib > max_library:
                blocked["library"] += 1
                continue
            if max_nt is not None and new_nt > max_nt:
                blocked["nt"] += 1
                continue

            # ---- rank: library slots consumed per natural sequence gained ---
            score = (new_lib - lib) / gain
            key = (score, nt, -gain, c)
            if best is None or key < best[0]:
                best = (key, c, plan)

        if best is None:
            # Why did we stop?  With several caps and only one of them in the
            # ranking, this line is what makes the output interpretable.
            stopped_by = (max(blocked, key=blocked.get)
                          if any(blocked.values()) else "no improving move")
            break
        _, c, plan = best
        for f in range(K):
            dv, _dn, kind, idx, unit = plan[f]
            if kind == "add":
                layers[f].append(unit)
            elif kind == "widen":
                layers[f][idx] = unit
            _absorb(layers[f], views, f, cov[f], n)
        covered = _intersect(cov, n)
        W = sum(weights[i] for i in covered)
        lib = library_size(layers)
        traj.append(_snapshot(layers, covered, W, lib))

    return {"layers": layers, "cov": cov, "covered": covered, "W": W,
            "library": lib, "trajectory": traj, "pins": pins,
            "cuts": list(cuts), "tokens": list(tokens), "K": K,
            "stopped_by": stopped_by,
            "nt_with_overhead": (sum(u.nt for units in layers for u in units)
                                 + oligo_overhead
                                 * sum(len(u) for u in layers))}


def _absorb(units, views, f, cov_f, n):
    """Recompute which cores this layer can produce.  Free recombinants and cores
    picked up incidentally by a widened codon are absorbed here at no cost."""
    for i in range(n):
        if i in cov_f:
            continue
        if _covered_in_layer(units, views[i][f]):
            cov_f.add(i)


def _hypothetical_coverage(layers, cov, plan, views, n, K):
    """Coverage if `plan` were applied, without mutating anything."""
    out = None
    for f in range(K):
        dv, _dn, kind, idx, unit = plan[f]
        cf = cov[f]
        if kind != "reuse":
            cf = set(cf)
            for i in range(n):
                if i in cf:
                    continue
                if unit.covers(*views[i][f]):
                    cf.add(i)
                elif kind == "widen":
                    # the replaced unit is gone, but every other unit remains
                    for j, u in enumerate(layers[f]):
                        if j != idx and u.covers(*views[i][f]):
                            cf.add(i)
                            break
        out = cf if out is None else (out & cf)
    return out if out is not None else set()


def _intersect(cov, n):
    out = set(range(n))
    for cf in cov:
        out &= cf
    return out


def _snapshot(layers, covered, W, lib):
    oligos = sum(len(u) for u in layers)
    return {"oligos": oligos, "library": lib, "junk": lib - len(covered),
            "junk_pct": 100.0 * (lib - len(covered)) / lib if lib else 0.0,
            "covered_cores": len(covered), "covered_weight": W,
            "nt": sum(u.nt for units in layers for u in units)}


def _pins_by_column(cuts, tokens, chemistry):
    """Junction overhangs are realised by PINNING the codons either side of the
    cut.  Returns one {column: codon} map per layer."""
    K = len(cuts) + 1
    pins = [dict() for _ in range(K)]
    if chemistry != "gg":
        return pins
    for i, p in enumerate(cuts):
        _oh, cl, cr = tokens[i]
        pins[i][p - 1] = cl
        pins[i + 1][p] = cr
    return pins


# =========================================================================== #
# STAGE 3 -- materialise and report.
# =========================================================================== #

def concrete_codon(iupac_cod, aa):
    """One concrete triplet inside `iupac_cod` that encodes `aa` (for assembling
    example full-length sequences)."""
    for n1 in IUPAC[iupac_cod[0]]:
        for n2 in IUPAC[iupac_cod[1]]:
            for n3 in IUPAC[iupac_cod[2]]:
                if AA_BY_CODON[n1 + n2 + n3] == aa:
                    return n1 + n2 + n3
    return None


def assemble_examples(design, aligned, L, k=3):
    """Concrete DNA for a few encoded cores -- and the proof obligation: none of
    them may contain a forbidden Type IIS site."""
    bounds = [0] + design["cuts"] + [L]
    K = design["K"]
    out, bad = [], 0
    for i in sorted(design["covered"])[:k]:
        dna, ok = [], True
        for f in range(K):
            cols, res = layer_view(aligned[i], bounds[f], bounds[f + 1])
            unit = next((u for u in design["layers"][f] if u.covers(cols, res)), None)
            if unit is None:
                ok = False
                break
            for cod, r in zip(unit.codons, res):
                cc = concrete_codon(cod, r)
                if cc is None:
                    ok = False
                    break
                dna.append(cc)
            if not ok:
                break
        if not ok:
            continue
        seq = "".join(dna)
        if any(s in seq for s in FORBIDDEN_SITES):
            bad += 1
        out.append((i, seq))
    return out, bad


def design_metrics(d, aligned, weights, L):
    total_w = sum(weights)
    oligos = sum(len(u) for u in d["layers"])
    nt = sum(u.nt for units in d["layers"] for u in units)
    return total_w, oligos, nt


def evaluate_K(aligned, weights, K, min_block, const, chemistry, arm_codons,
               reserved, L, max_junk_frac, widen_candidates,
               n_candidates=1, max_layer_cols=None, node_budget=2_000_000,
               max_library=None, max_nt=None, oligo_overhead=0,
               proxy_candidates=5, seed=0, exhaustive_max=150):
    """Run the greedy on EVERY candidate segmentation and keep the design that
    actually wins, rather than trusting the proxy's favourite.

    Selection is by encoded natural sequences, ties to fewer nucleotides.  That
    is deliberately NOT the proxy: the proxy ranks the full-coverage library,
    which under any cap is not the library we build."""
    # How to split the candidate budget between the two generators.  At K=2 the
    # number of segmentations equals the number of legal sites, so if that is
    # small we enumerate ALL of them and the question of "is N enough" does not
    # arise for that row.
    sites = legal_sites(L, min_block, const, chemistry, arm_codons, reserved)
    mode = "pool"
    if K == 1:
        n_proxy, n_random = 1, 0
        mode = "exhaustive"
    elif K == 2 and len(sites) <= exhaustive_max:
        n_proxy, n_random = len(sites), 0
        mode = "exhaustive"
    else:
        n_proxy = min(proxy_candidates, n_candidates)
        n_random = max(0, n_candidates - n_proxy)

    cands, truncated = place_cuts(aligned, L, K, min_block, const, chemistry,
                                  arm_codons, reserved, node_budget=node_budget,
                                  n_keep=n_proxy,
                                  max_layer_cols=max_layer_cols)
    if not cands and n_random == 0:
        return None
    n_pool = len(cands)
    if mode != "exhaustive":
        cands = diversify(cands, n_proxy, L)
    tagged = [("proxy", i, c) for i, c in enumerate(cands, start=1)]

    # Seeded, so a run is reproducible: the same seed reproduces the same design.
    if n_random:
        rng = random.Random(seed * 1000 + K)
        seen = {tuple(c[1]) for c in cands}
        rc = random_segmentations(aligned, L, K, min_block, const, chemistry,
                                  arm_codons, reserved, n_random, rng,
                                  max_layer_cols=max_layer_cols)
        rc = [c for c in rc if tuple(c[1]) not in seen]
        tagged += [("random", i, c) for i, c in enumerate(rc, start=1)]
    if not tagged:
        return None

    best = None
    n_site_rejected = 0
    for arm, arm_rank, (proxy_cost, cuts, tokens) in tagged:
        cd = greedy(aligned, weights, cuts, tokens, L, chemistry, max_junk_frac,
                    widen_candidates, max_library=max_library, max_nt=max_nt,
                    oligo_overhead=oligo_overhead)
        if cd is None:
            continue
        if cross_junction_site(cd["layers"]):
            n_site_rejected += 1
            continue
        nt_c = sum(u.nt for units in cd["layers"] for u in units)
        key = (-cd["W"], nt_c)             # most sequences, then fewest nt
        if best is None or key < best[0]:
            best = (key, cd, arm, arm_rank, proxy_cost)
    if best is None:
        return None
    _key, d, win_arm, win_rank, win_proxy = best
    d.update({
        "candidates_pooled": n_pool,
        "candidates_tried": len(tagged),
        "n_proxy": sum(1 for a, _, _ in tagged if a == "proxy"),
        "n_random": sum(1 for a, _, _ in tagged if a == "random"),
        "candidate_mode": mode,
        "candidates_rejected_junction_site": n_site_rejected,
        "seed": seed,
        "cut_search_truncated": truncated,
        "winner_arm": win_arm,
        "winner_arm_rank": win_rank,
        "winner_proxy_cost": win_proxy,
        "proxy_favourite_won": win_arm == "proxy" and win_rank == 1,
        "layer_widths": layer_widths(d["cuts"], L),
    })
    total_w = sum(weights)
    oligos = sum(len(u) for u in d["layers"])
    nt = sum(u.nt for units in d["layers"] for u in units)
    deg = sum(1 for units in d["layers"] for u in units
              for c in u.codons for b in c if len(IUPAC[b]) > 1)
    d.update({
        "n_cores_encoded": len(d["covered"]),
        "n_cores_total": len(aligned),
        "encoded_weight": d["W"],
        "total_weight": total_w,
        "coverage_pct": 100.0 * d["W"] / total_w if total_w else 0.0,
        "junk": d["library"] - len(d["covered"]),
        "junk_pct": (100.0 * (d["library"] - len(d["covered"])) / d["library"]
                     if d["library"] else 0.0),
        "oligos": oligos,
        "nt": nt,
        "degenerate_bases": deg,
        "seqs_per_oligo": d["W"] / max(1, oligos),
        "seqs_per_kb": 1000.0 * d["W"] / max(1, nt),
        "longest_oligo_nt": max(u.nt for units in d["layers"] for u in units),
        "stopped_by": d.get("stopped_by"),
        "nt_ordered": d.get("nt_with_overhead", nt),
    })
    return d


def recommend(results, ratio=False):
    """Under explicit caps the frontier is directly comparable, so the answer is
    simply the design encoding the most natural sequences, ties to fewer
    nucleotides ordered and then fewer fragments.

    This replaces ranking by seqs_per_kb, which is a RATIO and mis-fires whenever
    the arms it compares reach different coverage -- it repeatedly recommended a
    3%-coverage design over a full-coverage one.  `--rank-seqs-per-kb` restores
    the old behaviour for comparison."""
    if ratio:
        return max(results, key=lambda r: (round(r["seqs_per_kb"], 4),
                                           r["encoded_weight"], -r["K"]))
    return max(results, key=lambda r: (r["encoded_weight"],
                                       -r.get("nt_ordered", r["nt"]), -r["K"]))


def build_report(args, results, rec, aligned, weights, L, examples, bad,
                 overhead=0):
    n, total_w = len(aligned), sum(weights)
    out = []
    out.append(f"input: {args.aln_fasta}")
    out.append(f"{n} unique cores ({total_w} natural sequences), alignment width {L}")
    caps = []
    caps.append(f"library <= {args.max_library:,}" if args.max_library else "library: none")
    caps.append(f"nt <= {args.max_nt:,}" if args.max_nt else "nt: none")
    if args.max_junk_pct < 100:
        caps.append(f"junk <= {args.max_junk_pct}%")
    if args.max_oligo_nt:
        caps.append(f"oligo <= {args.max_oligo_nt} nt")
    out.append(f"chemistry: {args.chemistry}   caps: {'   '.join(caps)}")
    out.append(f"min block: {args.min_block_cols} cols   K: 1..{args.k_max}   "
               f"per-oligo assembly overhead: {overhead} nt")
    out.append(f"cut candidates: {args.cut_candidates}/K "
               f"({args.proxy_candidates} shortest-path + rest random)   "
               f"SEED {args.seed}")
    out.append("")
    out.append("FRONTIER over number of fragments K")
    out.append(f"{'K':>4} {'cores':>10} {'nat seqs':>12} {'library':>10} {'junk%':>7}"
               f" {'oligos':>7} {'nt ord':>8} {'longest':>8} {'stopped by':>14}"
               f" {'winner':>16}")
    for r in results:
        mark = "  <== recommended" if r is rec else ""
        trunc = "!" if r.get("cut_search_truncated") else " "
        out.append(f"{r['K']:>4} {r['n_cores_encoded']:>4}/{n:<5} "
                   f"{r['encoded_weight']:>5}/{total_w:<6} {r['library']:>10,} "
                   f"{r['junk_pct']:>6.1f}% {r['oligos']:>7} "
                   f"{r.get('nt_ordered', r['nt']):>8,} "
                   f"{r.get('longest_oligo_nt',0):>7}nt "
                   f"{str(r.get('stopped_by','?')):>14} "
                   f"{r.get('winner_arm','?')}#{r.get('winner_arm_rank','?')}"
                   f"/{r.get('candidates_tried',1)}{trunc}{mark}")
    nrej = sum(r.get("candidates_rejected_junction_site") or 0 for r in results)
    if nrej:
        out.append(f"  {nrej} candidate segmentation(s) were DISCARDED for a "
                   f"forbidden site spanning a fragment junction.")
    out.append("")
    if any(r.get("cut_search_truncated") for r in results):
        out.append("  ! = cut search hit the node budget; that row's segmentation "
                   "is the best FOUND, not provably the best.")
    nrand = sum(1 for r in results if r.get("winner_arm") == "random")
    out.append(f"  winner column = which generator produced the chosen "
               f"segmentation, its rank within that arm, and the total tried.")
    if nrand:
        out.append(f"  the random arm won {nrand}/{len(results)} rows -- "
                   f"evidence, from this run, that the shortest-path pool alone "
                   f"is not sufficient.")
    out.append("")
    out.append("=" * 70)
    out.append(f"RECOMMENDED: K = {rec['K']} fragment(s)")
    out.append("=" * 70)
    bounds = [0] + rec["cuts"] + [L]
    out.append("segments: " + "  |  ".join(f"[{bounds[i]},{bounds[i+1]})"
                                           for i in range(rec["K"])))
    out.append("")
    for f, units in enumerate(rec["layers"]):
        ndeg = sum(1 for u in units for c in u.codons for b in c if len(IUPAC[b]) > 1)
        widened = sum(1 for u in units if any(len(e) > 1 for e in u.exp))
        out.append(f"  fragment {f+1}: {len(units)} oligos "
                   f"({widened} carrying degenerate codons, {ndeg} degenerate bases), "
                   f"{layer_total(units)} producible pieces, "
                   f"{sum(u.nt for u in units):,} nt")
    out.append("")
    if rec["tokens"] and rec["tokens"][0][0] is not None:
        label = "overhang" if args.chemistry == "gg" else "homology arm"
        out.append(f"JUNCTIONS ({args.chemistry}) -- Level-1 valid, Level-2 mutually orthogonal:")
        for p, tok in zip(rec["cuts"], rec["tokens"]):
            if args.chemistry == "gg":
                out.append(f"  col {p:>4}:  {label} 5'-{tok[0]}-3'  (rc {revcomp(tok[0])})"
                           f"  pinned codons {tok[1]}|{tok[2]}")
            else:
                out.append(f"  col {p:>4}:  {label} {tok[0]}")
        if args.chemistry == "gg":
            out.append(f"  [backbone overhangs {' | '.join(sorted(BACKBONE_OVERHANGS))}"
                       f" reserved -- excluded from internal junctions]")
        out.append("")
    out.append("HEADLINE:")
    out.append(f"  {rec['n_cores_encoded']}/{n} unique cores "
               f"= {rec['encoded_weight']}/{total_w} natural sequences "
               f"({rec['coverage_pct']:.0f}%)")
    out.append(f"  library {rec['library']:,}   junk {rec['junk']:,} "
               f"({rec['junk_pct']:.1f}% of library, cap {args.max_junk_pct}%)")
    out.append(f"  order: {rec['oligos']} oligos, {rec['nt']:,} nt of coding "
               f"sequence, {rec.get('nt_ordered', rec['nt']):,} nt ordered "
               f"(incl. {overhead} nt/oligo assembly overhead), "
               f"{rec['degenerate_bases']} degenerate bases")
    out.append(f"  greedy stopped because: {rec.get('stopped_by','?')}")
    out.append(f"  segmentation from the {rec.get('winner_arm','?')} arm, "
               f"rank {rec.get('winner_arm_rank','?')} of "
               f"{rec.get('n_proxy',0)} proxy + {rec.get('n_random',0)} random "
               f"({rec.get('candidate_mode','pool')}), seed {args.seed}")
    out.append(f"  efficiency: {rec['seqs_per_oligo']:.2f} seqs/oligo, "
               f"{rec['seqs_per_kb']:.2f} seqs/kb")
    out.append("")
    out.append(f"forbidden-site check on {len(examples)} assembled full-length "
               f"example(s): {'FAIL' if bad else 'clean'}")
    out.append("")
    out.append("greedy trajectory (each row = one core bought):")
    out.append(f"{'oligos':>8} {'cores':>7} {'nat seqs':>9} {'library':>12} {'junk%':>8}")
    for t in rec["trajectory"]:
        out.append(f"{t['oligos']:>8} {t['covered_cores']:>7} {t['covered_weight']:>9} "
                   f"{t['library']:>12,} {t['junk_pct']:>7.1f}%")
    return "\n".join(out)


def save_run(out_root, stem, args, results, rec, aligned, L, report, examples):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run = os.path.join(out_root, f"{ts}_{stem}_{args.chemistry}_unified")
    os.makedirs(run, exist_ok=True)
    with open(os.path.join(run, "report.txt"), "w") as fh:
        fh.write(report + "\n")

    summary = {
        "input": os.path.abspath(args.aln_fasta),
        "args": vars(args),
        "seed": args.seed,
        "recommended_K": rec["K"],
        "cuts": rec["cuts"],
        "junctions": [list(t) for t in rec["tokens"]],
        "frontier": [{k: r.get(k) for k in
                      ("K", "n_cores_encoded", "encoded_weight", "library", "junk",
                       "junk_pct", "oligos", "nt", "degenerate_bases",
                       "seqs_per_oligo", "seqs_per_kb", "coverage_pct",
                       "longest_oligo_nt", "layer_widths", "candidates_pooled",
                       "candidates_tried", "cut_search_truncated",
                       "winner_arm", "winner_arm_rank", "n_proxy", "n_random",
                       "candidate_mode", "seed", "stopped_by", "nt_ordered",
                       "candidates_rejected_junction_site",
                       "proxy_favourite_won")}
                     for r in results],
        "fragments": [[{"oligo": u.oligo(), "variants": u.variants, "nt": u.nt,
                        "cols": list(u.cols),
                        "degenerate_columns": [i for i, e in enumerate(u.exp)
                                               if len(e) > 1],
                        "encodes": ["".join(sorted(e)) for e in u.exp]}
                       for u in units] for units in rec["layers"]],
    }
    with open(os.path.join(run, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    for f, units in enumerate(rec["layers"], start=1):
        with open(os.path.join(run, f"fragment{f}.fasta"), "w") as fh:
            for i, u in enumerate(units, start=1):
                fh.write(f">frag{f}_oligo{i}_var{u.variants}\n{u.oligo()}\n")
    if examples:
        with open(os.path.join(run, "examples_full_length_dna.fasta"), "w") as fh:
            for i, seq in examples:
                fh.write(f">core{i}\n{seq}\n")
    return run


# =========================================================================== #
# CLI
# =========================================================================== #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aln_fasta", help="ALIGNED core FASTA (gaps '-'), headers >coreN_n<k>")
    ap.add_argument("--chemistry", choices=["agnostic", "gg", "hr"], default="agnostic",
                    help="junction model (default agnostic = no assembly assumption)")
    ap.add_argument("--max-junk-pct", type=float, default=100.0,
                    help="junk cap, as %% of the produced library (default 100 = "
                         "OFF).  Redundant with --max-library: junk%% <= c is "
                         "exactly lib <= covered/(1-c), so this is a ratio "
                         "guarantee, not an independent budget")
    ap.add_argument("--max-library", type=int, default=None,
                    help="cap on |L_O|, the number of distinct proteins the order "
                         "can produce.  This is the wet-lab constraint and the "
                         "resource the greedy prices its moves against")
    ap.add_argument("--max-nt", type=int, default=None,
                    help="cap on total nucleotides ordered, INCLUDING per-oligo "
                         "assembly overhead; rations one cluster against the "
                         "whole multi-cluster budget")
    ap.add_argument("--oligo-overhead-nt", type=int, default=None,
                    help="nt added per oligo for Type IIS sites/spacers, counted "
                         "toward --max-nt (default 24 for gg, else 0)")
    ap.add_argument("--rank-seqs-per-kb", action="store_true",
                    help="rank the frontier by seqs/kb (the old ratio behaviour) "
                         "instead of most sequences within the caps")
    ap.add_argument("--k-max", type=int, default=6, help="max fragments to try")
    ap.add_argument("--min-block-cols", type=int, default=20,
                    help="minimum fragment width in alignment columns")
    ap.add_argument("--arm-codons", type=int, default=6,
                    help="hr only: constant residues each side of a cut")
    ap.add_argument("--widen-candidates", type=int, default=3,
                    help="how many nearest units to consider widening per layer")
    ap.add_argument("--cut-candidates", type=int, default=50,
                    help="TOTAL candidate segmentations per K to run the greedy "
                         "on, split between the two generators")
    ap.add_argument("--proxy-candidates", type=int, default=5,
                    help="how many of --cut-candidates come from the shortest-"
                         "path pool; the rest are random.  The pool collapses to "
                         "a handful of distinct designs, so a few saturate it")
    ap.add_argument("--seed", type=int, default=0,
                    help="seed for the random candidate arm.  RECORDED in the "
                         "report and summary.json -- the same seed reproduces "
                         "the same design")
    ap.add_argument("--exhaustive-max", type=int, default=150,
                    help="if a K has no more segmentations than this, enumerate "
                         "them all instead of sampling (applies at K<=2)")
    ap.add_argument("--max-oligo-nt", type=int, default=None,
                    help="reject any segmentation whose widest layer exceeds this "
                         "many nt (e.g. 300 for oligo pools); default unlimited")
    ap.add_argument("--cut-node-budget", type=int, default=2_000_000,
                    help="DFS node budget for cut placement; runs that hit it are "
                         "flagged '!' in the frontier")
    ap.add_argument("--no-codon-usage", action="store_true",
                    help="disable the E. coli codon-usage tie-break and fall "
                         "back to alphabetical order (reproduces pre-2026-08-26 "
                         "sequences; expect ~25%% rare codons)")
    ap.add_argument("--gg-enzyme", choices=["bsmbi", "esp3i", "bsai"], default="bsmbi",
                    help="Type IIS enzyme whose site is banned everywhere")
    ap.add_argument("--shared-backbone-overhangs", action="store_true",
                    help="gg: allow the reserved backbone overhangs internally")
    ap.add_argument("--out-dir", default="algoruns")
    args = ap.parse_args()

    global FORBIDDEN_SITES, USE_CODON_USAGE, DEG_TABLE
    if args.no_codon_usage:
        USE_CODON_USAGE = False
        DEG_TABLE = _build_degenerate_table()
        _codon_cache.clear()
    site = {"bsmbi": "CGTCTC", "esp3i": "CGTCTC", "bsai": "GGTCTC"}[args.gg_enzyme]
    FORBIDDEN_SITES = frozenset({site, revcomp(site)})
    reserved = (BACKBONE_OVERHANGS if args.chemistry == "gg"
                and not args.shared_backbone_overhangs else frozenset())

    seqs = read_aligned_cores(args.aln_fasta)
    aligned = [s for s, _ in seqs]
    weights = [w for _, w in seqs]
    L = len(aligned[0])
    const = constant_columns(aligned)
    max_junk_frac = args.max_junk_pct / 100.0
    overhead = (args.oligo_overhead_nt if args.oligo_overhead_nt is not None
                else (24 if args.chemistry == "gg" else 0))

    results = []
    for K in range(1, args.k_max + 1):
        if K > 1 and K * args.min_block_cols > L:
            break
        max_layer_cols = (args.max_oligo_nt // 3) if args.max_oligo_nt else None
        r = evaluate_K(aligned, weights, K, args.min_block_cols, const,
                       args.chemistry, args.arm_codons, reserved, L,
                       max_junk_frac, args.widen_candidates,
                       n_candidates=args.cut_candidates,
                       max_layer_cols=max_layer_cols,
                       node_budget=args.cut_node_budget,
                       max_library=args.max_library, max_nt=args.max_nt,
                       oligo_overhead=overhead,
                       proxy_candidates=args.proxy_candidates, seed=args.seed,
                       exhaustive_max=args.exhaustive_max)
        if r is not None:
            results.append(r)
    if not results:
        sys.exit("no valid segmentation for any K (try --min-block-cols, "
                 "--chemistry agnostic, or a smaller --arm-codons)")

    rec = recommend(results, ratio=args.rank_seqs_per_kb)
    examples, bad = assemble_examples(rec, aligned, L)
    report = build_report(args, results, rec, aligned, weights, L, examples, bad,
                          overhead=overhead)
    stem = os.path.splitext(os.path.basename(args.aln_fasta))[0]
    run = save_run(args.out_dir, stem, args, results, rec, aligned, L, report, examples)
    print(report)
    print(f"\nsaved run to {run}/")


if __name__ == "__main__":
    main()
