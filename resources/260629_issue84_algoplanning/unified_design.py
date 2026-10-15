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
  python unified_design.py ../../ninetypidorfs/cluster2.core.aln.fasta
  python unified_design.py ../../ninetypidorfs/cluster2.core.aln.fasta --chemistry gg
  python unified_design.py ... --max-junk-pct 50 --k-max 6
"""

from __future__ import annotations

import argparse
import json
import math
import os
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
                out.append((b1 + b2 + b3, frozenset(aas), ndeg, trip))
    # cheapest = fewest amino acids produced (least junk), then fewest degenerate
    # bases (Twist charges for them), then fewest concrete triplets.
    out.sort(key=lambda r: (len(r[1]), r[2], r[3]))
    return out


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
    self-ligate), GC content 25-75% (extreme GC ligates poorly), and not one of
    the destination vector's reserved overhangs.  Potapov 2018 / Pryor 2020."""
    if oh in reserved or revcomp(oh) in reserved:
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


def distinct_count(aligned, a, b, _memo={}):
    key = (id(aligned), a, b)
    hit = _memo.get(key)
    if hit is None:
        hit = len({piece_of(s, a, b) for s in aligned})
        _memo[key] = hit
    return hit


def place_cuts(aligned, L, K, min_block, const, chemistry, arm_codons, reserved,
               node_budget=400_000):
    """Best K-1 cuts for this K, or None.  'agnostic' has no set-level constraint
    so a DP suffices; gg/hr must keep the chosen junctions mutually orthogonal --
    a property of the whole SET, which does not decompose over edges -- so they
    use a bounded DFS that carries the chosen tokens along the path (this is
    GGAssembler's 'rainbow shortest path', done exactly because K and the number
    of legal sites are both small)."""
    if K == 1:
        return [], []

    sites = []
    for p in range(min_block, L - min_block + 1):
        toks = junction_options(p, const, chemistry, arm_codons, reserved)
        if toks:
            sites.append((p, toks))
    if len(sites) < K - 1:
        return None

    best = {"cost": math.inf, "cuts": None, "tokens": None}
    nodes = [0]

    def dfs(idx, last, cuts, tokens, cost):
        nodes[0] += 1
        if nodes[0] > node_budget:
            return
        if len(cuts) == K - 1:
            if L - last < min_block:
                return
            total = cost + math.log(distinct_count(aligned, last, L))
            if total < best["cost"]:
                best.update(cost=total, cuts=list(cuts), tokens=list(tokens))
            return
        if cost >= best["cost"]:
            return                                # log-costs are >= 0, so prune
        for ci in range(idx, len(sites)):
            p, toks = sites[ci]
            if p - last < min_block:
                continue
            if L - p < (K - len(cuts) - 1) * min_block:
                continue
            newcost = cost + math.log(distinct_count(aligned, last, p))
            if newcost >= best["cost"]:
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
    if best["cuts"] is None:
        return None
    return best["cuts"], best["tokens"]


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
           widen_candidates):
    """Seed from the heaviest real core, then repeatedly add the natural core with
    the cheapest junk per newly-encoded natural sequence, while the design's junk
    fraction stays <= the cap.

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

    while True:
        base_junk = lib - len(covered)
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
            if new_lib > 0 and junk / new_lib > max_junk_frac:
                continue
            score = (junk - base_junk) / gain
            nt = sum(p[1] for p in plan)
            key = (score, nt, -gain, c)
            if best is None or key < best[0]:
                best = (key, c, plan)

        if best is None:
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
            "cuts": list(cuts), "tokens": list(tokens), "K": K}


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


def evaluate_K(aligned, weights, K, min_block, const, chemistry, arm_codons,
               reserved, L, max_junk_frac, widen_candidates):
    placed = place_cuts(aligned, L, K, min_block, const, chemistry, arm_codons,
                        reserved)
    if placed is None:
        return None
    cuts, tokens = placed
    d = greedy(aligned, weights, cuts, tokens, L, chemistry, max_junk_frac,
               widen_candidates)
    if d is None:
        return None
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
    })
    return d


def recommend(results):
    """Factor #1 of the spec measures the order in INFORMATION, not in number of
    sequences, so the frontier is ranked by natural sequences encoded per
    kilobase ordered.  Ties break toward more coverage, then fewer fragments."""
    return max(results, key=lambda r: (round(r["seqs_per_kb"], 4),
                                       r["encoded_weight"], -r["K"]))


def build_report(args, results, rec, aligned, weights, L, examples, bad):
    n, total_w = len(aligned), sum(weights)
    out = []
    out.append(f"input: {args.aln_fasta}")
    out.append(f"{n} unique cores ({total_w} natural sequences), alignment width {L}")
    out.append(f"chemistry: {args.chemistry}   junk cap: {args.max_junk_pct}% of library"
               f"   min block: {args.min_block_cols} cols   K: 1..{args.k_max}")
    out.append("")
    out.append("FRONTIER over number of fragments K")
    out.append(f"{'K':>4} {'cores':>10} {'nat seqs':>12} {'library':>12} {'junk%':>8}"
               f" {'oligos':>8} {'nt':>9} {'seq/kb':>8}")
    for r in results:
        mark = "  <== recommended" if r is rec else ""
        out.append(f"{r['K']:>4} {r['n_cores_encoded']:>4}/{n:<5} "
                   f"{r['encoded_weight']:>5}/{total_w:<6} {r['library']:>12,} "
                   f"{r['junk_pct']:>7.1f}% {r['oligos']:>8} {r['nt']:>9,} "
                   f"{r['seqs_per_kb']:>8.2f}{mark}")
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
    out.append(f"  order: {rec['oligos']} oligos, {rec['nt']:,} nt, "
               f"{rec['degenerate_bases']} degenerate bases")
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
        "recommended_K": rec["K"],
        "cuts": rec["cuts"],
        "junctions": [list(t) for t in rec["tokens"]],
        "frontier": [{k: r[k] for k in
                      ("K", "n_cores_encoded", "encoded_weight", "library", "junk",
                       "junk_pct", "oligos", "nt", "degenerate_bases",
                       "seqs_per_oligo", "seqs_per_kb", "coverage_pct")}
                     for r in results],
        "fragments": [[{"oligo": u.oligo(), "variants": u.variants, "nt": u.nt,
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
    ap.add_argument("--max-junk-pct", type=float, default=80.0,
                    help="junk cap: max %% of the produced library that may be junk")
    ap.add_argument("--k-max", type=int, default=6, help="max fragments to try")
    ap.add_argument("--min-block-cols", type=int, default=20,
                    help="minimum fragment width in alignment columns")
    ap.add_argument("--arm-codons", type=int, default=6,
                    help="hr only: constant residues each side of a cut")
    ap.add_argument("--widen-candidates", type=int, default=3,
                    help="how many nearest units to consider widening per layer")
    ap.add_argument("--gg-enzyme", choices=["bsmbi", "esp3i", "bsai"], default="bsmbi",
                    help="Type IIS enzyme whose site is banned everywhere")
    ap.add_argument("--shared-backbone-overhangs", action="store_true",
                    help="gg: allow the reserved backbone overhangs internally")
    ap.add_argument("--out-dir", default="algoruns")
    args = ap.parse_args()

    global FORBIDDEN_SITES
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

    results = []
    for K in range(1, args.k_max + 1):
        if K > 1 and K * args.min_block_cols > L:
            break
        r = evaluate_K(aligned, weights, K, args.min_block_cols, const,
                       args.chemistry, args.arm_codons, reserved, L,
                       max_junk_frac, args.widen_candidates)
        if r is not None:
            results.append(r)
    if not results:
        sys.exit("no valid segmentation for any K (try --min-block-cols, "
                 "--chemistry agnostic, or a smaller --arm-codons)")

    rec = recommend(results)
    examples, bad = assemble_examples(rec, aligned, L)
    report = build_report(args, results, rec, aligned, weights, L, examples, bad)
    stem = os.path.splitext(os.path.basename(args.aln_fasta))[0]
    run = save_run(args.out_dir, stem, args, results, rec, aligned, L, report, examples)
    print(report)
    print(f"\nsaved run to {run}/")


if __name__ == "__main__":
    main()
