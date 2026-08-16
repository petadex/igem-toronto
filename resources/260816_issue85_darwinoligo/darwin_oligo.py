#!/usr/bin/env python3
"""
DARWIN downstream oligo designer.  Issue #85.

Given ONE base protein sequence and an explicit list of Top-K mutations from the
DARWIN algorithm, emit the Golden-Gate-ready DNA to order so that every
combination of those mutations is produced.

This is the SIMPLE sibling of `unified_design.py` (issue #84).  The difference is
that nothing is searched: we are told exactly which mutations to realise, so
there is no coverage objective and no junk cap.  The only questions are HOW to
realise each mutation and WHICH ones cannot be realised at all.

----------------------------------------------------------------------------- #
THE TWO OPERATIONS (same as the standard algorithm, used differently)
----------------------------------------------------------------------------- #
  SUBSTITUTION -> a DEGENERATE CODON covering {wild-type, mutant}.  Keeping the
      wild-type residue in the codon is what makes the library combinatorial:
      k substitutions give 2^k proteins, including the unmutated parent.

  INSERTION / DELETION -> a FRAGMENT.  A length change cannot be encoded in one
      oligo, so the fragment containing the indel is ordered TWICE (with and
      without it) and Golden Gate mixes the alternatives.  Each indel is
      isolated in its own fragment, which needs a valid junction between
      consecutive indels.

  FALLBACK: a few substitutions have no stop-free degenerate codon (e.g. {W,K}
      forces a TAG/TAA expansion).  Rather than dropping them, they are demoted
      to fragment variants -- the same 2-oligo treatment as an indel.

----------------------------------------------------------------------------- #
MUTATION NOTATION (1-based, '-' marks the gap side)
----------------------------------------------------------------------------- #
    A2N   substitution: base position 2 is A, also encode N
    -4R   insertion:    insert R immediately before base position 4
    G5-   deletion:     delete the G at base position 5

----------------------------------------------------------------------------- #
GOLDEN GATE IS ASSUMED (not an input)
----------------------------------------------------------------------------- #
  * the enzyme's recognition site (CGTCTC / GAGACG by default) must not occur
    ANYWHERE in any producible sequence -- enforced by synonymous codon choice,
    checked as a sliding window over the IUPAC oligos so it also holds for every
    expansion of a degenerate codon;
  * the destination backbone owns overhangs CGGA and GGTG, so internal junctions
    may not use them (unless --shared-backbone-overhangs);
  * a junction's 4-nt overhang must straddle two UNMUTATED residues, be
    non-palindromic and GC-balanced, and be mutually orthogonal to every other
    junction in the design.

Run:
  python darwin_oligo.py --seq-file base.fasta --mutations A2N -4R G5- W40K
  python darwin_oligo.py --seq MKV... --mutations-file muts.txt
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from itertools import product

# =========================================================================== #
# 1. Genetic code + degenerate codon table
# =========================================================================== #

_BASES = "TCAG"
_AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
AA_BY_CODON = dict(zip([a + b + c for a in _BASES for b in _BASES for c in _BASES],
                       _AA_STRING))
CODONS_BY_AA = defaultdict(list)
for _c, _a in AA_BY_CODON.items():
    if _a != "*":
        CODONS_BY_AA[_a].append(_c)

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def _build_degenerate_table():
    """Every stop-free IUPAC triplet -> (triplet, amino acids produced,
    #degenerate bases), cheapest first (fewest amino acids, then fewest
    degenerate bases)."""
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
                    continue
                ndeg = sum(1 for b in (b1, b2, b3) if len(IUPAC[b]) > 1)
                out.append((b1 + b2 + b3, frozenset(aas), ndeg))
    out.sort(key=lambda r: (len(r[1]), r[2]))
    return out


DEG_TABLE = _build_degenerate_table()
_cache: dict[frozenset, list] = {}


def codons_for(aa_set):
    """All stop-free triplets covering `aa_set`, cheapest first. May be empty --
    that is the case that triggers the fragment fallback."""
    key = frozenset(aa_set)
    hit = _cache.get(key)
    if hit is None:
        hit = [(t, aas, nd) for t, aas, nd in DEG_TABLE if key <= aas]
        _cache[key] = hit
    return hit


# =========================================================================== #
# 2. Forbidden Type IIS sites
# =========================================================================== #

FORBIDDEN_SITES = frozenset({"CGTCTC", "GAGACG"})
BACKBONE_OVERHANGS = frozenset({"CGGA", "GGTG"})


def iupac_may_contain_site(seq):
    """True if SOME expansion of an IUPAC string contains a forbidden site."""
    for site in FORBIDDEN_SITES:
        n = len(site)
        for i in range(len(seq) - n + 1):
            if all(site[k] in IUPAC[seq[i + k]] for k in range(n)):
                return True
    return False


# =========================================================================== #
# 3. Input: base sequence + mutation list
# =========================================================================== #

MUT_RE = re.compile(r"^([A-Z*-])(\d+)([A-Z*-])$")


class Mutation:
    """One parsed mutation.  `pos` is 0-based into the base sequence."""

    __slots__ = ("label", "kind", "pos", "ref", "alt", "status", "reason", "how")

    def __init__(self, label, kind, pos, ref, alt):
        self.label = label
        self.kind = kind          # 'sub' | 'ins' | 'del'
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.status = "pending"   # 'included' | 'dropped'
        self.reason = ""
        self.how = ""             # human-readable encoding used

    def drop(self, reason):
        self.status, self.reason = "dropped", reason


def read_base_sequence(path):
    """FASTA or raw text; returns an uppercase protein string with gaps removed."""
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    s = "".join(seq).replace("-", "").replace("*", "")
    if not s:
        sys.exit(f"no sequence found in {path}")
    return s


def parse_mutations(tokens, base):
    """Parse the DARWIN mutation strings.  Anything malformed, out of range, or
    disagreeing with the base sequence is dropped here with a reason -- it never
    reaches the design."""
    muts, seen = [], {}
    for tok in tokens:
        tok = tok.strip().upper()
        if not tok:
            continue
        m = MUT_RE.match(tok)
        if not m:
            bad = Mutation(tok, "sub", -1, "?", "?")
            bad.drop("unparseable (expected forms A2N, -4R, G5-)")
            muts.append(bad)
            continue
        ref, pos1, alt = m.group(1), int(m.group(2)), m.group(3)
        pos = pos1 - 1

        if ref == "-" and alt == "-":
            bad = Mutation(tok, "sub", pos, ref, alt)
            bad.drop("both sides are gaps")
            muts.append(bad)
            continue

        kind = "ins" if ref == "-" else ("del" if alt == "-" else "sub")
        mu = Mutation(tok, kind, pos, ref, alt)

        # An insertion may sit just past the last residue (append); the others
        # must land on a real residue.
        limit = len(base) if kind == "ins" else len(base) - 1
        if pos < 0 or pos > limit:
            mu.drop(f"position {pos1} outside base sequence (length {len(base)})")
            muts.append(mu)
            continue
        if kind != "ins" and base[pos] != ref:
            mu.drop(f"base sequence has {base[pos]} at position {pos1}, not {ref}")
            muts.append(mu)
            continue
        if kind != "del" and alt not in CODONS_BY_AA:
            mu.drop(f"'{alt}' is not an encodable amino acid")
            muts.append(mu)
            continue

        key = (kind, pos, alt)
        if key in seen:
            mu.drop(f"duplicate of {seen[key]}")
            muts.append(mu)
            continue
        seen[key] = tok
        muts.append(mu)
    return muts


# =========================================================================== #
# 4. Classify: which mutations become codons, which become fragment variants
# =========================================================================== #

def classify(muts, base):
    """Returns (sub_sets, split_sites).

    sub_sets      {pos: frozenset(amino acids)} -- realised as degenerate codons.
                  The wild-type residue is always included, so the library is
                  combinatorial and contains the unmutated parent.
    split_sites   ordered list of sites that force a fragment boundary: every
                  indel, plus any substitution with no stop-free codon."""
    wanted = defaultdict(set)
    for mu in muts:
        if mu.status == "dropped" or mu.kind != "sub":
            continue
        wanted[mu.pos].add(base[mu.pos])
        wanted[mu.pos].add(mu.alt)

    sub_sets, fallback_positions = {}, {}
    for pos, aas in wanted.items():
        if codons_for(aas):
            sub_sets[pos] = frozenset(aas)
        else:
            fallback_positions[pos] = frozenset(aas)

    sites = []
    for mu in muts:
        if mu.status == "dropped":
            continue
        if mu.kind in ("ins", "del"):
            sites.append({"pos": mu.pos, "kind": mu.kind, "muts": [mu],
                          "aas": None, "res": mu.alt if mu.kind == "ins" else None})
    for pos, aas in fallback_positions.items():
        owners = [m for m in muts
                  if m.status != "dropped" and m.kind == "sub" and m.pos == pos]
        sites.append({"pos": pos, "kind": "subsplit", "muts": owners, "aas": aas,
                      "res": None})
    sites.sort(key=lambda s: (s["pos"], s["kind"]))
    return sub_sets, sites


def mutated_positions(sub_sets, sites):
    """Residues a junction may NOT touch: anything carrying a degenerate codon or
    involved in a split site.  A junction pins the codons either side of the cut,
    and a pinned codon cannot also be degenerate or optional."""
    bad = set(sub_sets)
    for s in sites:
        bad.add(s["pos"])
        if s["kind"] == "ins":
            bad.add(s["pos"] - 1)      # the insertion lives between pos-1 and pos
    return bad


# =========================================================================== #
# 5. Junctions
# =========================================================================== #

def overhang_ok(oh, reserved):
    """Individually high-fidelity: non-palindromic (palindromes self-ligate),
    GC-balanced (extreme GC ligates poorly), and not a reserved backbone
    overhang.  Potapov 2018 / Pryor 2020."""
    if oh in reserved or revcomp(oh) in reserved:
        return False
    if oh == revcomp(oh):
        return False
    return 1 <= sum(1 for b in oh if b in "GC") <= 3


def gg_conflict(a, b):
    """Two overhangs cross-react if identical, or if one anneals to the other's
    reverse complement within a single mismatch (the dominant GGA failure)."""
    if a == b:
        return True
    rb = revcomp(b)
    return sum(1 for x, y in zip(a, rb) if x != y) <= 1


def junction_options(base, q, blocked, reserved):
    """Overhangs available for a cut immediately before residue index q.  The
    overhang straddles the cut: last 2 nt of residue q-1's codon + first 2 nt of
    residue q's codon, so both residues must be unmutated (their codons get
    pinned).  Returns [(overhang, left_codon, right_codon), ...]."""
    if q <= 0 or q >= len(base):
        return []
    if (q - 1) in blocked or q in blocked:
        return []
    seen, out = set(), []
    for cl in CODONS_BY_AA[base[q - 1]]:
        for cr in CODONS_BY_AA[base[q]]:
            oh = cl[1:] + cr[:2]
            if oh in seen or not overhang_ok(oh, reserved):
                continue
            seen.add(oh)
            out.append((oh, cl, cr))
    return out


def plan_cuts(base, sites, blocked, reserved, min_frag):
    """Isolate each split site in its own fragment.

    Walk the sites left to right; between two kept sites a cut must be placed at
    a position that offers a legal overhang orthogonal to the ones already
    chosen.  Where no such position exists the later site is DROPPED (this is the
    'impossible to include -- lack of a proper junction' output).  We search for
    the plan that drops the fewest sites."""
    n = len(sites)
    best = {"drops": n + 1, "cuts": None, "tokens": None, "kept": None}

    def dfs(i, prev_pos, cuts, tokens, kept, dropped):
        if len(dropped) >= best["drops"]:
            return                                   # prune: already worse
        if i == n:
            if len(base) - (cuts[-1] if cuts else 0) < min_frag and cuts:
                return
            best.update(drops=len(dropped), cuts=list(cuts), tokens=list(tokens),
                        kept=list(kept))
            return
        site = sites[i]

        # --- option A: keep this site ---------------------------------------
        if prev_pos is None:                          # first kept site: no cut yet
            dfs(i + 1, site["pos"], cuts, tokens, kept + [site], dropped)
        else:
            lo = max(prev_pos + 2, (cuts[-1] if cuts else 0) + min_frag)
            hi = site["pos"] - 1
            for q in range(lo, hi + 1):
                if len(base) - q < min_frag:
                    break
                for oh, cl, cr in junction_options(base, q, blocked, reserved):
                    if any(gg_conflict(oh, t[0]) for t in tokens):
                        continue
                    dfs(i + 1, site["pos"], cuts + [q], tokens + [(oh, cl, cr)],
                        kept + [site], dropped)
                    break            # one legal overhang per position is enough
        # --- option B: drop it ----------------------------------------------
        dfs(i + 1, prev_pos, cuts, tokens, kept, dropped + [site])

    dfs(0, None, [], [], [], [])
    if best["cuts"] is None:
        return [], [], [], list(sites)
    kept_ids = {id(s) for s in best["kept"]}
    dropped = [s for s in sites if id(s) not in kept_ids]
    return best["cuts"], best["tokens"], best["kept"], dropped


# =========================================================================== #
# 6. Build the oligos
# =========================================================================== #

def choose_codons(residues, sub_sets, offset, lead=None, tail=None):
    """One IUPAC codon per residue of a fragment:
      * pinned  -> the codon the junction overhang requires.  The pins are the
                   FIRST and LAST residue of the fragment, addressed by index
                   rather than by sequence position, because an indel inside the
                   fragment shifts every downstream position.
      * mutated -> the cheapest stop-free degenerate codon covering {wt, mutants}
      * else    -> a plain synonymous codon
    Then make sure no expansion of the whole oligo can spell a forbidden Type IIS
    site, repairing by swapping in alternative codons where the design allows it.
    Returns (codons, produced_sets) or None if no site-free assignment exists."""
    codons, produced, flexible = [], [], []
    last = len(residues) - 1
    for i, aa in enumerate(residues):
        pos = offset + i
        if (i == 0 and lead) or (i == last and tail):
            codons.append(lead if i == 0 else tail)
            produced.append(frozenset(aa))
            continue
        if pos in sub_sets:
            opts = codons_for(sub_sets[pos])
            codons.append(opts[0][0])
            produced.append(opts[0][1])
            flexible.append((i, [o[0] for o in opts[:8]],
                             [o[1] for o in opts[:8]]))
            continue
        alts = CODONS_BY_AA[aa]
        codons.append(alts[0])
        produced.append(frozenset(aa))
        flexible.append((i, alts, [frozenset(aa)] * len(alts)))

    if not iupac_may_contain_site("".join(codons)):
        return codons, produced

    for _round in range(len(codons)):
        seq = "".join(codons)
        if not iupac_may_contain_site(seq):
            return codons, produced
        hit = next((i for i in range(len(seq) - 5)
                    if iupac_may_contain_site(seq[i:i + 6])), None)
        if hit is None:
            return codons, produced
        lo_cod, hi_cod = hit // 3, min(len(codons) - 1, (hit + 5) // 3)
        fixed = False
        for idx, alts, prods in flexible:
            if not (lo_cod <= idx <= hi_cod):
                continue
            keep_c, keep_p = codons[idx], produced[idx]
            for cand, prod in zip(alts, prods):
                if cand == keep_c:
                    continue
                codons[idx], produced[idx] = cand, prod
                if not iupac_may_contain_site("".join(codons)):
                    return codons, produced
                if not iupac_may_contain_site("".join(codons)[max(0, hit - 6):hit + 12]):
                    fixed = True
                    break
                codons[idx], produced[idx] = keep_c, keep_p
            if fixed:
                break
        if not fixed:
            return None
    return None


def build_fragments(base, sub_sets, cuts, tokens, kept_sites, reserved):
    """One layer per fragment.  A fragment holding a split site is ordered twice
    (with and without the change); every other fragment is a single oligo.  All
    fragments carry the degenerate codons of the substitutions inside them."""
    bounds = [0] + list(cuts) + [len(base)]
    K = len(bounds) - 1
    # A junction pins the LAST codon of the fragment on its left and the FIRST
    # codon of the fragment on its right.
    lead = [None] * K
    tail = [None] * K
    for f, (oh, cl, cr) in enumerate(tokens):
        tail[f] = cl
        lead[f + 1] = cr

    frags = []
    for f in range(K):
        a, b = bounds[f], bounds[f + 1]
        inside = [s for s in kept_sites if a <= s["pos"] < b]

        # every alternative residue string this fragment must be able to produce
        variants = [("wild-type", list(base[a:b]), a)]
        for s in inside:
            new_variants = []
            for name, res, off in variants:
                seq = list(res)
                idx = s["pos"] - a
                if s["kind"] == "ins":
                    seq.insert(idx, s["res"])
                    lbl = f"+{s['muts'][0].label}"
                elif s["kind"] == "del":
                    del seq[idx]
                    lbl = f"+{s['muts'][0].label}"
                else:                                  # demoted substitution
                    lbl = "+" + "/".join(m.label for m in s["muts"])
                new_variants.append((name + lbl if name != "wild-type" else lbl.lstrip("+"),
                                     seq, off))
            variants += new_variants

        oligos = []
        for name, res, off in variants:
            # A length change shifts every downstream residue of this fragment,
            # so the substitution map is re-expressed in the variant's own
            # coordinates before codons are chosen.
            chosen = choose_codons(res, _shifted_subs(sub_sets, inside, name),
                                   a, lead[f], tail[f])
            if chosen is None:
                continue
            codons, produced = chosen
            nvar = 1
            for p in produced:
                nvar *= len(p)
            oligos.append({"name": name, "oligo": "".join(codons),
                           "codons": codons, "produced": produced,
                           "residues": "".join(res), "variants": nvar,
                           "nt": 3 * len(codons)})
        frags.append({"index": f + 1, "start": a, "end": b, "oligos": oligos,
                      "sites": inside})
    return frags


def _shifted_subs(sub_sets, s_list, variant):
    """Substitution positions inside a fragment, expressed in the coordinates of
    the variant being built.  An insertion/deletion inside the fragment shifts
    every substitution downstream of it by one residue."""
    if not s_list or variant == "wild-type":
        return sub_sets
    applied = [s for s in s_list if s["kind"] in ("ins", "del")
               and (s["muts"][0].label in variant)]
    if not applied:
        return sub_sets
    out = {}
    for pos, aas in sub_sets.items():
        shift = 0
        for s in applied:
            if s["pos"] <= pos:
                shift += 1 if s["kind"] == "ins" else -1
        out[pos + shift] = aas
    return out


# =========================================================================== #
# 7. Library accounting + concrete examples
# =========================================================================== #

def library_size(frags):
    total = 1
    for fr in frags:
        total *= max(1, sum(o["variants"] for o in fr["oligos"]))
    return total


def concrete_codon(iupac_cod, aa):
    for n1 in IUPAC[iupac_cod[0]]:
        for n2 in IUPAC[iupac_cod[1]]:
            for n3 in IUPAC[iupac_cod[2]]:
                if AA_BY_CODON[n1 + n2 + n3] == aa:
                    return n1 + n2 + n3
    return None


def assemble_examples(frags, base, sub_sets, limit=4):
    """Concrete full-length DNA for a few corner combinations (all wild-type, all
    mutant, and a couple of mixes) -- and the proof obligation that none of them
    contains a forbidden Type IIS site."""
    picks = []
    for fr in frags:
        picks.append(list(range(len(fr["oligos"]))))
    combos = list(product(*picks))[:limit] if picks else []
    out, bad = [], 0
    for combo in combos:
        dna, name = [], []
        for fr, oi in zip(frags, combo):
            o = fr["oligos"][oi]
            name.append(o["name"])
            for cod, aa in zip(o["codons"], o["residues"]):
                cc = concrete_codon(cod, aa)
                if cc is None:
                    dna = None
                    break
                dna.append(cc)
            if dna is None:
                break
        if dna is None:
            continue
        seq = "".join(dna)
        if any(s in seq for s in FORBIDDEN_SITES):
            bad += 1
        out.append(("|".join(name), seq))
    return out, bad


# =========================================================================== #
# 8. Report
# =========================================================================== #

def build_report(args, base, muts, sub_sets, frags, cuts, tokens, examples, bad):
    lib = library_size(frags)
    oligos = sum(len(fr["oligos"]) for fr in frags)
    nt = sum(o["nt"] for fr in frags for o in fr["oligos"])
    L = []
    L.append(f"base sequence: {len(base)} aa" +
             (f"  ({args.seq_file})" if args.seq_file else ""))
    L.append(f"mutations given: {len(muts)}   chemistry: Golden Gate "
             f"({args.gg_enzyme.upper()}, site {sorted(FORBIDDEN_SITES)[0]})")
    L.append("")

    L.append("MUTATIONS")
    L.append(f"  {'label':<10} {'type':<12} {'status':<9} how / why")
    for mu in muts:
        L.append(f"  {mu.label:<10} {mu.kind:<12} {mu.status:<9} "
                 f"{mu.how if mu.status == 'included' else mu.reason}")
    n_in = sum(1 for m in muts if m.status == "included")
    L.append(f"  -> {n_in}/{len(muts)} realised")
    L.append("")

    L.append(f"DESIGN: {len(frags)} fragment(s)")
    for fr in frags:
        degpos = sorted({fr['start'] + i + 1
                         for o in fr["oligos"]
                         for i, p in enumerate(o["produced"]) if len(p) > 1})
        L.append(f"  fragment {fr['index']}: residues {fr['start']+1}-{fr['end']} "
                 f"({fr['end']-fr['start']} aa), {len(fr['oligos'])} oligo(s), "
                 f"{sum(o['nt'] for o in fr['oligos'])} nt")
        for o in fr["oligos"]:
            L.append(f"      - {o['name']:<28} {o['variants']:>4} variant(s), {o['nt']} nt")
        if degpos:
            L.append(f"      degenerate codons at residue(s): "
                     f"{', '.join(str(p) for p in degpos)}")
    L.append("")

    if cuts:
        L.append("JUNCTIONS (Golden Gate) -- individually high-fidelity, mutually orthogonal:")
        for q, (oh, cl, cr) in zip(cuts, tokens):
            L.append(f"  before residue {q+1}: overhang 5'-{oh}-3' (rc {revcomp(oh)}), "
                     f"pinned codons {cl}|{cr}  [{base[q-1]}{q}|{base[q]}{q+1}]")
        L.append(f"  [backbone overhangs {' | '.join(sorted(BACKBONE_OVERHANGS))} "
                 f"reserved -- excluded from internal junctions]")
        L.append("")

    L.append("ORDER")
    L.append(f"  oligos to order : {oligos}")
    L.append(f"  total synthesis : {nt:,} nt")
    L.append(f"  proteins produced (all combinations): {lib:,}")
    L.append("")
    L.append(f"forbidden-site check on {len(examples)} assembled full-length "
             f"example(s): {'FAIL' if bad else 'clean'}")
    return "\n".join(L)


def save_run(out_root, args, report, base, muts, sub_sets, frags, cuts, tokens,
             examples):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = (os.path.splitext(os.path.basename(args.seq_file))[0]
            if args.seq_file else "seq")
    run = os.path.join(out_root, f"{ts}_{stem}_darwin")
    os.makedirs(run, exist_ok=True)

    with open(os.path.join(run, "report.txt"), "w") as fh:
        fh.write(report + "\n")

    summary = {
        "base_length": len(base),
        "args": {k: v for k, v in vars(args).items()},
        "mutations": [{"label": m.label, "kind": m.kind, "position_1based": m.pos + 1,
                       "status": m.status, "how": m.how, "reason": m.reason}
                      for m in muts],
        "degenerate_positions": {str(p + 1): "".join(sorted(a))
                                 for p, a in sorted(sub_sets.items())},
        "junctions": [{"before_residue_1based": q + 1, "overhang": oh,
                       "left_codon": cl, "right_codon": cr}
                      for q, (oh, cl, cr) in zip(cuts, tokens)],
        "fragments": [{"index": fr["index"],
                       "residues_1based": [fr["start"] + 1, fr["end"]],
                       "oligos": [{"name": o["name"], "oligo": o["oligo"],
                                   "variants": o["variants"], "nt": o["nt"]}
                                  for o in fr["oligos"]]}
                      for fr in frags],
        "library_size": library_size(frags),
        "total_oligos": sum(len(fr["oligos"]) for fr in frags),
        "total_nt": sum(o["nt"] for fr in frags for o in fr["oligos"]),
    }
    with open(os.path.join(run, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    with open(os.path.join(run, "oligos.fasta"), "w") as fh:
        for fr in frags:
            for i, o in enumerate(fr["oligos"], start=1):
                fh.write(f">frag{fr['index']}_oligo{i}_{o['name'].replace(' ', '')}"
                         f"_var{o['variants']}\n{o['oligo']}\n")
    if examples:
        with open(os.path.join(run, "examples_full_length_dna.fasta"), "w") as fh:
            for name, seq in examples:
                fh.write(f">{name}\n{seq}\n")
    return run


# =========================================================================== #
# 9. CLI
# =========================================================================== #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--seq-file", help="FASTA (or raw text) holding the base protein")
    src.add_argument("--seq", help="base protein sequence given inline")
    mut = ap.add_mutually_exclusive_group(required=True)
    mut.add_argument("--mutations",
                     help="comma- or space-separated list, e.g. \"A2N,-4R,G5-,W40K\". "
                          "Use the --mutations=\"...\" form: an insertion token "
                          "starts with '-', which a bare argument list would read "
                          "as a flag.")
    mut.add_argument("--mutations-file", help="one mutation per line (# comments ok)")
    ap.add_argument("--min-frag-aa", type=int, default=15,
                    help="minimum fragment length in residues (default 15)")
    ap.add_argument("--gg-enzyme", choices=["bsmbi", "esp3i", "bsai"], default="bsmbi",
                    help="Type IIS enzyme whose site is banned everywhere")
    ap.add_argument("--shared-backbone-overhangs", action="store_true",
                    help="allow the reserved backbone overhangs at internal junctions")
    ap.add_argument("--out-dir", default="darwinruns")
    args = ap.parse_args()

    global FORBIDDEN_SITES
    site = {"bsmbi": "CGTCTC", "esp3i": "CGTCTC", "bsai": "GGTCTC"}[args.gg_enzyme]
    FORBIDDEN_SITES = frozenset({site, revcomp(site)})
    reserved = frozenset() if args.shared_backbone_overhangs else BACKBONE_OVERHANGS

    base = (read_base_sequence(args.seq_file) if args.seq_file
            else args.seq.strip().upper().replace("-", ""))
    if args.mutations_file:
        with open(args.mutations_file) as fh:
            tokens = [ln.split("#")[0].strip() for ln in fh]
    else:
        tokens = re.split(r"[,\s]+", args.mutations.strip())

    muts = parse_mutations(tokens, base)
    sub_sets, sites = classify(muts, base)
    blocked = mutated_positions(sub_sets, sites)
    cuts, tokens_gg, kept, dropped = plan_cuts(base, sites, blocked, reserved,
                                               args.min_frag_aa)

    for s in dropped:
        for mu in s["muts"]:
            mu.drop("no orthogonal Golden Gate junction available to isolate it")
    # a demoted substitution that got dropped must also lose its degenerate codon
    frags = build_fragments(base, sub_sets, cuts, tokens_gg, kept, reserved)

    for mu in muts:
        if mu.status != "pending":
            continue
        if mu.kind == "sub" and mu.pos in sub_sets:
            aas = "".join(sorted(sub_sets[mu.pos]))
            cod = codons_for(sub_sets[mu.pos])[0][0]
            mu.status, mu.how = "included", (f"degenerate codon {cod} at residue "
                                             f"{mu.pos+1} (encodes {aas})")
        else:
            site_of = next((s for s in kept if mu in s["muts"]), None)
            if site_of is None:
                mu.drop("no fragment carries this change")
            else:
                mu.status = "included"
                mu.how = ("fragment variant (2 oligos: with and without) "
                          f"at residue {mu.pos+1}")

    examples, bad = assemble_examples(frags, base, sub_sets)
    report = build_report(args, base, muts, sub_sets, frags, cuts, tokens_gg,
                          examples, bad)
    run = save_run(args.out_dir, args, report, base, muts, sub_sets, frags, cuts,
                   tokens_gg, examples)
    print(report)
    print(f"\nsaved run to {run}/")


if __name__ == "__main__":
    main()
