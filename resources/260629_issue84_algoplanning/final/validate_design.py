#!/usr/bin/env python3
"""
Independent validity harness for a cutsearch_design.py run.  Issue #84.

DELIBERATELY IMPORTS NOTHING FROM cutsearch_design.py.  It reads the run
directory (summary.json + fragment*.fasta) and the input alignment, and
re-derives the codon table, the IUPAC map, the forbidden sites, the overhangs,
the translations and the coverage from scratch.  If the checker shared code with
the designer, a bug would be mirrored rather than caught -- so every constant
below is a second implementation on purpose.

Exits nonzero if any check FAILs, so it can gate an order.

  python validate_design.py runs/20260826_..._gg_unified
  python validate_design.py runs/... --check-reproducibility
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
from collections import defaultdict

# --------------------------------------------------------------------------- #
# Second implementations of everything the designer also defines.
# --------------------------------------------------------------------------- #

_B = "TTTTTCTTATTGCTTCTCCTACTGATTATCATAATGGTTGTCGTAGTG"  # unused, guards typos
BASES = "TCAG"
AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: aa for (a, b, c), aa in
               zip([(x, y, z) for x in BASES for y in BASES for z in BASES],
                   AA_STRING)}
STOPS = {c for c, aa in CODON_TABLE.items() if aa == "*"}

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

ENZYME_SITE = {"bsmbi": "CGTCTC", "esp3i": "CGTCTC", "bsai": "GGTCTC"}
BACKBONE_OVERHANGS = ("CGGA", "GGTG")
# Low-usage E. coli codons (Kane 1995 / the set Rosetta and RIL supplement).
RARE_ECOLI = {"AGA": "Arg", "AGG": "Arg", "CGA": "Arg", "CGG": "Arg",
              "ATA": "Ile", "CTA": "Leu", "CCC": "Pro", "GGA": "Gly",
              "GGG": "Gly", "TCA": "Ser"}
POOL_LIMIT_NT = 300          # informational only


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def iupac_expand_codon(trip):
    """Every concrete triplet an IUPAC triplet can spell."""
    return ["".join(t) for t in itertools.product(*(IUPAC[b] for b in trip))]


def aa_set_of(trip):
    """Amino acids an IUPAC triplet produces (may include '*')."""
    return {CODON_TABLE[c] for c in iupac_expand_codon(trip)}


def may_contain(seq, sites):
    """True if SOME expansion of an IUPAC string contains one of `sites`."""
    for site in sites:
        n = len(site)
        for i in range(len(seq) - n + 1):
            if all(site[k] in IUPAC.get(seq[i + k], "") for k in range(n)):
                return True
    return False


def codons(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #

def read_alignment(path):
    """[(header, aligned_seq, weight)] -- weight from the `_n<k>` suffix."""
    out, hdr, buf = [], None, []

    def flush():
        if hdr is not None and buf:
            w = 1
            if "_n" in hdr:
                t = hdr.rsplit("_n", 1)[1]
                if t.isdigit():
                    w = int(t)
            out.append((hdr, "".join(buf), w))

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                hdr, buf = line[1:], []
            else:
                buf.append(line.upper())
    flush()
    return out


class Result:
    def __init__(self):
        self.rows = []

    def add(self, name, status, detail=""):
        self.rows.append((name, status, detail))

    @property
    def failed(self):
        return any(s == "FAIL" for _, s, _ in self.rows)


# --------------------------------------------------------------------------- #
# CHECK 1 -- Golden Gate overhangs
# --------------------------------------------------------------------------- #

def check_overhangs(R, S, frags, cuts, notes):
    """Level 1 per site and Level 2 over the SET, both recomputed from the
    emitted DNA rather than from the saved token list -- checking the tokens
    would only re-read what the designer wrote down."""
    if S["args"]["chemistry"] != "gg":
        R.add("1. GG overhangs", "SKIP", "chemistry is not gg")
        return
    if not cuts:
        R.add("1. GG overhangs", "SKIP", "K=1, no internal junctions")
        return

    reserved = () if S["args"].get("shared_backbone_overhangs") else BACKBONE_OVERHANGS

    # -- uniformity: the overhang only exists if EVERY variant spells it ----- #
    bad_unif = []
    for f in range(len(cuts)):
        left_tails = {codons(u["oligo"])[-1] for u in frags[f]}
        right_heads = {codons(u["oligo"])[0] for u in frags[f + 1]}
        if len(left_tails) != 1:
            bad_unif.append(f"junction {f+1}: left oligos end with "
                            f"{len(left_tails)} different codons {sorted(left_tails)}")
        if len(right_heads) != 1:
            bad_unif.append(f"junction {f+1}: right oligos start with "
                            f"{len(right_heads)} different codons {sorted(right_heads)}")
    R.add("1a. overhang uniformity across variants",
          "FAIL" if bad_unif else "PASS",
          "; ".join(bad_unif) if bad_unif else
          f"all oligos agree at every one of {len(cuts)} junction(s)")

    # -- build the overhangs from the DNA ------------------------------------ #
    ohs = []
    for f in range(len(cuts)):
        cl = codons(frags[f][0]["oligo"])[-1]
        cr = codons(frags[f + 1][0]["oligo"])[0]
        ohs.append((cuts[f], cl[1:] + cr[:2], cl, cr))

    # -- Level 1 -------------------------------------------------------------- #
    l1 = []
    for p, oh, _cl, _cr in ohs:
        if len(oh) != 4 or any(b not in "ACGT" for b in oh):
            l1.append(f"col {p}: overhang {oh!r} is not 4 concrete bases")
            continue
        if oh == revcomp(oh):
            l1.append(f"col {p}: {oh} is self-complementary (palindrome)")
        gc = sum(1 for b in oh if b in "GC")
        if not 1 <= gc <= 3:
            l1.append(f"col {p}: {oh} GC={gc}/4 outside 25-75%")
        for r in reserved:
            if oh == r or revcomp(oh) == r:
                l1.append(f"col {p}: {oh} collides with reserved backbone {r}")
    R.add("1b. overhangs individually valid (Level 1)",
          "FAIL" if l1 else "PASS",
          "; ".join(l1) if l1 else
          "  ".join(f"col {p}:{oh}" for p, oh, _, _ in ohs))

    # -- Level 2, INCLUDING the backbone overhangs --------------------------- #
    # The designer excludes the backbone overhangs from internal junctions but
    # never tests mutual orthogonality against them: overhang_ok() rejects only
    # EXACT matches, and gg_conflict() is applied pairwise among chosen tokens
    # only.  An internal overhang one mismatch from CGGA would pass every check
    # in the designer and still mis-ligate into the vector.
    def conflict(a, b):
        if a == b:
            return True
        rb = revcomp(b)
        return sum(1 for x, y in zip(a, rb) if x != y) <= 1

    full = [(f"col {p}", oh) for p, oh, _, _ in ohs]
    full += [(f"backbone", r) for r in reserved]
    l2 = []
    for (na, a), (nb, b) in itertools.combinations(full, 2):
        if conflict(a, b):
            l2.append(f"{na} {a} conflicts with {nb} {b} "
                      f"(identical or <=1 mismatch to its revcomp)")
    R.add("1c. overhang set mutually orthogonal (Level 2, incl. backbone)",
          "FAIL" if l2 else "PASS",
          "; ".join(l2) if l2 else
          f"{len(full)} overhangs pairwise orthogonal "
          f"({len(ohs)} internal + {len(reserved)} backbone)")

    # -- geometric non-overlap ------------------------------------------------ #
    gaps = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
    bad = [g for g in gaps if g * 3 < 4]
    R.add("1d. junctions do not physically overlap",
          "FAIL" if bad else "PASS",
          f"min gap between cuts = {min(gaps)} cols "
          f"({min(gaps)*3} nt) vs 4 nt overhang" if gaps else "single junction")

    notes["overhangs"] = ohs


# --------------------------------------------------------------------------- #
# CHECK 5 -- reading frame.  RUNS FIRST and gates everything else: a truncated
# oligo produces a 2-base "codon", which makes every translation-based check
# below raise rather than report.  Frame is the precondition for the rest.
# --------------------------------------------------------------------------- #

def check_frame(R, frags):
    bad = [f"frag{f+1}#{i+1} len {len(u['oligo'])}"
           for f, units in enumerate(frags)
           for i, u in enumerate(units) if len(u["oligo"]) % 3]
    R.add("5. every oligo is a whole number of codons",
          "FAIL" if bad else "PASS",
          "; ".join(bad) if bad else
          f"all {sum(len(u) for u in frags)} oligos divisible by 3")
    return not bad


# --------------------------------------------------------------------------- #
# CHECK 2 -- forbidden Type IIS sites
# --------------------------------------------------------------------------- #

def check_forbidden(R, S, frags):
    enz = S["args"].get("gg_enzyme", "bsmbi")
    site = ENZYME_SITE.get(enz, "CGTCTC")
    sites = (site, revcomp(site))

    bad = [f"fragment {f+1} oligo {i+1}"
           for f, units in enumerate(frags)
           for i, u in enumerate(units) if may_contain(u["oligo"], sites)]
    R.add("2a. no forbidden site within any oligo",
          "FAIL" if bad else "PASS",
          "; ".join(bad) if bad else
          f"{sum(len(u) for u in frags)} oligos clear of {site}/{revcomp(site)}")

    # A 6-nt site spans at most one junction, so checking every ADJACENT PAIR
    # covers every assembled full-length sequence -- and it covers all of them,
    # not the 3 examples the designer spot-checks.
    if len(frags) > 1:
        bad2, npairs = [], 0
        for f in range(len(frags) - 1):
            for i, a in enumerate(frags[f]):
                for j, b in enumerate(frags[f + 1]):
                    npairs += 1
                    tail, head = a["oligo"][-6:], b["oligo"][:6]
                    if may_contain(tail + head, sites):
                        bad2.append(f"frag{f+1}#{i+1} | frag{f+2}#{j+1}")
        R.add("2b. no forbidden site ACROSS any fragment junction",
              "FAIL" if bad2 else "PASS",
              "; ".join(bad2[:6]) if bad2 else
              f"{npairs} adjacent oligo pairs checked (all junctions, "
              f"not a sample)")
    else:
        R.add("2b. no forbidden site ACROSS any fragment junction", "SKIP",
              "K=1, no junctions")


# --------------------------------------------------------------------------- #
# CHECK 3 -- the library really does encode the claimed natural sequences
# --------------------------------------------------------------------------- #

def check_coverage(R, S, frags, cuts, aln, notes):
    L = len(aln[0][1])
    bounds = [0] + list(cuts) + [L]
    K = len(bounds) - 1
    have_cols = all("cols" in u for units in frags for u in units)

    # translate each oligo back to per-column amino-acid sets, from the DNA
    exps = []
    for units in frags:
        layer = []
        for u in units:
            layer.append((tuple(u.get("cols", ())),
                          [aa_set_of(c) for c in codons(u["oligo"])]))
        exps.append(layer)

    covered, weight = [], 0
    for idx, (hdr, seq, w) in enumerate(aln):
        ok = True
        for f in range(K):
            cols = tuple(j for j in range(bounds[f], bounds[f + 1]) if seq[j] != "-")
            res = [seq[j] for j in cols]
            hit = False
            for ucols, uexp in exps[f]:
                if have_cols:
                    if ucols != cols:
                        continue
                elif len(uexp) != len(res):
                    continue
                if all(r in e for r, e in zip(res, uexp)):
                    hit = True
                    break
            if not hit:
                ok = False
                break
        if ok:
            covered.append(idx)
            weight += w

    rep = [r for r in S["frontier"] if r["K"] == S["recommended_K"]][0]
    exact = "exact (gap patterns from saved cols)" if have_cols else \
            "UPPER BOUND (cols missing; matched on column count)"
    good = (len(covered) == rep["n_cores_encoded"] and weight == rep["encoded_weight"])
    R.add("3. claimed cores are actually producible",
          "PASS" if good else "FAIL",
          f"recomputed {len(covered)} cores / weight {weight}; "
          f"reported {rep['n_cores_encoded']} / {rep['encoded_weight']}  [{exact}]")

    # The designer's cov[] sets only ever grow, so a widen that dropped an
    # over-included amino acid would leave coverage overstated and never
    # corrected.  The recount above is exactly the test for that.
    R.add("3b. coverage monotonicity assumption held",
          "PASS" if good else "FAIL",
          "independent recount matches, so no widen silently un-covered a core"
          if good else "recount disagrees -- see 3")
    notes["covered"] = covered


# --------------------------------------------------------------------------- #
# CHECKS 4, 5, 9 -- codon-level sanity
# --------------------------------------------------------------------------- #

def check_codons(R, S, frags, cuts, aln, notes):
    # -- 4. stop codons ------------------------------------------------------- #
    bad = []
    for f, units in enumerate(frags):
        for i, u in enumerate(units):
            for ci, trip in enumerate(codons(u["oligo"])):
                if "*" in aa_set_of(trip):
                    hits = sorted(set(iupac_expand_codon(trip)) & STOPS)
                    bad.append(f"frag{f+1}#{i+1} codon {ci} ({trip} -> {hits})")
    R.add("4. no expansion of any codon is a stop",
          "FAIL" if bad else "PASS",
          "; ".join(bad[:6]) if bad else
          f"{sum(len(codons(u['oligo'])) for units in frags for u in units)} "
          f"codons, none can spell TAA/TAG/TGA")

    # -- 9. pinned junction codons still encode the right residue ------------- #
    if S["args"]["chemistry"] != "gg" or not cuts:
        R.add("9. junction residues preserved", "SKIP", "no gg junctions")
        return
    L = len(aln[0][1])
    bad = []
    for f, p in enumerate(cuts):
        left_aa = {s[p - 1] for _h, s, _w in aln if s[p - 1] != "-"}
        right_aa = {s[p] for _h, s, _w in aln if s[p] != "-"}
        cl = codons(frags[f][0]["oligo"])[-1]
        cr = codons(frags[f + 1][0]["oligo"])[0]
        for label, cod, want in (("left", cl, left_aa), ("right", cr, right_aa)):
            got = aa_set_of(cod)
            if not want <= got:
                bad.append(f"col {p} {label}: pinned {cod} makes {sorted(got)}, "
                           f"cores need {sorted(want)}")
    R.add("9. pinned junction codons encode the residues the cores need",
          "FAIL" if bad else "PASS",
          "; ".join(bad) if bad else
          f"{2*len(cuts)} pinned codons all consistent with the alignment")


# --------------------------------------------------------------------------- #
# CHECK 6 -- true distinct-protein count vs the reported product
# --------------------------------------------------------------------------- #

def check_library_size(R, S, frags, cap=2_000_000):
    rep = [r for r in S["frontier"] if r["K"] == S["recommended_K"]][0]
    claimed = rep["library"]
    per_layer = []
    for units in frags:
        seqs = set()
        for u in units:
            sets = [sorted(aa_set_of(c)) for c in codons(u["oligo"])]
            n = 1
            for s in sets:
                n *= len(s)
            if n > cap:
                R.add("6. reported library size vs true distinct proteins",
                      "WARN", f"unit expands to {n:,} variants; skipped "
                              f"enumeration above {cap:,}")
                return
            for combo in itertools.product(*sets):
                seqs.add("".join(combo))
        per_layer.append(sorted(seqs))
    total = 1
    for s in per_layer:
        total *= len(s)
    if total > cap:
        R.add("6. reported library size vs true distinct proteins", "WARN",
              f"product is {total:,}; skipped enumeration above {cap:,}")
        return
    whole = {"".join(c) for c in itertools.product(*per_layer)}
    status = "PASS" if len(whole) <= claimed else "FAIL"
    R.add("6. reported library size vs true distinct proteins", status,
          f"reported |L_O| = {claimed:,} (product of per-layer counts, an upper "
          f"bound); true distinct proteins = {len(whole):,}"
          + ("" if len(whole) == claimed else
             f" -- {claimed - len(whole):,} collisions, so the cap is "
             f"conservative in the safe direction"))


# --------------------------------------------------------------------------- #
# CHECK 7 -- the caps were actually respected
# --------------------------------------------------------------------------- #

def check_caps(R, S, frags):
    a = S["args"]
    rep = [r for r in S["frontier"] if r["K"] == S["recommended_K"]][0]
    lines, bad = [], False

    if a.get("max_library"):
        ok = rep["library"] <= a["max_library"]
        bad |= not ok
        lines.append(f"library {rep['library']:,} <= {a['max_library']:,} "
                     f"{'OK' if ok else 'VIOLATED'}")
    if a.get("max_nt"):
        got = rep.get("nt_ordered") or rep["nt"]
        ok = got <= a["max_nt"]
        bad |= not ok
        lines.append(f"nt ordered {got:,} <= {a['max_nt']:,} "
                     f"{'OK' if ok else 'VIOLATED'}")
    if a.get("max_oligo_nt"):
        longest = max(len(u["oligo"]) for units in frags for u in units)
        ok = longest <= a["max_oligo_nt"]
        bad |= not ok
        lines.append(f"longest oligo {longest} <= {a['max_oligo_nt']} "
                     f"{'OK' if ok else 'VIOLATED'}")
    if a.get("max_junk_pct", 100) < 100:
        ok = rep["junk_pct"] <= a["max_junk_pct"] + 1e-9
        bad |= not ok
        lines.append(f"junk {rep['junk_pct']:.1f}% <= {a['max_junk_pct']}% "
                     f"{'OK' if ok else 'VIOLATED'}")
    if not lines:
        R.add("7. caps respected", "SKIP", "no caps were set")
        return
    R.add("7. caps respected", "FAIL" if bad else "PASS", "; ".join(lines))


# --------------------------------------------------------------------------- #
# CHECK 11 -- will a vendor actually make this, and will E. coli translate it?
#
# Nothing upstream guarantees either.  best_codon() ranks by amino-acid count
# then degeneracy; before the codon-usage tie-break was added it fell through to
# alphabetical order and chose GGA/ATA/CTA/AGA -- the four rarest E. coli codons
# -- for 24-29% of positions.  Maximising usage instead pushed 50-nt windows to
# 84% GC.  Both are caught here.
# --------------------------------------------------------------------------- #

def _gc(s):
    return 100.0 * sum(1 for b in s if b in "GC") / len(s) if s else 0.0


def _extreme_expansions(oligo):
    """An oligo with degenerate positions is a family; take the GC-richest and
    GC-poorest concrete members, which bracket every other one."""
    hi = "".join(max(IUPAC[b], key=lambda x: x in "GC") for b in oligo)
    lo = "".join(min(IUPAC[b], key=lambda x: x in "GC") for b in oligo)
    return hi, lo


def check_synthesis(R, frags, window=50, gc_lo=25.0, gc_hi=75.0, max_run=10):
    olig = [(f + 1, i + 1, u["oligo"])
            for f, us in enumerate(frags) for i, u in enumerate(us)]

    seen, dups = {}, []
    for f, i, o in olig:
        if o in seen:
            dups.append(f"frag{f}#{i} identical to frag{seen[o][0]}#{seen[o][1]}")
        seen[o] = (f, i)
    R.add("11a. no duplicate oligos in the order", "WARN" if dups else "PASS",
          "; ".join(dups[:4]) if dups else f"{len(olig)} oligos all distinct")

    bad_win, lo_all, hi_all, worst_run = [], 100.0, 0.0, 0
    for f, i, o in olig:
        for conc in _extreme_expansions(o):
            lo_all, hi_all = min(lo_all, _gc(conc)), max(hi_all, _gc(conc))
            run, prev = 1, ""
            for b in conc:
                run = run + 1 if b == prev else 1
                prev = b
                worst_run = max(worst_run, run)
            if len(conc) >= window:
                w = [_gc(conc[k:k + window])
                     for k in range(len(conc) - window + 1)]
                if max(w) > gc_hi or min(w) < gc_lo:
                    bad_win.append(f"frag{f}#{i} ({min(w):.0f}-{max(w):.0f}%)")
    R.add(f"11b. GC within {gc_lo:.0f}-{gc_hi:.0f}% in every {window}nt window",
          "WARN" if bad_win else "PASS",
          f"{len(set(bad_win))} oligos out of range: "
          + "; ".join(sorted(set(bad_win))[:4]) if bad_win else
          f"overall GC {lo_all:.1f}-{hi_all:.1f}% across all expansions")
    R.add(f"11c. no homopolymer run > {max_run}",
          "WARN" if worst_run > max_run else "PASS",
          f"longest run {worst_run} (worst-case expansion)")

    tot, forced = 0, defaultdict(int)
    for _f, _i, o in olig:
        for trip in codons(o):
            tot += 1
            exp = iupac_expand_codon(trip)
            if exp and all(c in RARE_ECOLI for c in exp):
                forced[trip] += 1
    n = sum(forced.values())
    pct = 100.0 * n / tot if tot else 0.0
    R.add("11d. rare E. coli codons (expression)",
          "WARN" if pct > 10.0 else "PASS",
          f"{n}/{tot} positions ({pct:.1f}%) forced to a rare codon"
          + (f": {dict(list(forced.items())[:4])}" if forced else "")
          + "   [native genes run 5-10%]")


# --------------------------------------------------------------------------- #
# CHECK 10 -- the recorded seed reproduces the design
# --------------------------------------------------------------------------- #

def check_reproducible(R, S, run_dir, python_exe):
    import glob
    import shutil
    import tempfile

    a = dict(S["args"])
    aln = a.get("aln_fasta")
    if not aln or not os.path.exists(aln):
        aln = S.get("input")
    if not aln or not os.path.exists(aln):
        R.add("10. seed reproduces the design", "SKIP",
              "input alignment not found at the recorded path")
        return

    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "cutsearch_design.py")
    tmp = tempfile.mkdtemp(prefix="revalidate_")
    cmd = [python_exe, script, aln, "--out-dir", tmp]
    flags = {
        "chemistry": "--chemistry", "max_junk_pct": "--max-junk-pct",
        "k_max": "--k-max", "min_block_cols": "--min-block-cols",
        "arm_codons": "--arm-codons", "widen_candidates": "--widen-candidates",
        "gg_enzyme": "--gg-enzyme", "cut_candidates": "--cut-candidates",
        "proxy_candidates": "--proxy-candidates", "seed": "--seed",
        "exhaustive_max": "--exhaustive-max", "max_oligo_nt": "--max-oligo-nt",
        "cut_node_budget": "--cut-node-budget", "max_library": "--max-library",
        "max_nt": "--max-nt", "oligo_overhead_nt": "--oligo-overhead-nt",
    }
    for k, flag in flags.items():
        if a.get(k) is not None:
            cmd += [flag, str(a[k])]
    if a.get("shared_backbone_overhangs"):
        cmd.append("--shared-backbone-overhangs")
    if a.get("rank_seqs_per_kb"):
        cmd.append("--rank-seqs-per-kb")

    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if p.returncode != 0:
            R.add("10. seed reproduces the design", "FAIL",
                  f"re-run exited {p.returncode}: {p.stderr.strip()[-300:]}")
            return
        newdirs = sorted(glob.glob(os.path.join(tmp, "*")))
        if not newdirs:
            R.add("10. seed reproduces the design", "FAIL", "re-run wrote nothing")
            return
        diffs = []
        for f in sorted(glob.glob(os.path.join(run_dir, "fragment*.fasta"))):
            other = os.path.join(newdirs[-1], os.path.basename(f))
            if not os.path.exists(other):
                diffs.append(f"{os.path.basename(f)} missing from re-run")
            elif open(f).read() != open(other).read():
                diffs.append(f"{os.path.basename(f)} differs")
        R.add("10. seed reproduces the design", "FAIL" if diffs else "PASS",
              "; ".join(diffs) if diffs else
              f"re-run with seed {a.get('seed')} produced byte-identical "
              f"fragment FASTAs")
    except subprocess.TimeoutExpired:
        R.add("10. seed reproduces the design", "WARN", "re-run timed out")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a run directory written by cutsearch_design.py")
    ap.add_argument("--aln", default=None,
                    help="input alignment (default: the path recorded in summary.json)")
    ap.add_argument("--check-reproducibility", action="store_true",
                    help="re-run the designer with the recorded seed and diff the "
                         "oligos (doubles runtime)")
    ap.add_argument("--python", default=sys.executable,
                    help="interpreter to use for the reproducibility re-run")
    args = ap.parse_args()

    with open(os.path.join(args.run_dir, "summary.json")) as fh:
        S = json.load(fh)
    frags = S["fragments"]
    cuts = S.get("cuts", [])
    # Take the first candidate that actually EXISTS.  `aln_fasta` is recorded as
    # it was typed, so it is usually relative to the directory the design was run
    # from -- which is not where the validator runs, especially once a run
    # directory has been copied elsewhere (negative_controls.py does exactly
    # that).  `input` is the absolute path saved alongside it.
    candidates = [args.aln, S["args"].get("aln_fasta"), S.get("input")]
    aln_path = next((c for c in candidates if c and os.path.exists(c)), None)
    if aln_path is None:
        tried = "\n  ".join(str(c) for c in candidates if c)
        sys.exit(f"input alignment not found; tried:\n  {tried}\n"
                 f"pass it explicitly with --aln")
    aln = read_alignment(aln_path)

    R, notes = Result(), {}
    print(f"run    : {args.run_dir}")
    print(f"input  : {aln_path}")
    print(f"design : K={S['recommended_K']}  cuts={cuts}  "
          f"chemistry={S['args']['chemistry']}  seed={S.get('seed')}")
    print(f"cores  : {len(aln)} unique / {sum(w for _h,_s,w in aln)} natural\n")

    if check_frame(R, frags):
        check_overhangs(R, S, frags, cuts, notes)
        check_forbidden(R, S, frags)
        check_coverage(R, S, frags, cuts, aln, notes)
        check_codons(R, S, frags, cuts, aln, notes)
        check_library_size(R, S, frags)
        check_synthesis(R, frags)
    else:
        for nm in ("1. GG overhangs", "2. forbidden sites", "3. coverage",
                   "4. stop codons", "6. true library size",
                   "9. junction residues", "11. synthesis complexity"):
            R.add(nm, "SKIP", "reading frame is broken; cannot translate")
    check_caps(R, S, frags)
    if args.check_reproducibility:
        check_reproducible(R, S, args.run_dir, args.python)

    w = max(len(n) for n, _, _ in R.rows)
    for name, status, detail in R.rows:
        print(f"  [{status:^4}] {name:<{w}}  {detail}")
    n_fail = sum(1 for _, s, _ in R.rows if s == "FAIL")
    n_warn = sum(1 for _, s, _ in R.rows if s == "WARN")
    print(f"\n{len(R.rows)} checks: {len(R.rows)-n_fail-n_warn} pass, "
          f"{n_fail} fail, {n_warn} warn")
    if R.failed:
        print("VALIDATION FAILED -- do not order this design.")
    sys.exit(1 if R.failed else 0)


if __name__ == "__main__":
    main()
