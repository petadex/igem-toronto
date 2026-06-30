#!/usr/bin/env python3
"""
Fragment-split oligo design -- ROUND-1 PROTOTYPE (3 fragments, 2 cuts).

This is the front end that `greedy_oligo.py` was missing. Instead of ONE
degenerate oligo across the whole core (which forces equal length and can only
encode a rectangular product of per-position choices), we cut the aligned cores
into 3 contiguous fragments at 2 chosen boundaries, synthesize each fragment's
variants separately, and assemble them combinatorially (Golden Gate / HR). A
producible full-length gene = one choice of fragment from each of the 3 layers.

WHY 2 cuts and not a per-fragment greedy: coverage is conjunctive over whole
cores (a core is realized only if ALL 3 of its fragment pieces are present and
co-assemble), so you cannot optimize fragments independently and glue. The cut
placement is the real decision variable; this script makes it explicitly.

----------------------------------------------------------------------------- #
RELATION TO GGAssembler (Hoch & Fleishman, Protein Science 2024)
----------------------------------------------------------------------------- #
GGAssembler models a combinatorial library as a graph where NODES are
restriction-enzyme cleavage points (candidate overhang positions) and EDGES are
the DNA fragment between two cleavage points; it then runs a (rainbow) shortest
path to find the cheapest segmentation, where an edge's cost is the number of
nucleotides needed to encode all the diversity in that fragment, and "rainbow"
enforces that the overhangs chosen along the path are mutually orthogonal /
high-fidelity (a hard, ~limited resource). Degenerate codons are exploited
because widening an amino-acid set within one codon is free.

This prototype is the SAME graph, specialized and turned inside out:

  GGAssembler concept            <->   here
  ---------------------------------    ----------------------------------------
  node = cleavage point          <->   a column boundary (candidate cut site)
  edge = fragment                <->   a column range [a, b)
  edge cost = nt to encode div.  <->   block_report(): nt to encode the block
  shortest path over edges       <->   enumerate the 2 cuts (K=3 special case;
                                       becomes the SwiftLib/RASPP DP for K>3)
  rainbow / fidelity threshold   <->   cut_cost(): the pluggable chemistry hook
  degenerate codon = free width  <->   per-block degenerate annotation (reuses
                                       greedy_oligo.best_codon_for)

KEY DIFFERENCE (our added contribution): GGAssembler takes the desired diversity
as GIVEN (a resfile: which amino acids at which positions) and only minimizes $.
It does NOT decide WHAT to put in each fragment and does NOT model covariation
between positions, so it cannot reduce phantom-recombinant JUNK by cut placement.
Here the diversity is DERIVED from the natural cluster, and the objective is to
place the 2 cuts so the combinatorial product (the junk) is smallest at full
coverage -- i.e. keep linked columns inside one fragment, cut only where natural
variation is ~independent. The cheap-$ encoding step is then exactly what you
hand to GGAssembler as the back end.

So the intended pipeline is:
    natural cores  --[this script: choose cuts + per-block content]-->
    per-fragment target diversity  --[GGAssembler: cheapest oligos + overhangs]-->
    physical order.

----------------------------------------------------------------------------- #
ENCODING CHOICE in this round-1 version
----------------------------------------------------------------------------- #
Each fragment is encoded DISCRETELY: we synthesize the distinct observed
sub-sequences for that column range. This (a) captures within-fragment linkage
exactly (no phantom recombinants inside a fragment), (b) gives 100% coverage of
the input cores by construction, and (c) handles LENGTH VARIATION for free -- a
sub-sequence of a different length (an indel/loop) is just another distinct
variant in that layer; because cuts land between alignment columns, every
variant in a layer shares the same end columns and plugs into the same
neighbours regardless of its interior length.

Junk (producible non-natural recombinants) = (product of per-block variant
counts) - (number of distinct natural cores). Minimizing the product over the 2
cuts is the covariation-aware objective. Blocks whose variants are all the same
length can later be COMPRESSED with a degenerate codon (GGAssembler's free-width
trick); we annotate that but keep the discrete count as the junk metric.

Run:
  # cores must be ALIGNED first (gaps as '-'); e.g. in WSL:
  #   mafft --auto cluster1.core.fasta > cluster1.core.aln.fasta
  py fragment_split.py cluster1.core.aln.fasta
  py fragment_split.py cluster1.core.aln.fasta --min-block-cols 25 --chemistry gg
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict

# Reuse the single-oligo prototype as the PER-FRAGMENT encoder: within each
# chosen fragment we run its greedy degenerate set-cover -- but in LOSSLESS mode
# (cover everything; the budget only decides what is degenerate-compressed vs
# added as a discrete oligo, nothing is dropped). That is what keeps it
# consistent with the impossibility result: the per-fragment greedy makes no
# independent coverage/dropping decision, so whole-core coverage stays 100%.
import json
import os
from datetime import datetime

try:
    from greedy_oligo import (
        best_codon_for, _PLAIN, design_oligo,
        expand_degenerate_codon, AA_BY_CODON,
    )
    _HAVE_GREEDY = True
except Exception:  # pragma: no cover - falls back to discrete-only encoding
    best_codon_for = None
    _PLAIN = set("ACGT")
    design_oligo = None
    _HAVE_GREEDY = False


# --------------------------------------------------------------------------- #
# 1. Input: an ALIGNED core FASTA (gaps '-'), headers like >coreN_n<k>
# --------------------------------------------------------------------------- #

def read_aligned_cores(path):
    """Return [(aligned_seq, weight), ...]; all aligned_seq must be equal length."""
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
            "Align the cores first (gaps become '-'), e.g. in WSL:\n"
            f"  mafft --auto <cores.fasta> > <cores.aln.fasta>\n"
            "then re-run this script on the .aln.fasta."
        )
    return seqs


# --------------------------------------------------------------------------- #
# 2. Per-block content = the distinct ungapped sub-sequences over a column range
# --------------------------------------------------------------------------- #

def block_variants(seqs, a, b):
    """
    For columns [a, b) collapse each core to its real (ungapped) fragment and
    group. Returns dict {fragment_aa_seq: total_weight} -- the discrete library
    for that fragment layer.
    """
    variants = defaultdict(int)
    for s, w in seqs:
        frag = s[a:b].replace("-", "")
        variants[frag] += w
    return dict(variants)


def block_distinct_count(seqs, a, b):
    """CHEAP score used during cut enumeration: number of distinct ungapped
    sub-sequences in columns [a, b). The PRODUCT of these across the 3 fragments
    is the producible-library lower bound (= the covariation/junk signal that
    decides where to cut). Kept light so the O(L^2) enumeration stays fast; the
    expensive greedy encoding runs only on the winning split."""
    return len(block_variants(seqs, a, b))


# --------------------------------------------------------------------------- #
# 2b. Global budgeted greedy: pick WHICH whole cores to cover under a hard cap
#     on library size (the "junk budget"), dropping the rest. This is the direct
#     analogue of greedy_oligo's budgeted set-cover, but the budget is the
#     producible library size (junk) instead of the number of degenerate bases.
# --------------------------------------------------------------------------- #

def select_cores_under_budget(seqs, c1, c2, L, junk_budget):
    """
    Greedily grow a set S of whole cores to cover, keeping the producible library
    size (product over the 3 fragments of the number of DISTINCT pieces in S)
    at or below `junk_budget`. A core whose three pieces are all already present
    costs ZERO extra junk (a free recombinant); a core that introduces new pieces
    multiplies the library. We add the feasible core with the best weight gained
    per multiplicative junk spent, and stop when nothing more fits under the cap.

    Multi-starts from every seed core (n small) and keeps the run with the most
    covered weight. Because S is one consistent set of cores, every core in S has
    all three of its pieces present simultaneously -> it really is assembled
    (this is the global, coordinated dropping; per-fragment dropping is the trap
    we must avoid).
    """
    pieces = [(s[0:c1].replace("-", ""),
               s[c1:c2].replace("-", ""),
               s[c2:L].replace("-", "")) for s, _ in seqs]
    weights = [w for _, w in seqs]
    n = len(seqs)

    def run_from(seed):
        present = [set(), set(), set()]
        in_s = [False] * n
        for f in range(3):
            present[f].add(pieces[seed][f])
        in_s[seed] = True
        covered_w = weights[seed]
        while True:
            base = len(present[0]) * len(present[1]) * len(present[2])
            best = None  # (score, gain, -p, c)
            for c in range(n):
                if in_s[c]:
                    continue
                p = 1
                for f in range(3):
                    p *= len(present[f]) + (0 if pieces[c][f] in present[f] else 1)
                if p > junk_budget:
                    continue
                gain = weights[c]
                score = gain * base / p          # weight per multiplicative junk
                key = (score, gain, -p, -c)
                if best is None or key > best[0]:
                    best = (key, c, p)
            if best is None:
                break
            _, c, _ = best
            in_s[c] = True
            for f in range(3):
                present[f].add(pieces[c][f])
            covered_w += weights[c]
        S = [i for i in range(n) if in_s[i]]
        lib = len(present[0]) * len(present[1]) * len(present[2])
        return {"S": S, "present": present, "covered_w": covered_w, "lib": lib}

    best = None
    for seed in range(n):
        r = run_from(seed)
        key = (r["covered_w"], -r["lib"], len(r["S"]))
        if best is None or key > best[0]:
            best = (key, r)
    return best[1]


def encode_fragment(pieces, max_degenerate):
    """
    Turn ONE fragment's chosen distinct pieces into orderable oligos. With
    max_degenerate == 0 (default) every distinct piece is its own discrete oligo,
    so the fragment's library contribution = number of distinct pieces and the
    junk budget stays exactly honest. With max_degenerate > 0 we try to
    degenerate-compress equal-length pieces (fewer oligos) and report the codon
    expansion so the caller can check it still fits the budget.
    """
    by_len = defaultdict(list)
    for p in pieces:
        by_len[len(p)].append(p)

    oligos, lib, nt, deg_bases = [], 0, 0, 0
    for width, rows in sorted(by_len.items()):
        if len(rows) == 1 or max_degenerate <= 0 or not _HAVE_GREEDY:
            for r in rows:
                oligos.append({"kind": "discrete", "aa": r})
                nt += width * 3
                lib += 1
            continue
        res = design_oligo([(r, 1) for r in rows], max_degenerate=max_degenerate)
        oligo = res["oligo"]
        codons = [oligo[i:i + 3] for i in range(0, len(oligo), 3)]
        aa_sets = [frozenset(AA_BY_CODON[t] for t in expand_degenerate_codon(c))
                   for c in codons]
        covered = {r for r in rows if all(r[j] in aa_sets[j] for j in range(width))}
        oligos.append({"kind": "degenerate", "dna": oligo,
                       "ndeg": res["degenerate_bases"],
                       "lib": res["library_size"], "covers": sorted(covered)})
        nt += len(oligo)
        deg_bases += res["degenerate_bases"]
        lib += res["library_size"]
        for r in rows:
            if r not in covered:
                oligos.append({"kind": "discrete", "aa": r})
                nt += width * 3
                lib += 1
    return {"oligos": oligos, "lib": lib, "nt": nt, "deg_bases": deg_bases,
            "n_oligos": len(oligos)}


# --------------------------------------------------------------------------- #
# 3. The pluggable chemistry hook -- the ONLY thing GG vs HR changes.
#    No assumptions baked in: default is agnostic (every boundary is legal).
# --------------------------------------------------------------------------- #

def constant_run_around(seqs, p, want):
    """How many consecutive constant, gap-free alignment columns straddle the
    boundary at p (i.e. columns ... p-1 | p ...). Used by the chemistry hooks to
    decide whether a real overhang / homology arm can sit there."""
    n = len(seqs)

    def is_const(j):
        col = {s[j] for s, _ in seqs}
        return len(col) == 1 and "-" not in col

    left = 0
    j = p - 1
    while j >= 0 and is_const(j):
        left += 1
        j -= 1
    right = 0
    j = p
    L = len(seqs[0][0])
    while j < L and is_const(j):
        right += 1
        j += 1
    return left, right


def make_cut_cost(chemistry):
    """
    Returns cut_cost(seqs, p) -> finite penalty or float('inf') (illegal here).
    'agnostic' makes NO assumption about where cuts may go (default); the
    objective's junk term alone steers cuts toward low-linkage columns. 'gg' and
    'hr' add the physical anchor requirement as a soft/hard cost instead of a
    pre-filter, so the same optimizer serves both chemistries.
    """
    if chemistry == "agnostic":
        return lambda seqs, p: 0.0

    if chemistry == "gg":
        # Golden Gate: a fixed 4-nt overhang needs ~2 constant codons straddling
        # the cut. (Overhang-orthogonality / the limited high-fidelity set is
        # GGAssembler's 'rainbow' constraint -- left for the back end.)
        def cc(seqs, p):
            left, right = constant_run_around(seqs, p)
            return 0.0 if (left >= 1 and right >= 1) else float("inf")
        return cc

    if chemistry == "hr":
        # Homologous recombination / Gibson: needs a longer shared, uniquifiable
        # constant window (~20 bp ~= 7 codons total straddling the cut).
        def cc(seqs, p):
            left, right = constant_run_around(seqs, p)
            return 0.0 if (left + right >= 7 and left >= 2 and right >= 2) else float("inf")
        return cc

    raise ValueError(f"unknown chemistry {chemistry!r}")


# --------------------------------------------------------------------------- #
# 4. Enumerate the 2 cuts (K=3 special case of the RASPP/SwiftLib shortest path)
# --------------------------------------------------------------------------- #

def choose_cuts(seqs, min_block_cols, chemistry):
    """PHASE A -- pick the 2 cuts. Enumerate every legal boundary pair and keep
    the split that minimizes the product of per-fragment distinct-variant counts
    (the producible-library lower bound = the covariation/junk objective). This
    is the K=3 special case of the RASPP/SwiftLib shortest path; for K>3 it
    becomes the Bellman DP over the same boundary->fragment graph."""
    L = len(seqs[0][0])
    cut_cost = make_cut_cost(chemistry)

    legal = [p for p in range(min_block_cols, L - min_block_cols + 1)
             if cut_cost(seqs, p) != float("inf")]

    best = None  # (key, c1, c2)
    for c1 in legal:
        if c1 > L - 2 * min_block_cols:
            continue
        n1 = block_distinct_count(seqs, 0, c1)
        for c2 in legal:
            if c2 < c1 + min_block_cols or c2 > L - min_block_cols:
                continue
            n2 = block_distinct_count(seqs, c1, c2)
            n3 = block_distinct_count(seqs, c2, L)
            product = n1 * n2 * n3
            balance = max(c1, c2 - c1, L - c2) - min(c1, c2 - c1, L - c2)
            key = (product, n1 + n2 + n3, balance)
            if best is None or key < best[0]:
                best = (key, c1, c2)

    if best is None:
        sys.exit("no legal 2-cut split under the given constraints "
                 "(try a smaller --min-block-cols or --chemistry agnostic).")
    return best[1], best[2], L


def design_three_fragments(seqs, min_block_cols, chemistry, junk_budget, max_degenerate):
    """PHASE A (choose cuts) + PHASE B (budgeted greedy: cover as much weight as
    fits under the junk budget, dropping the rest)."""
    c1, c2, L = choose_cuts(seqs, min_block_cols, chemistry)
    n_unique_cores = len(seqs)
    total_w = sum(w for _, w in seqs)

    sel = select_cores_under_budget(seqs, c1, c2, L, junk_budget)
    present = sel["present"]
    cuts = [(0, c1), (c1, c2), (c2, L)]

    # Encode each fragment's chosen pieces. Default (max_degenerate=0) = discrete,
    # so library == product of distinct piece counts and the budget is exact.
    # If compression is requested, accept it only while the realized library
    # stays under the junk budget (compression spends junk to save oligos).
    frags = []
    realized_lib = 1
    for (a, b), pset in zip(cuts, present):
        pieces = sorted(pset)
        enc = encode_fragment(pieces, max_degenerate)
        if max_degenerate > 0:
            disc_lib = len(pieces)
            # would committing this fragment's compression bust the budget?
            tentative = realized_lib
            for other in present[len(frags) + 1:]:
                tentative *= len(other)        # remaining frags assumed discrete
            if tentative * enc["lib"] > junk_budget:
                enc = encode_fragment(pieces, 0)   # revert to discrete
        lengths = sorted({len(p) for p in pieces})
        frags.append({
            "a": a, "b": b, "ncols": b - a,
            "n_pieces": len(pieces),
            "lengths": lengths,
            "has_indel": len(lengths) > 1,
            **enc,
        })
        realized_lib *= enc["lib"]

    covered_w = sel["covered_w"]
    return {
        "c1": c1, "c2": c2, "L": L,
        "frags": frags,
        "junk_budget": junk_budget,
        "library_size": realized_lib,
        "n_unique_cores": n_unique_cores,
        "n_covered_cores": len(sel["S"]),
        "covered_weight": covered_w,
        "total_weight": total_w,
        "coverage_pct": 100.0 * covered_w / total_w,
        "junk": realized_lib - len(sel["S"]),
        "total_oligos": sum(f["n_oligos"] for f in frags),
        "total_nt": sum(f["nt"] for f in frags),
    }


# --------------------------------------------------------------------------- #
# 5. CLI
# --------------------------------------------------------------------------- #

def build_report(seqs, args, res):
    """Human-readable report (also what gets printed)."""
    total_w = sum(w for _, w in seqs)
    lines = []
    lines.append(f"input: {args.aln_fasta}")
    lines.append(f"loaded {len(seqs)} unique aligned cores ({total_w} natural sequences), "
                 f"alignment width {len(seqs[0][0])} columns")
    lines.append(f"chemistry: {args.chemistry}   min block cols: {args.min_block_cols}   "
                 f"junk budget (max library): {res['junk_budget']:,}   "
                 f"degenerate budget: {args.max_degenerate}")
    lines.append("")
    lines.append(f"best 2-cut split: columns [0,{res['c1']}) | "
                 f"[{res['c1']},{res['c2']}) | [{res['c2']},{res['L']})")
    lines.append("")
    for i, f in enumerate(res["frags"], 1):
        lvar = (f"{f['lengths'][0]}-{f['lengths'][-1]} aa  <-- LENGTH VARIATION (indel)"
                if f["has_indel"] else f"{f['lengths'][0]} aa")
        lines.append(f"  fragment {i}: cols [{f['a']},{f['b']})  width {f['ncols']}")
        lines.append(f"      {f['n_pieces']} distinct pieces kept, length {lvar}")
        lines.append(f"      encoded as {f['n_oligos']} oligos, "
                     f"{f['deg_bases']} degenerate nt, ~{f['nt']} nt")
    lines.append("")
    lines.append(f"COVERAGE (cores whose 3 pieces are all in the library):")
    lines.append(f"  covered:   {res['n_covered_cores']}/{res['n_unique_cores']} unique cores  "
                 f"= {res['covered_weight']}/{res['total_weight']} natural seqs "
                 f"({res['coverage_pct']:.0f}%)")
    lines.append(f"  dropped:   {res['n_unique_cores'] - res['n_covered_cores']} cores "
                 f"(didn't fit under the junk budget)")
    lines.append("")
    lines.append(f"producible library size: {res['library_size']:,}  "
                 f"(budget {res['junk_budget']:,})")
    lines.append(f"phantom recombinants (junk): {res['junk']:,}")
    lines.append(f"total oligos to order: {res['total_oligos']}")
    lines.append(f"total synthesis (GGAssembler-style): ~{res['total_nt']} nt")
    lines.append("")
    lines.append("next step: hand each fragment's kept pieces to GGAssembler as the "
                 "back end to pick the cheapest degenerate oligos + orthogonal overhangs.")
    return "\n".join(lines)


def save_run(out_root, stem, chemistry, seqs, args, res, report):
    """Write this run to algoruns/<timestamp>_<stem>_<chemistry>/."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(out_root, f"{ts}_{stem}_{chemistry}")
    os.makedirs(run_dir, exist_ok=True)

    with open(os.path.join(run_dir, "report.txt"), "w") as fh:
        fh.write(report + "\n")

    summary = {
        "input": os.path.abspath(args.aln_fasta),
        "n_unique_cores": res["n_unique_cores"],
        "n_natural_sequences": res["total_weight"],
        "alignment_width": len(seqs[0][0]),
        "chemistry": chemistry,
        "min_block_cols": args.min_block_cols,
        "junk_budget": res["junk_budget"],
        "max_degenerate": args.max_degenerate,
        "cuts": [res["c1"], res["c2"]],
        "library_size": res["library_size"],
        "junk": res["junk"],
        "n_covered_cores": res["n_covered_cores"],
        "covered_weight": res["covered_weight"],
        "total_weight": res["total_weight"],
        "coverage_pct": res["coverage_pct"],
        "total_oligos": res["total_oligos"],
        "total_nt": res["total_nt"],
        "fragments": res["frags"],
    }
    with open(os.path.join(run_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # one FASTA per fragment listing the KEPT distinct pieces (what you order).
    for i, f in enumerate(res["frags"], 1):
        pieces = []
        for o in f["oligos"]:
            if o["kind"] == "discrete":
                pieces.append(o["aa"])
            else:                                  # degenerate: list what it covers
                pieces.extend(o.get("covers", []))
        with open(os.path.join(run_dir, f"fragment{i}.fasta"), "w") as fh:
            for k, frag in enumerate(sorted(set(pieces)), 1):
                fh.write(f">frag{i}_p{k}_len{len(frag)}\n{frag}\n")

    return run_dir


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aln_fasta", help="ALIGNED core FASTA (gaps '-'), headers >coreN_n<k>")
    ap.add_argument("--min-block-cols", type=int, default=20,
                    help="minimum fragment width in alignment columns (default 20)")
    ap.add_argument("--chemistry", choices=["agnostic", "gg", "hr"], default="agnostic",
                    help="cut-site feasibility model; 'agnostic' (default) makes no "
                         "assumption about where cuts may go")
    ap.add_argument("--junk-budget", type=int, default=10000,
                    help="HARD cap on producible library size; cores that don't fit "
                         "under it are dropped (default 10000). Lower = fewer oligos, "
                         "less junk, less coverage.")
    ap.add_argument("--max-degenerate", type=int, default=0,
                    help="per-fragment degenerate-nt budget for optional oligo "
                         "compression (default 0 = pure discrete, exact junk budget); "
                         "compression is applied only while it stays under --junk-budget")
    ap.add_argument("--out-dir", default="algoruns",
                    help="parent folder for per-run output subfolders (default algoruns/)")
    args = ap.parse_args()

    seqs = read_aligned_cores(args.aln_fasta)
    res = design_three_fragments(seqs, args.min_block_cols, args.chemistry,
                                 args.junk_budget, args.max_degenerate)
    report = build_report(seqs, args, res)
    print(report)

    stem = os.path.splitext(os.path.basename(args.aln_fasta))[0]
    run_dir = save_run(args.out_dir, stem, args.chemistry, seqs, args, res, report)
    print(f"\nsaved run to {run_dir}/  (report.txt, summary.json, fragment[1-3].fasta)")


if __name__ == "__main__":
    main()
