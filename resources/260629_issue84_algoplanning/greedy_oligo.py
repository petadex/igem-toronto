#!/usr/bin/env python3
"""
Greedy degenerate-codon oligo design -- FIRST PROTOTYPE.

Given a cluster of protein cores (amino-acid sequences), design ONE degenerate
DNA oligo. When that oligo is synthesized as a degenerate library it expands to
many distinct DNA sequences; we want it to capture as much of the cluster's
natural variation as possible using only a small, fixed number of degenerate
nucleotide positions.

What this prototype does (and the simplifying assumptions it makes):

  1. Greedy degenerate-codon insertion under a HARD CAP on the number of
     degenerate nucleotide positions (default 10). Each MSA column gets one
     degenerate codon; making a base position degenerate (an IUPAC ambiguity
     code instead of a plain A/C/G/T) "spends" from the budget. We greedily
     spend the budget on the positions that capture the most variation.
  2. Nucleotide <-> amino-acid conversion (standard genetic code + IUPAC codes).

  ASSUMPTION (to be relaxed later with the fragment idea): preprocessing is done
  and the cores are NOT all the same length. We sidestep indels for now by
  designing over the LARGEST equal-length subset of cores, so every position is a
  clean alignment column. Cores of other lengths are dropped for this prototype.

Run:  py ninetypidorfs/greedy_oligo.py ninetypidorfs/cluster1.core.fasta
      py ninetypidorfs/greedy_oligo.py ninetypidorfs/cluster1.core.fasta --max-degenerate 10
"""

from __future__ import annotations

import argparse
from collections import defaultdict

# --------------------------------------------------------------------------- #
# 1. Genetic code + IUPAC degenerate-base tables
# --------------------------------------------------------------------------- #

_BASES = "TCAG"
# Standard genetic code (NCBI table 1), in TCAG x TCAG x TCAG codon order.
_AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
AA_BY_CODON = dict(zip(_CODONS, _AA_STRING))  # "ATG" -> "M", "TAA" -> "*"

# IUPAC ambiguity codes: symbol -> set of concrete bases it stands for.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_PLAIN = set("ACGT")  # symbols that are NOT degenerate (cost 0)


def expand_degenerate_codon(deg: str) -> list[str]:
    """A 3-symbol IUPAC codon -> list of the concrete DNA triplets it encodes."""
    s1, s2, s3 = (IUPAC[c] for c in deg)
    return [a + b + c for a in s1 for b in s2 for c in s3]


def n_degenerate_bases(deg: str) -> int:
    """How many of a codon's 3 positions are degenerate (an IUPAC ambiguity code)."""
    return sum(1 for b in deg if b not in _PLAIN)


# --------------------------------------------------------------------------- #
# 2. Precompute every degenerate codon
# --------------------------------------------------------------------------- #
# A degenerate codon = 3 IUPAC symbols => 15*15*15 = 3375 possibilities.
# For each: (codon, frozenset_of_amino_acids, n_degenerate_bases, n_triplets, has_stop)
#   - n_degenerate_bases  = what it costs from the budget (the thing we cap)
#   - n_triplets          = its contribution to library size (the "junk" measure)

_ALL_CODONS: list[tuple[str, frozenset, int, int, bool]] = []
for _c1 in IUPAC:
    for _c2 in IUPAC:
        for _c3 in IUPAC:
            _deg = _c1 + _c2 + _c3
            _triplets = expand_degenerate_codon(_deg)
            _aas = frozenset(AA_BY_CODON[t] for t in _triplets)
            _ALL_CODONS.append(
                (_deg, _aas, n_degenerate_bases(_deg), len(_triplets), "*" in _aas)
            )

# Rank by (degenerate bases, then library junk): the first match in a scan is the
# cheapest codon in budget terms, breaking ties toward the least junk.
_ALL_CODONS.sort(key=lambda r: (r[2], r[3]))

_best_codon_cache: dict[frozenset, tuple | None] = {}


def best_codon_for(required_aas: frozenset):
    """
    Cheapest stop-free degenerate codon whose amino-acid set covers ALL of
    `required_aas`, where "cheapest" = fewest degenerate bases, then least junk.
    Returns (codon, covered_aa_set, n_degenerate_bases, n_triplets) or None if no
    single stop-free codon can cover that set of amino acids.
    """
    if required_aas in _best_codon_cache:
        return _best_codon_cache[required_aas]
    result = None
    for deg, aas, ndeg, ntrip, has_stop in _ALL_CODONS:  # cheapest-first
        if has_stop:
            continue
        if required_aas <= aas:
            result = (deg, aas, ndeg, ntrip)
            break
    _best_codon_cache[required_aas] = result
    return result


# --------------------------------------------------------------------------- #
# 3. FASTA input + subset selection
# --------------------------------------------------------------------------- #

def read_cores(path: str) -> list[tuple[str, int]]:
    """
    Read a *.core.fasta. Headers look like `>core12_n3`; the `_n<k>` suffix is the
    number of natural sequences that collapsed onto this unique core (its weight).
    Returns [(sequence, weight), ...].
    """
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
                buf.append(line)
    flush()
    return seqs


def largest_equal_length_subset(seqs):
    """Pick the modal length (most weighted sequences) and return that subset."""
    by_len = defaultdict(list)
    for s, w in seqs:
        by_len[len(s)].append((s, w))
    best_len = max(by_len, key=lambda L: sum(w for _, w in by_len[L]))
    return best_len, by_len[best_len]


# --------------------------------------------------------------------------- #
# 4. Build the oligo for a chosen set of cores
# --------------------------------------------------------------------------- #
# An oligo that "includes" a set T of cores must, at every column, use a codon
# whose amino-acid set covers every residue that the cores in T show there.
# Columns where all of T agree cost 0 degenerate bases; disagreements cost more.
# The same oligo automatically also covers any *other* core whose residues all
# happen to fall inside the chosen codons (a free bonus -- or junk, depending).

def build_for_cores(seqs, members, fallback):
    """
    Codons for an oligo that includes cores `members` (a list of row indices).
    `fallback` is the consensus, used per column when `members` is empty so the
    empty design is just the consensus sequence (0 degenerate bases).
    Returns (codons, aa_sets, total_degenerate_bases, library_size) or None if
    some column's required residues can't be covered by one stop-free codon.
    """
    ncols = len(fallback)
    codons, aa_sets, ndeg, lib = [], [], 0, 1
    for j in range(ncols):
        req = frozenset(seqs[i][j] for i in members) if members else frozenset(fallback[j])
        bc = best_codon_for(req)
        if bc is None:
            return None
        deg, aas, nd, ntrip = bc
        codons.append(deg)
        aa_sets.append(aas)
        ndeg += nd
        lib *= ntrip
    return codons, aa_sets, ndeg, lib


# --------------------------------------------------------------------------- #
# 5. Greedy set-cover over cores (budget = number of degenerate nt positions)
# --------------------------------------------------------------------------- #

def design_oligo(subset, max_degenerate=10):
    """
    Design ONE degenerate oligo that covers as many whole natural cores as
    possible using at most `max_degenerate` degenerate nucleotide positions.

    This is a set-cover: a core is "covered" only if EVERY one of its residues
    is in the corresponding codon, so degenerate bases must be concentrated to
    complete specific cores, not spread to capture isolated point variants.

    Greedy: grow the included-core set, each step adding the not-yet-covered core
    with the best (newly-covered weight) / (extra degenerate bases) that stays
    under the cap. Because the first core added is free (0 degenerate bases), we
    multi-start from every possible seed core and keep the best run (n is small).
    """
    seqs = [s for s, _ in subset]
    weights = [w for _, w in subset]
    nseq = len(seqs)
    total_weight = sum(weights)
    ncols = len(seqs[0])

    consensus = []
    for j in range(ncols):
        col = defaultdict(int)
        for i in range(nseq):
            col[seqs[i][j]] += weights[i]
        consensus.append(max(col, key=lambda a: (col[a], a)))

    def covered_rows(aa_sets):
        return [i for i in range(nseq)
                if all(seqs[i][j] in aa_sets[j] for j in range(ncols))]

    def wt(rows):
        return sum(weights[i] for i in rows)

    def greedy_from(seed):
        members = [seed]
        built = build_for_cores(seqs, members, consensus)
        if built is None or built[2] > max_degenerate:
            return None
        codons, aa_sets, ndeg, lib = built
        covered = covered_rows(aa_sets)
        traj = [(ndeg, wt(covered))]
        while True:
            best = None  # (score, built, covered)
            for c in range(nseq):
                if c in covered:
                    continue
                trial = build_for_cores(seqs, members + [c], consensus)
                if trial is None or trial[2] > max_degenerate:
                    continue
                new_covered = covered_rows(trial[1])
                gain = wt(new_covered) - wt(covered)
                if gain <= 0:
                    continue
                added = trial[2] - ndeg
                score = gain / added if added > 0 else float("inf")
                if best is None or score > best[0]:
                    best = (score, c, trial, new_covered)
            if best is None:
                break
            _, c, trial, new_covered = best
            members.append(c)
            codons, aa_sets, ndeg, lib = trial
            covered = new_covered
            traj.append((ndeg, wt(covered)))
        return {"members": members, "codons": codons, "aa_sets": aa_sets,
                "ndeg": ndeg, "lib": lib, "covered": covered, "traj": traj}

    best_run = None
    for seed in range(nseq):
        run = greedy_from(seed)
        if run is None:
            continue
        key = (wt(run["covered"]), -run["ndeg"])  # most coverage, then fewest bases
        if best_run is None or key > best_run[0]:
            best_run = (key, run)
    run = best_run[1]

    oligo = "".join(run["codons"])
    degenerate_cols = sum(1 for c in run["codons"] if any(b not in _PLAIN for b in c))
    return {
        "oligo": oligo,
        "degenerate_bases": run["ndeg"],
        "degenerate_cols": degenerate_cols,
        "library_size": run["lib"],
        "whole_seq_coverage": wt(run["covered"]),
        "n_cores_covered": len(run["covered"]),
        "total_weight": total_weight,
        "ncols": ncols,
        "trajectory": run["traj"],
    }


# --------------------------------------------------------------------------- #
# 6. CLI
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="path to a *.core.fasta")
    ap.add_argument("--max-degenerate", type=int, default=10,
                    help="hard cap on degenerate nucleotide positions (default 10)")
    args = ap.parse_args()

    seqs = read_cores(args.fasta)
    best_len, subset = largest_equal_length_subset(seqs)
    sub_weight = sum(w for _, w in subset)
    print(f"loaded {len(seqs)} unique cores "
          f"({sum(w for _, w in seqs)} natural sequences)")
    print(f"designing over largest equal-length subset: length {best_len}, "
          f"{len(subset)} unique cores ({sub_weight} natural sequences)")
    print(f"degenerate-base budget: {args.max_degenerate}\n")

    res = design_oligo(subset, max_degenerate=args.max_degenerate)
    tot = res["total_weight"]

    print("greedy trajectory (degenerate bases spent -> cores covered):")
    print(f"  {'deg bases':>9}  {'cores covered':>14}")
    last = None
    for nd, wc in res["trajectory"]:
        if last is None or wc != last:
            print(f"  {nd:>9}  {wc:>3}/{tot:<3} ({100*wc/tot:>3.0f}%)")
            last = wc
    print()

    print(f"final oligo: {res['ncols']} codons, "
          f"{res['degenerate_bases']} degenerate nt positions "
          f"across {res['degenerate_cols']} codons")
    print(f"library size (junk): {res['library_size']:,}")
    print(f"whole-sequence coverage: {res['whole_seq_coverage']}/{tot} cores "
          f"({100*res['whole_seq_coverage']/tot:.0f}%)")
    print(f"\n{res['oligo']}")


if __name__ == "__main__":
    main()
