"""sample_library.py -- pull a foldable sample out of a finished design.

A design does not encode 54 proteins, it encodes a LIBRARY.  Cluster 1 at K=3
builds 2,116 distinct proteins out of 61 oligos: 54 are the natural cores the
greedy set out to cover, and 2,062 are recombinants nobody asked for -- the
"junk" the library-size cap is charged for.  Nothing so far has looked at whether
those recombinants are plausible proteins.

This enumerates the whole library from the emitted oligos and samples a set to
fold, defaulting to 5 covered cores and 95 junk.  The 5 are the control: if they
fold and the junk does not, that is a real signal rather than a folding artefact.

Junk comes in two kinds, and the manifest labels each one:

  codon variant   the three oligos are exactly the ones some natural core uses,
                  but a degenerate codon went the other way.  Differs from a real
                  core at one or two residues.
  chimera         the three oligos come from different cores, so the protein is a
                  recombinant that no natural sequence spans.

Sequences are the EXPRESSED proteins: the designed core plus the four residues
the backbone contributes so its overhangs sit in frame (GC+CGGA reads Ala-Gly,
GGTG+GC reads Gly-Gly).  --no-flanks folds the bare core instead.

Read-only.  Nothing here touches the design.

    python sample_library.py <run-dir> [--covered 5] [--junk 95] [--seed 0]

Writes into <run-dir>/fold/:
    sample_proteins.fasta   what to fold, bare ids in the first token
    sample_dna.fasta        the coding sequence each protein comes from
    sample_manifest.csv     provenance, class, nearest natural core and identity
    sample_report.txt       what the library holds and what was drawn from it
"""

import argparse
import csv
import itertools
import json
import os
import random
import sys

from _paths import HERE                      # noqa: F401  (puts HERE on sys.path)
import cutsearch_design as U


# =========================================================================== #
# Reading
# =========================================================================== #

def read_named_cores(path):
    """[(name, aligned_seq)].  read_aligned_cores drops the header, and here the
    header is the answer to "which core is this"."""
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and buf:
                    out.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line.upper())
    if name is not None and buf:
        out.append((name, "".join(buf)))
    return out


def unit_options(unit):
    """Every protein this one oligo can make: (variant_index, aa_tuple).

    `encodes` is the amino-acid set per column AFTER codon choice -- what the
    oligo actually produces, not what was asked for -- so the product over
    columns is exactly this oligo's share of the library."""
    return list(enumerate(itertools.product(*unit["encodes"])))


# =========================================================================== #
# The library
# =========================================================================== #

def aligned_string(L, parts):
    """An assembly laid back onto the alignment, so it can be compared to a core
    column by column with no alignment step."""
    row = ["-"] * L
    for cols, aas in parts:
        for col, aa in zip(cols, aas):
            row[col] = aa
    return "".join(row)


def identity(a, b):
    """Fraction of occupied columns where two aligned rows agree.  Columns where
    both are gaps are not evidence of anything and are excluded."""
    same = occupied = 0
    for x, y in zip(a, b):
        if x == "-" and y == "-":
            continue
        occupied += 1
        if x == y:
            same += 1
    return same / occupied if occupied else 0.0


def cores_per_unit(unit, cores, lo, hi):
    """Which natural cores this oligo can make, within its own fragment position.

    A core is only buildable from a unit when the GAP PATTERN matches: a
    degenerate codon encodes "V or I", never "residue or nothing"."""
    cols = tuple(unit["cols"])
    hit = set()
    for i, (_name, seq) in enumerate(cores):
        block = seq[lo:hi]
        core_cols = tuple(lo + j for j, c in enumerate(block) if c != "-")
        if core_cols != cols:
            continue
        if all(seq[c] in e for c, e in zip(cols, unit["encodes"])):
            hit.add(i)
    return hit


def build_library(frags, bounds, cores, L):
    """Every distinct protein the design can produce, with its provenance.

    Returns (rows, by_aligned) where a row is
        (aligned_seq, [(fragment position, oligo_index, variant_index)], single_core_set)
    and `single_core_set` is the cores whose oligos this assembly reuses in ALL
    fragment positions -- empty for a chimera."""
    fragment_positions = []
    for f, units in enumerate(frags):
        lo, hi = bounds[f], bounds[f + 1]
        opts = []
        for ui, u in enumerate(units):
            owned = cores_per_unit(u, cores, lo, hi)
            for vi, aas in unit_options(u):
                opts.append((ui, vi, tuple(u["cols"]), aas, owned))
        fragment_positions.append(opts)

    rows = []
    for combo in itertools.product(*fragment_positions):
        parts = [(cols, aas) for _ui, _vi, cols, aas, _own in combo]
        owned = set(combo[0][4])
        for c in combo[1:]:
            owned &= c[4]
        rows.append((aligned_string(L, parts),
                     [(f, c[0], c[1]) for f, c in enumerate(combo)],
                     owned))
    return rows


def dna_for(frags, provenance, frags_opts):
    """Concrete coding sequence for one assembly, codon by codon.

    The oligo carries an IUPAC triplet per column; `concrete_codon` picks one
    expansion of it that encodes the residue this variant chose."""
    out = []
    for (f, ui, vi) in provenance:
        unit = frags[f][ui]
        aas = frags_opts[f][ui][vi]
        oligo = unit["oligo"]
        for j, aa in enumerate(aas):
            trip = oligo[3 * j:3 * j + 3]
            cod = U.concrete_codon(trip, aa)
            if cod is None:
                raise AssertionError(
                    "oligo %d/%d column %d cannot encode %s" % (f + 1, ui, j, aa))
            out.append(cod)
    return "".join(out)


def translate(dna):
    return "".join(U.AA_BY_CODON[dna[i:i + 3]] for i in range(0, len(dna), 3))


# =========================================================================== #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a finished cutsearch_design.py run directory")
    ap.add_argument("--covered", type=int, default=5,
                    help="how many natural cores to include as the control")
    ap.add_argument("--junk", type=int, default=95,
                    help="how many recombinants to include")
    ap.add_argument("--seed", type=int, default=0,
                    help="sampling seed; recorded in the report")
    ap.add_argument("--prefix", default="",
                    help="prepended to every sequence id, e.g. 'c1_'")
    ap.add_argument("--backbone5", default="CGGA",
                    help="overhang at the backbone->fragment 1 junction")
    ap.add_argument("--backbone3", default="GGTG",
                    help="overhang at the fragment K->backbone junction")
    ap.add_argument("--backbone-pre", default="GC",
                    help="nt the backbone supplies immediately before backbone5, "
                         "putting its overhang in frame")
    ap.add_argument("--backbone-post", default="GC",
                    help="nt the backbone supplies immediately after backbone3")
    ap.add_argument("--no-flanks", action="store_true",
                    help="fold the designed core alone, without the residues the "
                         "backbone contributes at each end")
    ap.add_argument("--aln", default=None,
                    help="the alignment (default: the path in summary.json)")
    ap.add_argument("--out", default=None, help="output dir (default <run-dir>/fold)")
    args = ap.parse_args()

    with open(os.path.join(args.run_dir, "summary.json")) as fh:
        S = json.load(fh)
    frags = S["fragments"]
    cuts = S.get("cuts", [])

    aln = args.aln or S["input"]
    if not os.path.exists(aln):
        sys.exit("alignment not found: %s (pass --aln)" % aln)
    cores = read_named_cores(aln)
    L = len(cores[0][1])
    bounds = [0] + list(cuts) + [L]
    out_dir = args.out or os.path.join(args.run_dir, "fold")
    os.makedirs(out_dir, exist_ok=True)

    # --- the whole library ------------------------------------------------- #
    frags_opts = [[dict(unit_options(u)) for u in units] for units in frags]
    rows = build_library(frags, bounds, cores, L)

    expected = len(rows)
    rec = [r for r in S["frontier"] if r["K"] == S["recommended_K"]][0]
    assert expected == rec["library"], (
        "enumerated %d proteins, summary.json says the library is %d"
        % (expected, rec["library"]))

    by_core = {seq: name for name, seq in cores}
    covered, junk = [], []
    for aligned, prov, owned in rows:
        if aligned in by_core:
            covered.append((by_core[aligned], aligned, prov, owned))
        else:
            junk.append((aligned, prov, owned))

    assert len(covered) == rec["n_cores_encoded"], (
        "%d assemblies match a natural core, the design claims %d"
        % (len(covered), rec["n_cores_encoded"]))

    n_variant = sum(1 for _a, _p, owned in junk if owned)
    n_chimera = len(junk) - n_variant

    # --- draw the sample --------------------------------------------------- #
    rng = random.Random(args.seed)
    n_cov = min(args.covered, len(covered))
    n_junk = min(args.junk, len(junk))
    if n_cov < args.covered or n_junk < args.junk:
        print("  ! asked for %d + %d, the library only holds %d + %d"
              % (args.covered, args.junk, len(covered), len(junk)))
    pick_cov = rng.sample(covered, n_cov)
    pick_junk = rng.sample(junk, n_junk)

    # --- materialise ------------------------------------------------------- #
    records = []
    for name, aligned, prov, _owned in sorted(pick_cov):
        records.append(("%s%s" % (args.prefix, name.split("_")[0]),
                        "covered", name, aligned, prov, set()))
    for i, (aligned, prov, owned) in enumerate(pick_junk, start=1):
        kind = "codon_variant" if owned else "chimera"
        records.append(("%sjunk%03d" % (args.prefix, i), kind, "", aligned,
                        prov, owned))

    # What the ribosome actually sees.  The backbone supplies two nucleotides
    # either side so that its overhangs land in frame: GC + CGGA reads GCC GGA,
    # and GGTG + GC reads GGT GGC.  Those four residues are on the expressed
    # protein whether or not the designer put them there, so they are folded
    # unless --no-flanks says otherwise.
    pre = "" if args.no_flanks else args.backbone_pre + args.backbone5
    post = "" if args.no_flanks else args.backbone3 + args.backbone_post
    assert (len(pre) + len(post)) % 3 == 0, (
        "flanks add %d nt, which shifts the frame" % (len(pre) + len(post)))

    manifest, prot_lines, dna_lines = [], [], []
    for ident, kind, core_name, aligned, prov, owned in records:
        core_protein = aligned.replace("-", "")
        dna = pre + dna_for(frags, prov, frags_opts) + post
        assert len(dna) % 3 == 0, "%s: %d nt is not a whole number of codons" % (
            ident, len(dna))
        protein = translate(dna)
        assert core_protein in protein, (
            "%s: the designed core is not inside the expressed protein" % ident)
        assert "*" not in protein, "%s: internal stop" % ident

        best_name, best_id = "", 0.0
        for cname, cseq in cores:
            v = identity(aligned, cseq)
            if v > best_id:
                best_name, best_id = cname, v

        oligos = "+".join("f%d:o%d.v%d" % (f + 1, ui + 1, vi + 1)
                          for f, ui, vi in prov)
        manifest.append({
            "id": ident,
            "class": kind,
            "natural_core": core_name,
            "length_aa": len(protein),
            "designed_core_aa": len(core_protein),
            "flank_aa": len(protein) - len(core_protein),
            "oligos": oligos,
            "nearest_core": best_name,
            "identity_to_nearest_pct": round(100 * best_id, 2),
            "reuses_core_oligos": ";".join(sorted(cores[i][0] for i in owned)),
            "dna_nt": len(dna),
        })
        prot_lines.append(">%s %s len=%d nearest=%s id=%.1f%%\n%s"
                          % (ident, kind, len(protein), best_name or "-",
                             100 * best_id, protein))
        dna_lines.append(">%s %s\n%s" % (ident, kind, dna))

    with open(os.path.join(out_dir, "sample_proteins.fasta"), "w") as fh:
        fh.write("\n".join(prot_lines) + "\n")
    with open(os.path.join(out_dir, "sample_dna.fasta"), "w") as fh:
        fh.write("\n".join(dna_lines) + "\n")
    with open(os.path.join(out_dir, "sample_manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)

    # --- report ------------------------------------------------------------ #
    ids = [m["identity_to_nearest_pct"] for m in manifest if m["class"] != "covered"]
    lens = [m["length_aa"] for m in manifest]
    rep = []
    rep.append("LIBRARY SAMPLE FOR FOLDING")
    rep.append("=" * 70)
    rep.append("run        : %s" % os.path.abspath(args.run_dir))
    rep.append("design     : K=%d, cuts %s, %d oligos"
               % (S["recommended_K"], cuts, rec["oligos"]))
    rep.append("seed       : %d" % args.seed)
    rep.append("")
    rep.append("THE LIBRARY")
    rep.append("  %6d distinct proteins the 61 oligos can build" % len(rows))
    rep.append("  %6d natural cores the design set out to cover" % len(covered))
    rep.append("  %6d junk recombinants, of which" % len(junk))
    rep.append("         %6d codon variants (a real core's oligos, one "
               "degenerate codon flipped)" % n_variant)
    rep.append("         %6d chimeras (oligos from different cores)" % n_chimera)
    rep.append("")
    rep.append("THE SAMPLE  (%d sequences)" % len(records))
    rep.append("  %d covered cores, as the control" % n_cov)
    for m in manifest:
        if m["class"] == "covered":
            rep.append("      %-14s %s, %d aa" % (m["id"], m["natural_core"],
                                                  m["length_aa"]))
    rep.append("  %d junk: %d chimeras, %d codon variants"
               % (n_junk,
                  sum(1 for m in manifest if m["class"] == "chimera"),
                  sum(1 for m in manifest if m["class"] == "codon_variant")))
    if ids:
        ids_sorted = sorted(ids)
        rep.append("      identity to the nearest natural core: "
                   "%.1f%% min, %.1f%% median, %.1f%% max"
                   % (ids_sorted[0], ids_sorted[len(ids_sorted) // 2],
                      ids_sorted[-1]))
    rep.append("      length %d-%d aa" % (min(lens), max(lens)))
    rep.append("")
    rep.append("CHECKS")
    rep.append("  [PASS] enumerated library = %d, the size the design reports"
               % len(rows))
    rep.append("  [PASS] %d assemblies reproduce a natural core exactly"
               % len(covered))
    rep.append("  [PASS] every sampled protein back-translates from its own DNA")
    rep.append("  [PASS] no internal stop codon in any sampled protein")
    rep.append("")
    if args.no_flanks:
        rep.append("NOTE  --no-flanks: these are the designed cores alone.  The "
                   "residues the")
        rep.append("      backbone contributes at each end are NOT here, so this "
                   "is not quite")
        rep.append("      the protein that gets expressed.")
    else:
        rep.append("NOTE  These are the EXPRESSED proteins: the designed core "
                   "plus what the")
        rep.append("      backbone contributes.  %s reads %s at the N terminus "
                   "and %s reads %s"
                   % (pre, translate(pre), post, translate(post)))
        rep.append("      at the C, so every sequence carries %d extra residues "
                   "the designer"
                   % (len(translate(pre)) + len(translate(post))))
        rep.append("      never chose.  --no-flanks folds the bare core instead.")
    text = "\n".join(rep)
    with open(os.path.join(out_dir, "sample_report.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)
    print("\nwrote %s/" % out_dir)


if __name__ == "__main__":
    main()
