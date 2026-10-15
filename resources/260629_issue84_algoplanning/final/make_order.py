"""Turn a finished design into order-ready constructs, and prove they reassemble.

Adds the BsmBI landing pads, the spacers, the backbone overhangs and the frame
padding, then simulates digestion and ligation to check the result is exactly the
sequence the designer intended.  See POSTPROCESS_PLAN.md for the reasoning.

Per ordered oligo:

    CGTCTC [sp] [ ---- insert ---- ] [sp] GAGACG
                  ^ first 4 nt = left overhang
                            last 4 nt = right overhang ^

No universal primer tails: wet lab confirmed 2026-09-09 that the fragments are
not amplified, so a tail would be dead weight that is paid for and immediately
cut off.

Adjacent inserts OVERLAP by their shared 4-nt overhang: the internal overhangs
already live in the coding sequence, straddling a codon boundary, so each insert
reaches 2 nt into each neighbour.  After ligation the overlap merges and the
overhang appears once.

  python make_order.py <run-dir> [--out DIR]
"""

import argparse
import csv
import itertools
import json
import os
import sys

from _paths import HERE                      # noqa: F401  (puts HERE on sys.path)
import cutsearch_design as U

BSMBI = "CGTCTC"
BSMBI_RC = "GAGACG"
RARE = {"AGA", "AGG", "CGA", "CGG", "ATA", "CTA", "CCC", "GGA", "GGG", "TCA"}


# =========================================================================== #
# Construction
# =========================================================================== #

def pick_spacer(insert, choices="ACGT"):
    """One nucleotide either side of the insert, chosen so the finished construct
    contains no Type IIS site beyond the two intended ones.

    Biologically the spacer is arbitrary -- digestion discards it -- but it sits
    directly against the recognition site, so a careless choice can spell a THIRD
    site at the seam and the fragment would be cut in half."""
    for sp in choices:
        for sp2 in choices:
            cand = BSMBI + sp + insert + sp2 + BSMBI_RC
            if count_sites(cand) == 2:
                return sp, sp2
    return None, None


def count_sites(seq):
    """How many Type IIS sites some expansion of `seq` could contain.

    IUPAC-aware on purpose: a var2 oligo spells two sequences, and a site present
    in only one of them is still a site that will cut."""
    n = 0
    for site in (BSMBI, BSMBI_RC, "GGTCTC", "GAGACC"):
        m = len(site)
        for i in range(len(seq) - m + 1):
            if all(site[k] in U.IUPAC[seq[i + k]] for k in range(m)):
                n += 1
    return n


def build_inserts(design, five_pad, backbone5, backbone3):
    """The insert for every oligo, keyed by fragment index.

    Flanks are uniform within a fragment: the junction codons are PINNED, so
    every fragment-f oligo ends with the same two bases.  Taken from the junction
    tokens rather than by inspecting oligos, then asserted against the oligos."""
    frags = design["fragments"]
    tokens = design["junctions"]          # [(overhang, left_codon, right_codon)]
    K = len(frags)

    left_add, right_add = [], []
    for f in range(K):
        if f == 0:
            left_add.append(backbone5 + five_pad)     # CGGA + frame pad
        else:
            left_add.append(tokens[f - 1][1][-2:])    # last 2 nt of left codon
        if f == K - 1:
            right_add.append(backbone3)               # GGTG
        else:
            right_add.append(tokens[f][2][:2])        # first 2 nt of right codon

    # The pinned codons must actually be what the oligos carry.
    for f, units in enumerate(frags):
        if f > 0:
            got = {u["oligo"][:2] for u in frags[f - 1]}
            # previous fragment's last 2 nt must equal what we prepend here
            got = {u["oligo"][-2:] for u in frags[f - 1]}
            assert got == {left_add[f]}, (
                f"fragment {f} left flank {left_add[f]} does not match "
                f"fragment {f} oligo endings {got}")
        if f < K - 1:
            got = {u["oligo"][:2] for u in frags[f + 1]}
            assert got == {right_add[f]}, (
                f"fragment {f+1} right flank {right_add[f]} does not match "
                f"fragment {f+2} oligo starts {got}")

    out = []
    for f, units in enumerate(frags):
        out.append([left_add[f] + u["oligo"] + right_add[f] for u in units])
    return out, left_add, right_add


# =========================================================================== #
# Digestion / ligation -- the proof that the constructs work
# =========================================================================== #

def digest(construct, site_i, site_j):
    """The piece BsmBI releases.  CGTCTC(1/5): the top strand is cut 1 nt past
    the site and the bottom strand 5 nt past, so the released piece begins 7 nt
    after the left site starts and ends 1 nt before the right site starts."""
    return construct[site_i + 7:site_j - 1]


def ligate(pieces):
    """Join by overhang matching.  Each piece's last 4 nt must equal the next
    piece's first 4 -- that shared window is the overhang, and it survives once."""
    out = pieces[0]
    for p in pieces[1:]:
        if out[-4:] != p[:4]:
            return None, f"overhang mismatch: {out[-4:]} vs {p[:4]}"
        out += p[4:]
    return out, None


def translate(dna):
    aa = []
    for i in range(0, len(dna) - 2, 3):
        c = dna[i:i + 3]
        aa.append(U.AA_BY_CODON.get(c, "X"))
    return "".join(aa)


# =========================================================================== #
# Main
# =========================================================================== #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="a finished cutsearch_design.py run directory")
    ap.add_argument("--backbone5", default="CGGA",
                    help="overhang at the backbone->fragment 1 junction")
    ap.add_argument("--backbone3", default="GGTG",
                    help="overhang at the fragment K->backbone junction")
    ap.add_argument("--five-pad", default="",
                    help="nt inserted between backbone5 and the first codon. "
                         "EMPTY by default: wet lab confirmed 2026-09-09 that the "
                         "backbone already carries the frame offset, so fragment 1 "
                         "begins immediately after CGGA. Setting one nt here gives "
                         "the ?CG|GA? framing instead, which is the only one that "
                         "avoids a rare codon at this junction")
    ap.add_argument("--backbone-pre", default="GC",
                    help="nt the BACKBONE supplies immediately before backbone5. "
                         "Not part of our order; only its LENGTH mod 3 matters for "
                         "frame (must be 2 when --five-pad is empty). The bases "
                         "themselves are wet lab's and are assumed here purely so "
                         "the frame and translation checks have something to run on")
    ap.add_argument("--backbone-post", default="GC",
                    help="nt the BACKBONE supplies immediately after backbone3. "
                         "Not part of our order; needed to check frame")
    ap.add_argument("--spacer", default="auto",
                    help="the 1 nt between each BsmBI site and the overhang. "
                         "'auto' takes the first of A/C/G/T that leaves no extra "
                         "Type IIS site; give a base to match a house standard "
                         "(the iGEM spec map uses C). A forced spacer that would "
                         "create a third site is refused, not silently replaced")
    ap.add_argument("--max-roundtrip", type=int, default=2000,
                    help="cap on assemblies simulated (all, if fewer exist)")
    ap.add_argument("--out", default=None, help="output dir (default <run-dir>/order)")
    args = ap.parse_args()

    # Specs often write tails in lowercase and the functional parts in uppercase
    # (Benchling's convention).  Case is display only, so normalise before any
    # sequence comparison -- the IUPAC table is uppercase-keyed and a lowercase
    # base would otherwise raise mid-check rather than simply not matching.
    for k in ("backbone5", "backbone3", "five_pad",
              "backbone_pre", "backbone_post"):
        setattr(args, k, getattr(args, k).upper())

    with open(os.path.join(args.run_dir, "summary.json")) as fh:
        d = json.load(fh)
    frags = d["fragments"]
    K = len(frags)
    outdir = args.out or os.path.join(args.run_dir, "order")
    os.makedirs(outdir, exist_ok=True)

    inserts, left_add, right_add = build_inserts(
        d, args.five_pad, args.backbone5, args.backbone3)

    # ---- build every construct ------------------------------------------- #
    rows, constructs, problems = [], [], []
    for f in range(K):
        choices = "ACGT" if args.spacer == "auto" else args.spacer.upper()
        sp, sp2 = pick_spacer(inserts[f][0], choices)
        if sp is None:
            problems.append(
                f"fragment {f+1}: spacer "
                + ("'auto' found no base of A/C/G/T that leaves the construct "
                   "free of extra Type IIS sites"
                   if args.spacer == "auto" else
                   f"'{args.spacer}' creates an extra Type IIS site -- pick "
                   f"another base or use --spacer auto"))
            continue
        # Spell the additions out in words, per fragment, so a reader of the CSV
        # can see what was bolted onto the designed oligo without decoding it.
        la, ra = left_add[f], right_add[f]
        if f == 0:
            l_desc = (f"{BSMBI}(BsmBI) + {sp}(spacer) + "
                      f"{args.backbone5}(backbone overhang)"
                      + (f" + {args.five_pad}(frame pad)" if args.five_pad else ""))
        else:
            l_desc = (f"{BSMBI}(BsmBI) + {sp}(spacer) + {la}"
                      f"(last 2nt of fragment {f} to complete overhang "
                      f"{d['junctions'][f-1][0]})")
        if f == K - 1:
            r_desc = (f"{args.backbone3}(backbone overhang) + {sp2}(spacer) + "
                      f"{BSMBI_RC}(BsmBI, reversed)")
        else:
            r_desc = (f"{ra}(first 2nt of fragment {f+2} to complete overhang "
                      f"{d['junctions'][f][0]}) + {sp2}(spacer) + "
                      f"{BSMBI_RC}(BsmBI, reversed)")

        for i, ins in enumerate(inserts[f], start=1):
            u = frags[f][i - 1]
            prefix = BSMBI + sp + la
            suffix = ra + sp2 + BSMBI_RC
            con = prefix + u["oligo"] + suffix
            site_i = 0
            site_j = 6 + 1 + len(ins) + 1
            assert con[site_i:site_i + 6] == BSMBI
            assert con[site_j:site_j + 6] == BSMBI_RC
            assert con == BSMBI + sp + ins + sp2 + BSMBI_RC
            name = f"frag{f+1}_oligo{i}_var{u['variants']}"
            constructs.append((f, name, con, site_i, site_j, ins))
            rows.append({
                "name": name,
                "fragment": f + 1,
                "oligo": i,
                "sequence_to_order": con,
                "length_nt": len(con),
                "added_left": prefix,
                "added_left_detail": l_desc,
                "added_right": suffix,
                "added_right_detail": r_desc,
                "added_nt_total": len(con) - len(u["oligo"]),
                "variants": u["variants"],
                "degenerate_sites": sum(1 for c in u["oligo"] if c not in "ACGT"),
            })

    # ---- checks ----------------------------------------------------------- #
    checks = []

    n_bad = [n for _f, n, c, _i, _j, _ins in constructs if count_sites(c) != 2]
    checks.append(("exactly 2 Type IIS sites per construct",
                   not n_bad,
                   f"{len(constructs)} constructs, all with exactly 2"
                   if not n_bad else f"{len(n_bad)} wrong: {n_bad[:5]}"))

    bad_dig = []
    for f, n, c, i, j, ins in constructs:
        if digest(c, i, j) != ins:
            bad_dig.append(n)
    checks.append(("simulated digest releases the intended insert",
                   not bad_dig,
                   f"{len(constructs)} inserts recovered exactly"
                   if not bad_dig else f"{len(bad_dig)} wrong: {bad_dig[:5]}"))

    exp_oh = [(left_add[f] if f else args.backbone5 + args.five_pad,
               right_add[f]) for f in range(K)]
    bad_oh = []
    for f, n, c, i, j, ins in constructs:
        want_l = (args.backbone5 if f == 0 else
                  d["junctions"][f - 1][0])
        want_r = (args.backbone3 if f == K - 1 else d["junctions"][f][0])
        if ins[:4] != want_l or ins[-4:] != want_r:
            bad_oh.append(f"{n} {ins[:4]}/{ins[-4:]} want {want_l}/{want_r}")
    checks.append(("overhangs match the design's junction tokens",
                   not bad_oh,
                   "every insert starts and ends on its designed overhang"
                   if not bad_oh else "; ".join(bad_oh[:4])))

    # ---- round trip: digest, ligate, compare ----------------------------- #
    per_frag = [[x for x in constructs if x[0] == f] for f in range(K)]
    total = 1
    for pf in per_frag:
        total *= len(pf)
    combos = itertools.islice(itertools.product(*per_frag), args.max_roundtrip)
    n_rt, bad_rt = 0, []
    for combo in combos:
        pieces = [digest(c, i, j) for _f, _n, c, i, j, _ins in combo]
        prod, err = ligate(pieces)
        n_rt += 1
        if err:
            bad_rt.append(err)
            continue
        want = (args.backbone5 + args.five_pad
                + "".join(frags[f][int(n.split("_oligo")[1].split("_")[0]) - 1]
                          ["oligo"] for f, n, _c, _i, _j, _ins in combo)
                + args.backbone3)
        if prod != want:
            bad_rt.append(f"{[n for _f, n, *_ in combo]}: product != design")
    checks.append((f"round trip: digest + ligate reproduces the design",
                   not bad_rt,
                   f"{n_rt} of {total:,} assemblies rebuilt exactly"
                   if not bad_rt else "; ".join(bad_rt[:3])))

    # ---- frame and translation ------------------------------------------- #
    sample = [digest(c, i, j) for _f, _n, c, i, j, _ins in
              [per_frag[f][0] for f in range(K)]]
    prod, _ = ligate(sample)
    orf = args.backbone_pre + prod + args.backbone_post
    in_frame = len(orf) % 3 == 0
    prot = translate(orf.replace("Y", "T").replace("S", "C").replace("R", "A")
                     .replace("M", "A").replace("W", "A").replace("K", "G"))
    checks.append(("assembled ORF is in frame",
                   in_frame,
                   f"{len(orf)} nt with backbone padding "
                   f"({args.backbone_pre}|...|{args.backbone_post}) = "
                   f"{len(orf)//3} codons" if in_frame else
                   f"{len(orf)} nt is not a multiple of 3"))
    checks.append(("no internal stop codon",
                   "*" not in prot[:-1],
                   f"{len(prot)} residues, first 4 = {prot[:4]}, "
                   f"last 4 = {prot[-4:]}" if "*" not in prot[:-1]
                   else f"stop at residue {prot.index('*')+1}"))

    def codons_of(seq):
        return [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]

    added_codons = codons_of(args.backbone_pre + args.backbone5 + args.five_pad)
    tail_codons = codons_of(args.backbone3 + args.backbone_post)
    rare_hit = [c for c in added_codons + tail_codons if c in RARE]
    # ADVISORY, not a failure.  A rare codon costs some translation speed; a stop
    # costs the whole protein.  With --five-pad empty the codon boundary is forced
    # 2 nt into CGGA, which makes GGA (rare Gly) unavoidable at that junction --
    # a real consequence of the no-padding decision, worth surfacing but not worth
    # blocking an order wet lab has already signed off.
    checks.append(("junction codons translate without a stop",
                   "*" not in translate("".join(added_codons + tail_codons)),
                   f"N-term {'-'.join(added_codons)}="
                   f"{translate(''.join(added_codons))}, "
                   f"C-term {'-'.join(tail_codons)}="
                   f"{translate(''.join(tail_codons))}"
                   + (f"   NOTE rare in E. coli: {', '.join(rare_hit)}"
                      if rare_hit else "   none rare")))

    # ---- write ------------------------------------------------------------ #
    fa = os.path.join(outdir, "constructs.fasta")
    with open(fa, "w", newline="\n") as fh:
        for _f, n, c, _i, _j, ins in constructs:
            fh.write(f">{n} left={ins[:4]} right={ins[-4:]} len={len(c)}\n{c}\n")
    cs = os.path.join(outdir, "constructs.csv")
    with open(cs, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Bare two-column version for pasting straight into a vendor's upload form,
    # where any extra column is something to strip by hand.
    cs2 = os.path.join(outdir, "constructs_simple.csv")
    with open(cs2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "sequence_to_order"])
        w.writeheader()
        w.writerows({"name": r["name"],
                     "sequence_to_order": r["sequence_to_order"]} for r in rows)

    lines = []
    lines.append("ORDER-READY CONSTRUCTS")
    lines.append("=" * 70)
    lines.append(f"run          : {args.run_dir}")
    lines.append(f"fragments    : {K}   constructs: {len(constructs)}")
    lines.append("primer tails : none (fragments are not amplified)")
    lines.append(f"length range : {min(len(c) for _f,_n,c,*_ in constructs)}"
                 f"-{max(len(c) for _f,_n,c,*_ in constructs)} nt")
    lines.append("")
    lines.append("PER FRAGMENT")
    for f in range(K):
        ins0 = inserts[f][0]
        lines.append(f"  fragment {f+1}: {len(per_frag[f])} constructs, "
                     f"left {ins0[:4]}  right {ins0[-4:]}")
    lines.append("")
    lines.append("CHECKS")
    npass = 0
    for name, ok, detail in checks:
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<52} {detail}")
        npass += bool(ok)
    lines.append("")
    lines.append(f"{npass}/{len(checks)} checks pass")
    if problems:
        lines.append("")
        for p in problems:
            lines.append("  PROBLEM: " + p)
    lines.append("")
    lines.append("constructs.csv COLUMNS")
    lines.append("  name                 fragment, oligo index, and how many "
                 "proteins that oligo makes")
    lines.append("  sequence_to_order    ORDER THIS. The complete construct, "
                 "flanks included")
    lines.append("  length_nt            length of sequence_to_order")
    lines.append("  added_left/_right    exactly what was prepended and appended "
                 "to the designed oligo,")
    lines.append("                       with a _detail column naming each piece "
                 "and why it is there")
    lines.append("  added_nt_total       added_left + added_right, in nt")
    lines.append("  variants             proteins this one oligo produces "
                 "(2 if it carries a degenerate codon)")
    lines.append("  degenerate_sites     IUPAC positions in it; see "
                 "degenerate_codons.tsv for what they encode")
    lines.append("")
    lines.append("constructs_simple.csv is the same 61 sequences with only the "
                 "name and sequence_to_order")
    lines.append("  columns, for pasting into a vendor upload form.")
    lines.append("")
    lines.append("THE BACKBONE MUST SUPPLY")
    lines.append(f"  immediately BEFORE {args.backbone5}: {args.backbone_pre!r} "
                 f"(completes codon {added_codons[0]})")
    lines.append(f"  immediately AFTER  {args.backbone3}: {args.backbone_post!r} "
                 f"(completes codon {tail_codons[1]})")

    rep = "\n".join(lines)
    with open(os.path.join(outdir, "order_report.txt"), "w", newline="\n") as fh:
        fh.write(rep + "\n")
    print(rep)
    print()
    print("wrote", fa)
    print("wrote", cs)
    print("wrote", cs2)
    return 0 if npass == len(checks) and not problems else 1


if __name__ == "__main__":
    sys.exit(main())
