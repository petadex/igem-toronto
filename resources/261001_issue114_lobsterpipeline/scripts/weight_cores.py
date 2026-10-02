#!/usr/bin/env python3
"""weight_cores.py -- HMM-extracted cores  ->  design-ready weighted cores.

Issue algocontrols (IsPETase / LCC positive controls).  Self-contained.

    hmmsearch -A -> sto_to_cores.py -> THIS -> cutsearch_design.py

This is Stage-0's `finalize` step for the HMM route: it does what the old
assemble_cores.py did (drop truncated cores, then collapse identical cores into
weights) but starts from an alignment that hmmsearch already made against the
C001 model, so there is no separate mafft pass.

    python3 weight_cores.py <in.aln.fa> <out_prefix> [--min-cov 0.80]

INPUT   the aligned cores from sto_to_cores.py: one row per ORF, uniform width,
        headers `>{orf_id}|...|/{start-end}`.

OUTPUT  <out_prefix>.weighted.aln.fa   aligned, `>core<N>_n<k>`  <- cutsearch reads THIS
        <out_prefix>.weighted.fa       ungapped representative of each unique core
        <out_prefix>.weight_audit.tsv  per-ORF: kept/dropped + which core it fell on

WHY DROP TRUNCATED CORES FIRST (order matters, same as assemble_cores.py)
------------------------------------------------------------------------
A partial core is mostly gaps across the alignment; left in the pool it smears
the per-column occupancy the design relies on to place junctions, and it inflates
a core's weight with a sequence that does not actually carry the whole domain.
So completeness is judged BEFORE weights are counted: a dropped ORF is invisible
as both a target and a weight.  Completeness here = fraction of alignment columns
this row actually occupies (non-gap / width); the C001 model is ~221 match states,
so a genuine core sits near 1.0 and a fragment sits well below --min-cov.

WHY IDENTICAL CORES COLLAPSE
----------------------------
w(s) = how many natural ORFs share one unique core.  The clusters are dominated
by re-sequenced lab constructs (His-tagged, fused IsPETase/LCC variants), so many
ORFs reduce to the same catalytic domain once the flanks are gone.  Each distinct
aligned core becomes one `>core<N>_n<k>` record with k = that count; the design
then weights coverage by k.  Cores are numbered heaviest-first, matching the
greedy's "seed from the heaviest real core".
"""
import argparse
import sys

GAP = set("-.")


def read_aln(path):
    """[(header, aligned_seq)] -- uppercased, '.' normalised to '-'."""
    recs, header, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(buf)))
                header, buf = line[1:], []
            else:
                buf.append(line.upper().replace(".", "-"))
    if header is not None:
        recs.append((header, "".join(buf)))
    if not recs:
        sys.exit(f"no sequences read from {path}")
    widths = {len(s) for _, s in recs}
    if len(widths) != 1:
        sys.exit(f"input is not aligned: {len(widths)} distinct widths {sorted(widths)}")
    return recs, widths.pop()


def orf_id_of(header):
    return header.split("|", 1)[0]


def coord_of(header):
    return header.rsplit("/", 1)[1] if "/" in header else ""


def trim_all_gap_columns(rows):
    """Drop columns that are gaps in EVERY surviving row (they only existed to
    host a fragment's insert).  Same column set removed from all rows, so width
    stays uniform.  Cannot merge distinct cores: a column that separates two rows
    has a non-gap in at least one, so it is never all-gap."""
    if not rows:
        return rows
    width = len(rows[0])
    keep = [c for c in range(width)
            if any(row[c] not in GAP for row in rows)]
    return ["".join(row[c] for c in keep) for row in rows]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_aln")
    ap.add_argument("out_prefix")
    ap.add_argument("--min-cov", type=float, default=0.80,
                    help="keep a core if non-gap/width >= this (default 0.80)")
    args = ap.parse_args()

    recs, width = read_aln(args.in_aln)

    # phase 1 -- completeness filter (before any weighting)
    kept, audit = [], []      # kept: (header, aligned_seq)
    for header, seq in recs:
        occ = sum(1 for ch in seq if ch not in GAP)
        cov = occ / width
        status = "kept" if cov >= args.min_cov else "dropped_truncated"
        audit.append({"orf_id": orf_id_of(header), "header": header,
                      "coord": coord_of(header), "occ": occ,
                      "cov": round(cov, 3), "status": status, "core": ""})
        if status == "kept":
            kept.append((header, seq, audit[-1]))

    if not kept:
        sys.exit("every core failed the coverage filter -- lower --min-cov?")

    # phase 2 -- collapse identical aligned cores into weights (heaviest first)
    groups = {}               # aligned_seq -> [audit rows]
    for _h, seq, arow in kept:
        groups.setdefault(seq, []).append(arow)
    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))

    reps = [seq for seq, _ in ordered]
    reps = trim_all_gap_columns(reps)      # tighten to survivors' occupied columns

    aln_path = args.out_prefix + ".weighted.aln.fa"
    fa_path = args.out_prefix + ".weighted.fa"
    tsv_path = args.out_prefix + ".weight_audit.tsv"

    with open(aln_path, "w") as fa_aln, open(fa_path, "w") as fa_un:
        for i, ((_seq, members), rep) in enumerate(zip(ordered, reps), start=1):
            k = len(members)
            cid = f"core{i}"
            for arow in members:
                arow["core"] = cid
            fa_aln.write(f">{cid}_n{k}\n{rep}\n")
            fa_un.write(f">{cid}_n{k}\n{rep.replace('-', '')}\n")

    with open(tsv_path, "w") as tsv:
        tsv.write("orf_id\tcoord\tocc\tcov\tstatus\tcore\theader\n")
        for a in audit:
            tsv.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\n"
                      % (a["orf_id"], a["coord"], a["occ"], a["cov"],
                         a["status"], a["core"], a["header"]))

    n_in = len(recs)
    n_drop = sum(1 for a in audit if a["status"] != "kept")
    n_kept = len(kept)
    weights = sorted((len(m) for _s, m in ordered), reverse=True)
    print(f"{args.in_aln}")
    print(f"  in {n_in} | dropped truncated {n_drop} | kept {n_kept}"
          f" -> {len(ordered)} unique cores, total weight {sum(weights)}")
    print(f"  aligned width {width} -> {len(reps[0])} after trimming all-gap columns")
    print(f"  weights (heaviest first): {weights}")
    print(f"  wrote {aln_path}\n         {fa_path}\n         {tsv_path}")


if __name__ == "__main__":
    main()
