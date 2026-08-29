"""
assemble_cores.py -- Stage 0, LAST step:  trimmed ORFs  ->  design-ready core FASTA.

Issue #84.  Self-contained; imports nothing from this repo.

    fetch_clusters.py -> fold -> domain_from_pae.py -> [mafft] -> THIS -> unified_design.py

WHAT THIS PRODUCES (the Stage 1 input contract)
-----------------------------------------------
Two files, sharing `>core<N>_n<k>` headers, matching ninetypidorfs/cluster*.core*:

    <prefix>.core.aln.fasta   ALIGNED, uniform width, gaps kept  <- unified_design reads THIS
    <prefix>.core.fasta       ungapped, ragged lengths           <- the biological cores
    <prefix>.assembly.tsv     per-ORF audit: kept/dropped + reason + which core

`k` in `_n<k>` is w(s): how many natural ORFs collapsed onto that unique core.
This script is where w(s) comes into existence -- weighted coverage in Stage 2
is computed against these counts, so a sequence dropped here is invisible to
every downstream metric.

PHASES (the MSA is a terminal command, so it sits between them)
---------------------------------------------------------------
    trim      cluster.fa + boundaries.tsv  ->  cores.raw.fa   (ungapped, trimmed)
    [ mafft --auto cores.raw.fa > cores.aln.fa ]
    finalize  cores.aln.fa                 ->  the three files above

ORDER MATTERS IN `finalize`
---------------------------
Truncated ORFs are dropped BEFORE per-column occupancy is computed.  A fragment
left in the pool smears the occupancy profile and makes conserved columns look
like indels, which mislabels the domain block.  Occupancy is therefore recomputed
on survivors only.  A second, finer completeness check then runs on the TRIMMED
coordinates -- raw ORF length cannot distinguish "extra material" from "extra
material plus missing material" (an ORF carrying a fusion partner is long even
when its enzyme half is truncated).

A DELIBERATE OMISSION
---------------------
A truncated ORF is excluded here as both a target AND a weight.  The arguably
better treatment -- exclude it as a target but let its weight flow to the
full-length core it matches over its covered span -- is NOT implemented, because
it is a scoping decision, not an oversight.  With ~85/213 fragments in cluster2
it is not a rounding error; `assembly.tsv` records every drop so the shift can be
measured before deciding.

USAGE
-----
    python assemble_cores.py trim --cluster-fa out/clusters/13/cluster.fa \\
        --boundaries out/clusters/13/boundaries.tsv -o out/clusters/13/cores.raw.fa
    mafft --auto out/clusters/13/cores.raw.fa > out/clusters/13/cores.aln.fa
    python assemble_cores.py finalize --aln out/clusters/13/cores.aln.fa \\
        --orf-meta out/orf_meta.tsv --prefix out/clusters/13/cluster13
"""

import argparse
import os
import sys
from collections import Counter, OrderedDict

GAP = "-"


# --------------------------------------------------------------------------- #
# io
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """[(header, seq)] preserving file order."""
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    if not out:
        sys.exit("no sequences read from %s" % path)
    return out


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(">%s\n" % h)
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def orf_id_of(header):
    """Field 1 of the PETadex ORF header: >{orf_id}|{acc}|{run}|..."""
    return header.split("|", 1)[0].strip()


# --------------------------------------------------------------------------- #
# phase: trim
# --------------------------------------------------------------------------- #

def read_boundaries(path):
    """TSV: orf_id <tab> start <tab> end,  1-based INCLUSIVE amino-acid coords.
    This is the contract domain_from_pae.py must emit."""
    bounds = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n\r").split("\t")
        if head and head[0].lower() not in ("orf_id", "orfid"):
            fh.seek(0)                      # no header row
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 3 or not (f[1].isdigit() and f[2].isdigit()):
                continue
            bounds[f[0].strip()] = (int(f[1]), int(f[2]))
    return bounds


def phase_trim(args):
    records = read_fasta(args.cluster_fa)
    bounds = read_boundaries(args.boundaries) if args.boundaries else {}
    if not bounds:
        print("WARNING: no boundaries given -- passing sequences through UNTRIMMED. "
              "This is a pipeline smoke-test path, not a real run.", file=sys.stderr)

    out, n_trimmed, missing = [], 0, []
    for h, s in records:
        oid = orf_id_of(h)
        if oid in bounds:
            a, b = bounds[oid]
            a = max(1, a)
            b = min(len(s), b)
            if a > b:
                missing.append(oid)
                continue
            s = s[a - 1:b]
            n_trimmed += 1
        elif bounds:
            missing.append(oid)
            continue
        out.append((h, s.replace(GAP, "")))

    write_fasta(args.out, out)
    print("%d sequences in, %d trimmed, %d written -> %s"
          % (len(records), n_trimmed, len(out), args.out))
    if missing:
        print("  %d ORF(s) dropped: no usable boundary (first few: %s)"
              % (len(missing), missing[:5]), file=sys.stderr)


# --------------------------------------------------------------------------- #
# phase: finalize
# --------------------------------------------------------------------------- #

def longest_block(occ, threshold):
    """Longest contiguous run of columns with occupancy >= threshold."""
    best = cur = None
    for j, o in enumerate(occ):
        if o >= threshold:
            cur = (cur[0], j + 1) if cur else (j, j + 1)
            if best is None or (cur[1] - cur[0]) > (best[1] - best[0]):
                best = cur
        else:
            cur = None
    return best


def read_orf_meta(path):
    meta = {}
    with open(path) as fh:
        cols = fh.readline().rstrip("\n\r").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        for line in fh:
            f = line.rstrip("\n\r").split("\t")
            if len(f) < len(cols):
                continue
            meta[f[idx["orf_id"]]] = f
    return meta, idx


def phase_finalize(args):
    records = read_fasta(args.aln)
    widths = {len(s) for _, s in records}
    if len(widths) != 1:
        sys.exit("input is not aligned: %d distinct lengths %s -- run mafft first"
                 % (len(widths), sorted(widths)[:5]))
    L = widths.pop()

    meta, midx = ({}, {})
    if args.orf_meta:
        meta, midx = read_orf_meta(args.orf_meta)

    audit = OrderedDict()          # orf_id -> [status, reason, core]
    for h, _ in records:
        audit[orf_id_of(h)] = ["kept", "", ""]

    # -- 1. drop truncations FIRST, so they cannot smear the occupancy profile --
    ungapped = [(h, s, len(s.replace(GAP, ""))) for h, s in records]
    lens = sorted(n for _, _, n in ungapped)
    med = lens[len(lens) // 2]
    min_len = args.min_len_frac * med

    survivors = []
    for h, s, n in ungapped:
        oid = orf_id_of(h)
        if n < min_len:
            audit[oid] = ["dropped", "short: %d aa < %.0f (%.2f x median %d)"
                          % (n, min_len, args.min_len_frac, med), ""]
            continue
        if args.drop_contig_edge and meta.get(oid) and "contig_edge_5p" in midx:
            if meta[oid][midx["contig_edge_5p"]] == "1":
                audit[oid] = ["dropped", "contig_edge_5p (start_nt == 0)", ""]
                continue
        survivors.append((h, s))

    if not survivors:
        sys.exit("every sequence was dropped -- loosen --min-len-frac")

    # -- 2. occupancy on survivors only, then the domain block --
    occ = [sum(1 for _, s in survivors if s[j] != GAP) / len(survivors)
           for j in range(L)]
    block = longest_block(occ, args.min_occupancy)
    if block is None:
        sys.exit("no column reaches occupancy %.2f -- loosen --min-occupancy"
                 % args.min_occupancy)
    a, b = block

    # -- 3. completeness check on TRIMMED coordinates --
    kept = []
    for h, s in survivors:
        piece = s[a:b]
        cov = sum(1 for c in piece if c != GAP) / (b - a)
        if cov < args.min_block_coverage:
            audit[orf_id_of(h)] = ["dropped", "block coverage %.2f < %.2f"
                                   % (cov, args.min_block_coverage), ""]
            continue
        kept.append((h, piece))

    if not kept:
        sys.exit("every sequence failed the block-coverage check -- "
                 "loosen --min-block-coverage")

    # -- 4. dedup on the UNGAPPED slice; that is what "the same protein" means --
    groups = OrderedDict()
    for h, piece in kept:
        key = piece.replace(GAP, "")
        groups.setdefault(key, []).append((h, piece))

    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    aln_out, raw_out = [], []
    for i, (key, members) in enumerate(ordered, 1):
        name = "core%d_n%d" % (i, len(members))
        aln_out.append((name, members[0][1]))
        raw_out.append((name, key))
        for h, _ in members:
            audit[orf_id_of(h)] = ["kept", "", name]

    prefix = args.prefix
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
    write_fasta(prefix + ".core.aln.fasta", aln_out)
    write_fasta(prefix + ".core.fasta", raw_out)

    with open(prefix + ".assembly.tsv", "w") as fh:
        fh.write("orf_id\tstatus\treason\tcore\n")
        for oid, row in audit.items():
            fh.write("%s\t%s\t%s\t%s\n" % (oid, row[0], row[1], row[2]))

    n_drop = sum(1 for r in audit.values() if r[0] == "dropped")
    total_w = sum(len(m) for _, m in ordered)
    print("input          %d aligned sequences, width %d" % (len(records), L))
    print("dropped        %d (%d short/edge, %d low block coverage)"
          % (n_drop,
             sum(1 for r in audit.values() if r[0] == "dropped"
                 and not r[1].startswith("block coverage")),
             sum(1 for r in audit.values() if r[1].startswith("block coverage"))))
    print("domain block   columns [%d,%d)  width %d  (occupancy >= %.2f)"
          % (a, b, b - a, args.min_occupancy))
    print("unique cores   %d   total weight %d   heaviest n=%d"
          % (len(ordered), total_w, len(ordered[0][1])))
    print("\nwrote %s.core.aln.fasta   <- Stage 1 input" % prefix)
    print("wrote %s.core.fasta" % prefix)
    print("wrote %s.assembly.tsv" % prefix)


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="phase", required=True)

    s = sub.add_parser("trim", help="apply per-ORF domain boundaries")
    s.add_argument("--cluster-fa", required=True)
    s.add_argument("--boundaries", help="TSV: orf_id, start, end (1-based inclusive). "
                                        "Omit only for a smoke test.")
    s.add_argument("-o", "--out", required=True)
    s.set_defaults(func=phase_trim)

    s = sub.add_parser("finalize", help="aligned FASTA -> deduped core files")
    s.add_argument("--aln", required=True, help="output of mafft")
    s.add_argument("--orf-meta", help="orf_meta.tsv from fetch_clusters.py")
    s.add_argument("--prefix", required=True, help="output path prefix")
    s.add_argument("--min-len-frac", type=float, default=0.8,
                   help="drop sequences shorter than this fraction of the median (0.8)")
    s.add_argument("--min-occupancy", type=float, default=0.5,
                   help="column occupancy defining the domain block (0.5)")
    s.add_argument("--min-block-coverage", type=float, default=0.9,
                   help="min non-gap fraction within the block (0.9)")
    s.add_argument("--drop-contig-edge", action="store_true",
                   help="also drop ORFs with start_nt == 0. OFF by default: the "
                        "signal is real but catches only the 5' edge and is not a "
                        "clean separator.")
    s.set_defaults(func=phase_finalize)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
