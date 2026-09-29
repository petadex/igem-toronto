#!/usr/bin/env python3
"""sto_to_cores.py -- Stockholm alignment (from `hmmsearch -A`) -> core FASTAs.

A dependency-free stand-in for `esl-reformat` when the Easel miniapps are not
installed.  Issue algocontrols (IsPETase / LCC positive controls).

    python3 sto_to_cores.py <in.sto> <out.aln.fa> <out.cores.fa>

out.aln.fa    cores aligned to the HMM match states (gaps kept, uniform width,
              uppercased, '.' insert-gaps normalised to '-')
out.cores.fa  ungapped core sequences (all gap chars stripped)

hmmsearch -A writes an interleaved Stockholm file: each sequence's alignment is
split across several blocks separated by blank lines, so the residues for one
name must be concatenated across the whole file (not just the first block).
Sequence names carry a `/start-end` suffix marking where the core sat in the
original ORF; that is preserved verbatim.
"""
import sys


def read_stockholm(path):
    """Return [(name, aligned_seq), ...] in first-seen order, concatenating the
    interleaved blocks.  Skips markup (#=GC/#=GS/#=GR, comments) and the // end."""
    seqs, order = {}, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line == "//":
                continue
            parts = line.split()
            if len(parts) != 2:           # a real seq line is exactly: name residues
                continue
            name, block = parts
            if name not in seqs:
                seqs[name] = []
                order.append(name)
            seqs[name].append(block)
    return [(n, "".join(seqs[n])) for n in order]


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    src, aln_out, core_out = sys.argv[1:]
    records = read_stockholm(src)
    if not records:
        sys.exit("no sequences parsed from %s -- did hmmsearch find any hits?" % src)

    with open(aln_out, "w") as fa:
        for name, seq in records:
            fa.write(">%s\n%s\n" % (name, seq.upper().replace(".", "-")))
    with open(core_out, "w") as fc:
        for name, seq in records:
            ungapped = seq.upper().replace("-", "").replace(".", "")
            fc.write(">%s\n%s\n" % (name, ungapped))

    width = len(records[0][1])
    print("%d cores | aln width %d -> %s , %s"
          % (len(records), width, aln_out, core_out))


if __name__ == "__main__":
    main()
