"""Make degenerate codons OBVIOUS in an order sheet.

cutsearch_design.py writes IUPAC codes inline in the fragment FASTAs, which is
correct but easy to miss: one letter in a 500-nt string looks like a sequencing
ambiguity rather than a deliberate design choice.  This reads a finished run
directory and emits, beside it:

  degenerate_codons.tsv          one row per degenerate site, with expansions
  fragment<N>.annotated.fasta    same sequences, headers naming every site

Reads only; the design is never recomputed and the original files are untouched.

  python annotate_degenerate.py <run-dir>
"""
import os
import sys
import glob
from _paths import HERE                      # noqa: F401  (puts HERE on sys.path)
import cutsearch_design as U


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def sites(seq):
    """Every degenerate position, as a dict describing the codon it falls in."""
    out = []
    for i, ch in enumerate(seq):
        if ch in "ACGT":
            continue
        ci = i // 3                              # 0-based codon index
        codon = seq[ci * 3:ci * 3 + 3]
        exps = sorted({a + b + c
                       for a in U.IUPAC[codon[0]]
                       for b in U.IUPAC[codon[1]]
                       for c in U.IUPAC[codon[2]]})
        aas = sorted({U.AA_BY_CODON[e] for e in exps})
        out.append({"nt_pos": i + 1, "codon_index": ci + 1, "iupac_base": ch,
                    "base_expands": "/".join(U.IUPAC[ch]), "codon": codon,
                    "codon_expands": "/".join(exps), "amino_acids": "/".join(aas)})
    return out


def main(run_dir):
    frags = sorted(glob.glob(os.path.join(run_dir, "fragment*.fasta")))
    frags = [f for f in frags if ".annotated." not in f]
    if not frags:
        sys.exit("no fragment*.fasta found in " + run_dir)

    rows, n_deg, n_oligos = [], 0, 0
    for path in frags:
        frag = os.path.basename(path).replace(".fasta", "")
        ann = []
        for name, seq in read_fasta(path):
            n_oligos += 1
            st = sites(seq)
            n_deg += len(st)
            if st:
                tags = ["nt{nt_pos} {iupac_base}({base_expands}) "
                        "codon{codon_index} {codon}>{amino_acids}".format(**s)
                        for s in st]
                hdr = (name + "  |  *** " + str(len(st)) +
                       " DEGENERATE SITE(S) ***  " + "; ".join(tags))
            else:
                hdr = name + "  |  no degenerate sites"
            ann.append((hdr, seq))
            for s in st:
                rows.append([frag, name, s["nt_pos"], s["codon_index"],
                             s["iupac_base"], s["base_expands"], s["codon"],
                             s["codon_expands"], s["amino_acids"]])
        out = path.replace(".fasta", ".annotated.fasta")
        with open(out, "w") as fh:
            for hdr, seq in ann:
                fh.write(">" + hdr + chr(10) + seq + chr(10))

    tsv = os.path.join(run_dir, "degenerate_codons.tsv")
    cols = ["fragment", "oligo", "nt_pos", "codon_index", "iupac_base",
            "base_expands", "codon", "codon_expands", "amino_acids"]
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + chr(10))
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + chr(10))

    print("run       : " + run_dir)
    print("oligos    : " + str(n_oligos))
    print("degenerate: " + str(n_deg) + " site(s) in " +
          str(len({r[1] for r in rows})) + " oligo(s)")
    print("wrote     : degenerate_codons.tsv + " + str(len(frags)) +
          " annotated FASTA(s)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
