"""Trim an aligned FASTA to its conserved domain (the structured enzyme core).

Pipeline:
  1. drop fragments by ungapped length (BEFORE occupancy, so they don't smear it)
  2. recompute per-column occupancy on full-length seqs only
  3. domain = longest contiguous block of columns with occ >= OCC_T
  4. trim each seq to that block, ungap, dedup (keep counts)

Pure stdlib -- run with Windows `py`.

Usage:
    py trim_to_core.py cluster2.aln.fasta cluster2.core.fasta
"""
import sys
from collections import Counter

infile  = sys.argv[1]
outfile = sys.argv[2]
FRAG_FRAC = 0.75   # seq is a fragment if ungapped len < FRAG_FRAC * median
OCC_T     = 0.50   # column is "domain" if >= this frac of full-len seqs non-gap
CONS_T    = 0.90   # column is "conserved" (vs variable loop) if occ >= this

def read_fasta(path):
    names, seqs, name, buf = [], [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                seqs.append("".join(buf))
            name, buf = line[1:], []
            names.append(name)
        elif line:
            buf.append(line)
    if name is not None:
        seqs.append("".join(buf))
    return names, seqs

names, aln = read_fasta(infile)
N, L = len(aln), len(aln[0])

# --- 1. drop fragments by ungapped length ---
ungapped = [sum(1 for ch in s if ch != "-") for s in aln]
med = sorted(ungapped)[N // 2]
full = [(nm, s) for nm, s, u in zip(names, aln, ungapped) if u >= FRAG_FRAC * med]
frags = [nm for nm, s, u in zip(names, aln, ungapped) if u < FRAG_FRAC * med]
F = len(full)
print(f"alignment: {N} seqs x {L} cols | median ungapped len = {med}")
print(f"fragments dropped (len < {FRAG_FRAC:.0%} median = {FRAG_FRAC*med:.0f}): {len(frags)}")
print(f"full-length seqs kept for profiling: {F}")

# --- 2. occupancy on full-length only ---
occ = [sum(1 for _, s in full if s[c] != "-") / F for c in range(L)]

# coarse profile (each char = a 1/25-of-length bucket, max occ in bucket)
buckets = 50
prof = []
for b in range(buckets):
    lo, hi = b * L // buckets, (b + 1) * L // buckets
    m = max(occ[lo:hi]) if hi > lo else 0
    prof.append(" .:-=+*#@"[min(8, int(m * 8.999))])
print(f"occupancy profile (cols 0..{L-1}, ' '=empty '@'=full):\n  {''.join(prof)}")

# --- 3. longest contiguous domain block ---
best, start = (0, 0, 0), None
for c in range(L + 1):
    inrun = c < L and occ[c] >= OCC_T
    if inrun and start is None:
        start = c
    elif not inrun and start is not None:
        if c - start > best[0]:
            best = (c - start, start, c)
        start = None
blk_len, cs, ce = best
cols = list(range(cs, ce))
var_cols = [c for c in cols if occ[c] < CONS_T]
print(f"domain block: cols {cs}..{ce-1} ({blk_len} cols) | scaffolding cols removed: {L-blk_len}")
print(f"  conserved cols (occ>={CONS_T}): {blk_len-len(var_cols)} | variable/loop cols: {len(var_cols)}")

# --- 4. trim, dedup ---
kept = [(nm, "".join(s[c] for c in cols).replace("-", "")) for nm, s in full]
uniq = Counter(seq for _, seq in kept)
with open(outfile, "w") as out:
    for i, (seq, cnt) in enumerate(uniq.most_common(), 1):
        out.write(f">core{i}_n{cnt}\n{seq}\n")

lens = Counter(len(seq) for _, seq in kept)
print(f"trimmed-core length histogram: {dict(sorted(lens.items()))}")
print(f"unique cores after dedup: {len(uniq)} (wrote -> {outfile})")
