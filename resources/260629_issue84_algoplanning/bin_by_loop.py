"""Find the big loop indel in an aligned cluster and split into length bins.

  1. drop fragments by ungapped length
  2. domain block (longest occ>=OCC_T run); variable cols = occ<CONS_T inside it
  3. find the largest contiguous variable run = the loop
  4. bin each seq by loop presence (non-gap fraction over loop cols)
  5. within each bin, re-trim to that bin's occupied cols; report length constancy

Usage:  py bin_by_loop.py cluster2.aln.fasta
"""
import sys
from collections import Counter

infile = sys.argv[1]
FRAG_FRAC, OCC_T, CONS_T = 0.75, 0.50, 0.90
MERGE = 3       # merge variable-col runs separated by <= this many conserved cols

def read_fasta(path):
    names, seqs, name, buf = [], [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None: seqs.append("".join(buf))
            name, buf = line[1:], []; names.append(name)
        elif line: buf.append(line)
    if name is not None: seqs.append("".join(buf))
    return names, seqs

names, aln = read_fasta(infile)
N, L = len(aln), len(aln[0])
ungapped = [sum(ch != "-" for ch in s) for s in aln]
med = sorted(ungapped)[N // 2]
full = [(nm, s) for nm, s, u in zip(names, aln, ungapped) if u >= FRAG_FRAC * med]
F = len(full)

occ = [sum(s[c] != "-" for _, s in full) / F for c in range(L)]
# domain block
best, start = (0, 0, 0), None
for c in range(L + 1):
    inrun = c < L and occ[c] >= OCC_T
    if inrun and start is None: start = c
    elif not inrun and start is not None:
        if c - start > best[0]: best = (c - start, start, c)
        start = None
_, cs, ce = best

# variable columns inside domain, merged into runs
var = [c for c in range(cs, ce) if occ[c] < CONS_T]
runs, cur = [], []
for c in var:
    if cur and c - cur[-1] - 1 > MERGE:
        runs.append((cur[0], cur[-1])); cur = []
    cur.append(c)
if cur: runs.append((cur[0], cur[-1]))
runs.sort(key=lambda r: -(r[1] - r[0]))
print(f"full-length seqs: {F} | domain cols {cs}..{ce-1} | {len(var)} variable cols")
print(f"variable runs (start,end,width): {[(a,b,b-a+1) for a,b in runs[:6]]}")

# significant indel sites = variable runs at least this wide
sites = [(a, b) for a, b in runs if b - a + 1 >= 8]
sites.sort()       # left-to-right
print(f"significant indel sites (>=8 wide): {[(a,b,b-a+1) for a,b in sites]}")

def present(s, a, b):
    return sum(s[c] != "-" for c in range(a, b + 1)) / (b - a + 1) >= 0.5

# classify each seq by its presence pattern across the sites
groups = {}
for nm, s in full:
    pat = tuple(present(s, a, b) for a, b in sites)
    groups.setdefault(pat, []).append((nm, s))

def core_len(s, group):
    g = len(group)
    bocc = [sum(t[c] != "-" for _, t in group) / g for c in range(cs, ce)]
    keepcols = [cs + i for i, o in enumerate(bocc) if o >= 0.5]
    return len("".join(s[c] for c in keepcols).replace("-", ""))

print(f"\nlength classes by indel pattern (site order = {[a for a,_ in sites]}):")
for pat, group in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    lens = Counter(core_len(s, group) for _, s in group)
    uniq = len({"".join(s[c] for c in range(cs, ce)).replace("-", "") for _, s in group})
    lo, hi = min(lens), max(lens)
    tag = "+".join("N" if p and i == 0 else "C" if p and i == 1 else "-"
                   for i, p in enumerate(pat))
    print(f"  pattern {pat} [{tag}]: {len(group):3d} seqs | "
          f"len {lo}..{hi} (spread {hi-lo:2d}) | unique {uniq} | {dict(sorted(lens.items()))}")
