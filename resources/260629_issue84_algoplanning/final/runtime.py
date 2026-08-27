"""How does one greedy scale with cluster size?  Times real inputs plus
subsamples, so 50 candidates per K can be costed for any cluster.

  python runtime.py                       # whichever default inputs are present
  python runtime.py cores.aln.fasta ...   # explicit
"""
import sys, time, random
from _paths import inputs_from_argv, CLUSTER1
import cutsearch_design as U


def timeone(aligned, weights, K, cap=300, mb=20):
    L = len(aligned[0]); const = U.constant_columns(aligned)
    t0 = time.time()
    c, _ = U.place_cuts(aligned, L, K, mb, const, "gg", 6,
                        U.BACKBONE_OVERHANGS, n_keep=1)
    tp = time.time() - t0
    if not c: return None
    t1 = time.time()
    d = U.greedy(aligned, weights, c[0][1], c[0][2], L, "gg", 1.0, 3,
                 max_library=cap)
    tg = time.time() - t1
    return tp, tg, (d["W"] if d else None), L


pairs = inputs_from_argv(sys.argv)
print(f"{'input':>26} {'n':>5} {'cols':>6} {'K':>3} {'place':>7} {'greedy':>8} {'W':>5}")
for label, path in pairs:
    seqs = U.read_aligned_cores(path)
    al = [s for s, _ in seqs]; we = [w for _, w in seqs]
    for K in (3, 5):
        r = timeone(al, we, K)
        if r: print(f"{label:>26} {len(al):>5} {r[3]:>6} {K:>3} "
                    f"{r[0]:>6.2f}s {r[1]:>7.2f}s {r[2]:>5}")

# n-scaling curve: subsample the LARGEST input we were given
big = max(pairs, key=lambda lp: len(U.read_aligned_cores(lp[1])))
seqs = U.read_aligned_cores(big[1])
print(f"\nsubsampled {big[0]} (K=3), one greedy:")
step = max(1, len(seqs) // 5)
for n in list(range(step, len(seqs), step)) + [len(seqs)]:
    rg = random.Random(3)
    idx = sorted(rg.sample(range(len(seqs)), min(n, len(seqs))))
    al = [seqs[i][0] for i in idx]; we = [seqs[i][1] for i in idx]
    r = timeone(al, we, 3)
    if r: print(f"   n={n:>3}: greedy {r[1]:>6.2f}s   (place {r[0]:.2f}s)  W={r[2]}")
