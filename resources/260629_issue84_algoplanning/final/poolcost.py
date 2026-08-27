"""place_cuts dominates at high K.  How much does enlarging the proxy pool cost?
(Random candidates need no DFS at all, so this prices the proxy half only.)

  python poolcost.py                       # whichever default inputs are present
  python poolcost.py cores.aln.fasta ...   # explicit
"""
import sys, time
from _paths import inputs_from_argv
import cutsearch_design as U

for label, path in inputs_from_argv(sys.argv):
    seqs = U.read_aligned_cores(path)
    al = [s for s, _ in seqs]; L = len(al[0]); const = U.constant_columns(al)
    print(f"\n{label}: n={len(al)}, {L} cols")
    print(f"   {'K':>3} {'n_keep':>7} {'place_cuts':>11} {'pool':>6} {'trunc':>6}")
    for K in (3, 5, 6):
        for nk in (1, 5, 50):
            t = time.time()
            c, tr = U.place_cuts(al, L, K, 20, const, "gg", 6,
                                 U.BACKBONE_OVERHANGS, n_keep=nk)
            dt = time.time() - t
            print(f"   {K:>3} {nk:>7} {dt:>10.2f}s {len(c):>6} {str(tr):>6}")
