"""How many candidate segmentations per K do we actually need?  Evaluate 200 of
each arm, then read off best-of-first-N to find where it saturates."""
import sys, math, random, statistics as st
from _paths import inputs_from_argv
import cutsearch_design as U

ALN = inputs_from_argv(sys.argv)[0][1]
LIBCAP, POOL = 300, 120
seqs = U.read_aligned_cores(ALN)
aligned = [s for s, _ in seqs]; weights = [w for _, w in seqs]
L = len(aligned[0]); const = U.constant_columns(aligned)
RES, MB = U.BACKBONE_OVERHANGS, 20
SITES = [(p, t) for p in range(MB, L - MB + 1)
         for t in [U.junction_options(p, const, "gg", 6, RES)] if t]

def rand_segs(K, n, rng):
    out, seen, tries = [], set(), 0
    while len(out) < n and tries < 100000:
        tries += 1
        pick = sorted(rng.sample(range(len(SITES)), K - 1))
        ps = tuple(SITES[i][0] for i in pick)
        if ps in seen: continue
        b = [0] + list(ps) + [L]
        if any(b[i+1] - b[i] < MB for i in range(len(b) - 1)): continue
        toks, ok = [], True
        for i in pick:
            cand = [t for t in SITES[i][1]
                    if not any(U.tokens_conflict(t, u, "gg") for u in toks)]
            if not cand: ok = False; break
            toks.append(rng.choice(cand))
        if not ok: continue
        seen.add(ps); out.append((0.0, list(ps), toks))
    return out

def run(c):
    d = U.greedy(aligned, weights, c[1], c[2], L, "gg", 1.0, 3, max_library=LIBCAP)
    if d is None: return None
    return d["W"], sum(u.nt for units in d["layers"] for u in units)

print(f"library cap {LIBCAP}, {len(SITES)} legal sites, "
      f"{sum(weights)} natural sequences total")
for K in (5,):
    pool, _t = U.place_cuts(aligned, L, K, MB, const, "gg", 6, RES,
                            n_keep=POOL, pool_factor=1)
    prox = [r for r in (run(c) for c in pool[:POOL]) if r]
    rnd = [r for r in (run(c) for c in rand_segs(K, POOL, random.Random(7))) if r]
    print(f"\nK={K}: {len(prox)} proxy candidates, {len(rnd)} random")
    print(f"   distinct outcomes: proxy {len(set(prox))}, random {len(set(rnd))}")
    print(f"   {'N':>5} {'proxy best-of-N':>16} {'random mean':>12} {'random p10':>11}")
    rg = random.Random(11)
    for N in (1, 2, 5, 10, 20, 50, 100, 200):
        if N > len(prox): break
        pb = max(w for w, _ in prox[:N])
        sims = []
        for _ in range(400):
            sims.append(max(w for w, _ in rg.sample(rnd, min(N, len(rnd)))))
        sims.sort()
        print(f"   {N:>5} {pb:>16} {st.mean(sims):>12.1f} "
              f"{sims[int(0.10*len(sims))]:>11}")
