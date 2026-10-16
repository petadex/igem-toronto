"""Re-split the remapped PETai labels on 30% clusters, stratified.

Issue 117 step 2. Runs on the step-1 output, never on the v1 release.

The unit of assignment is the cluster_30, which keeps the split homology-clean
at 30% identity by construction. The difficulty is that the negatives are not
spread across clusters: all 200 sit in 13 of 196 clusters, and 160 of those in
the four largest clusters in the set. Any fold that gets a usable inactive
count also gets a large block of sequences, so inactive count in test trades
directly against fold size.

The search makes that trade explicit. The 13 negative-bearing clusters are
enumerated exhaustively over the three folds (3^13). The remaining 183 are
filled greedily, largest first, into whichever fold is furthest below its
sequence target. Candidates are scored and the frontier is reported so the
choice of floor is a decision rather than a default.

Usage: python3 resplit.py <remapped_dir> [out_dir]
"""
import csv, json, sys, os, collections, itertools

FOLDS = ("train", "val", "test")
TARGET = {"train": 0.70, "val": 0.15, "test": 0.15}

def load(src):
    labels = list(csv.DictReader(open(os.path.join(src, "labels.csv"), newline="")))
    seqs = list(csv.DictReader(open(os.path.join(src, "sequences.csv"), newline="")))
    by_seq = collections.defaultdict(list)
    for l in labels:
        by_seq[l["seq_md5"]].append(l)

    clusters = collections.defaultdict(lambda: collections.Counter())
    members = collections.defaultdict(list)
    for s in seqs:
        c = s["cluster_30"]
        L = by_seq[s["seq_md5"]]
        neg = sum(1 for l in L if l["is_negative"] == "True")
        halo = sum(1 for l in L if l["source"] == "halo")
        llm = sum(1 for l in L if l["source"] == "llm")
        a = clusters[c]
        a["nseq"] += 1
        a["nlab"] += len(L)
        a["neg"] += neg
        a["negseq"] += neg > 0
        a["halo"] += halo
        a["haloseq"] += halo > 0
        a["llm"] += llm
        a["human"] += len(L) - llm
        members[c].append(s["seq_md5"])
    return labels, seqs, dict(clusters), dict(members)

def totals(clusters):
    t = collections.Counter()
    for a in clusters.values():
        t.update(a)
    return t

def evaluate(assign, clusters, T):
    agg = {f: collections.Counter() for f in FOLDS}
    for c, f in assign.items():
        agg[f].update(clusters[c])
    size_pen = sum(abs(agg[f]["nseq"] / T["nseq"] - TARGET[f]) for f in FOLDS)
    lab_pen = sum(abs(agg[f]["nlab"] / T["nlab"] - TARGET[f]) for f in FOLDS)
    # Provenance: the v1 failure was LLM-heavy evaluation folds. Push LLM into
    # train and keep human-curated labels in val and test.
    llm_test = agg["test"]["llm"] / max(1, agg["test"]["nlab"])
    llm_val = agg["val"]["llm"] / max(1, agg["val"]["nlab"])
    halo_pen = sum(abs(agg[f]["haloseq"] / T["haloseq"] - TARGET[f]) for f in FOLDS)
    score = 2.0 * size_pen + 1.0 * lab_pen + 1.5 * (llm_test + llm_val) + 0.5 * halo_pen
    return score, agg


def search_all(clusters, members, floors, floor_val=10, size_tol=0.5, keep=4000):
    """One pass over the negative-cluster enumeration, best candidate per floor.

    The greedy fill is the expensive part, so it runs only on the most
    promising candidates. Candidates are ranked by a cheap proxy (the
    sequence-budget imbalance the negative clusters already impose, plus the
    LLM share they put into the evaluation folds) and only the top `keep` per
    floor get a full fill and score.
    """
    T = totals(clusters)
    neg_cl = sorted((c for c, a in clusters.items() if a["neg"] > 0),
                    key=lambda c: -clusters[c]["neg"])
    rest = sorted((c for c in clusters if clusters[c]["neg"] == 0),
                  key=lambda c: -clusters[c]["nseq"])
    # Vector form for the inner loop.
    vec = [(clusters[c]["negseq"], clusters[c]["nseq"], clusters[c]["nlab"],
            clusters[c]["llm"]) for c in neg_cl]
    budget = {f: TARGET[f] * T["nseq"] for f in FOLDS}
    cap = {f: (1 + size_tol) * budget[f] for f in FOLDS}
    fi = {f: i for i, f in enumerate(FOLDS)}

    pools = {fl: [] for fl in floors}
    for combo in itertools.product((0, 1, 2), repeat=len(vec)):
        negseq = [0, 0, 0]; nseq = [0, 0, 0]; nlab = [0, 0, 0]; llm = [0, 0, 0]
        ok = True
        for k, f in enumerate(combo):
            n, s, L, m = vec[k]
            negseq[f] += n; nseq[f] += s; nlab[f] += L; llm[f] += m
        for f in FOLDS:
            if nseq[fi[f]] > cap[f]:
                ok = False
                break
        if not ok or negseq[fi["val"]] < floor_val:
            continue
        t = negseq[fi["test"]]
        proxy = (sum(abs(nseq[fi[f]] - budget[f]) for f in FOLDS) / T["nseq"]
                 + 1.5 * (llm[fi["test"]] + llm[fi["val"]])
                 / max(1, nlab[fi["test"]] + nlab[fi["val"]]))
        for fl in floors:
            if t >= fl:
                pools[fl].append((proxy, combo))

    best = {}
    for fl, pool in pools.items():
        if not pool:
            best[fl] = None
            continue
        pool.sort(key=lambda x: x[0])
        top = pool[:keep]
        bb = None
        for _, combo in top:
            full = {c: FOLDS[f] for c, f in zip(neg_cl, combo)}
            cur = collections.Counter()
            for c, f in full.items():
                cur[f] += clusters[c]["nseq"]
            for c in rest:
                f = max(FOLDS, key=lambda g: budget[g] - cur[g])
                full[c] = f
                cur[f] += clusters[c]["nseq"]
            score, agg = evaluate(full, clusters, T)
            if bb is None or score < bb[0]:
                bb = (score, full, agg)
        best[fl] = bb
    return best, T

def search(clusters, members, floor_test, floor_val, size_tol=0.5):
    T = totals(clusters)
    neg_cl = sorted((c for c, a in clusters.items() if a["neg"] > 0),
                    key=lambda c: -clusters[c]["neg"])
    rest = sorted((c for c in clusters if clusters[c]["neg"] == 0),
                  key=lambda c: -clusters[c]["nseq"])
    best = None
    for combo in itertools.product(FOLDS, repeat=len(neg_cl)):
        assign = dict(zip(neg_cl, combo))
        negseq = collections.Counter()
        nseq = collections.Counter()
        for c, f in assign.items():
            negseq[f] += clusters[c]["negseq"]
            nseq[f] += clusters[c]["nseq"]
        if negseq["test"] < floor_test or negseq["val"] < floor_val:
            continue
        # A fold cannot already be over its sequence budget before filling.
        if any(nseq[f] > (1 + size_tol) * TARGET[f] * T["nseq"] for f in FOLDS):
            continue
        full = dict(assign)
        cur = dict(nseq)
        for c in rest:
            f = max(FOLDS, key=lambda g: TARGET[g] * T["nseq"] - cur.get(g, 0))
            full[c] = f
            cur[f] = cur.get(f, 0) + clusters[c]["nseq"]
        score, agg = evaluate(full, clusters, T)
        if best is None or score < best[0]:
            best = (score, full, agg)
    return best, T

def report(agg, T, label=""):
    print(label)
    hdr = ("fold", "seqs", "%", "labels", "%", "neg", "negseq", "halo_seq", "llm%")
    print("  %-6s %5s %6s %7s %6s %5s %7s %9s %6s" % hdr)
    for f in FOLDS:
        a = agg[f]
        print("  %-6s %5d %5.1f%% %7d %5.1f%% %5d %7d %9d %5.1f%%" % (
            f, a["nseq"], 100 * a["nseq"] / T["nseq"], a["nlab"],
            100 * a["nlab"] / T["nlab"], a["neg"], a["negseq"], a["haloseq"],
            100 * a["llm"] / max(1, a["nlab"])))

def main():
    src = sys.argv[1]
    labels, seqs, clusters, members = load(src)
    T = totals(clusters)
    print("clusters %d | sequences %d | labels %d | negatives %d | neg-bearing seqs %d"
          % (len(clusters), T["nseq"], T["nlab"], T["neg"], T["negseq"]))
    print("negative-bearing clusters: %d\n" % sum(1 for a in clusters.values() if a["neg"] > 0))

    print("=== frontier: best achievable at each test inactive-sequence floor ===")
    floors = (10, 20, 30, 40, 50, 60, 70)
    best, _ = search_all(clusters, members, floors)
    frontier = []
    for fl in floors:
        b = best[fl]
        if b is None:
            print("floor_test=%-3d INFEASIBLE" % fl)
            continue
        score, full, agg = b
        frontier.append((fl, score, full, agg))
        print("floor_test=%-3d score=%.4f  test: %d seqs (%.1f%%), %d labels (%.1f%%), "
              "%d negseq, llm %.1f%%" % (
                  fl, score, agg["test"]["nseq"],
                  100 * agg["test"]["nseq"] / T["nseq"], agg["test"]["nlab"],
                  100 * agg["test"]["nlab"] / T["nlab"], agg["test"]["negseq"],
                  100 * agg["test"]["llm"] / max(1, agg["test"]["nlab"])))
    print()
    seen = set()
    for fl, score, full, agg in frontier:
        key = tuple(sorted(full.items()))
        if key in seen:
            continue
        seen.add(key)
        report(agg, T, "--- floor_test=%d (score %.4f) ---" % (fl, score))
        print()

    if len(sys.argv) > 2:
        out = sys.argv[2]
        os.makedirs(out, exist_ok=True)
        json.dump([{"floor_test": f, "score": s,
                    "assignment": a,
                    "folds": {k: dict(v) for k, v in g.items()}}
                   for f, s, a, g in frontier],
                  open(os.path.join(out, "frontier.json"), "w"), indent=2)
        print("frontier written to", os.path.join(out, "frontier.json"))

if __name__ == "__main__":
    main()
