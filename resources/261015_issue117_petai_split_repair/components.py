"""Build leakage-safe assignment units for the PETai re-split.

Issue 117 step 2, prerequisite.

The release ships three clusterings (90%, 50%, 30% identity) that are not
nested: 14 of the 318 50%-identity clusters straddle two different 30%
clusters, so they were produced by independent runs rather than a cascade.
Splitting on cluster_30 alone therefore leaves sequences that share a 50%
cluster, and so are more similar than the 30% criterion implies, sitting in
different folds. The v1 release has four such pairs.

The fix is to assign connected components of the graph that links any two
sequences sharing a cluster at any threshold. A component is clean at every
threshold simultaneously, by construction.

Usage: python3 components.py <seq_csv> <out_csv>
"""
import csv, sys, collections

LEVELS = ("cluster_90", "cluster_50", "cluster_30")

def main(src, out):
    seqs = list(csv.DictReader(open(src, newline="")))
    parent = {s["seq_md5"]: s["seq_md5"] for s in seqs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for level in LEVELS:
        groups = collections.defaultdict(list)
        for s in seqs:
            groups[s[level]].append(s["seq_md5"])
        for members in groups.values():
            for m in members[1:]:
                union(members[0], m)

    comps = collections.defaultdict(list)
    for s in seqs:
        comps[find(s["seq_md5"])].append(s["seq_md5"])

    # Name each component by its lexicographically smallest member, so the id
    # is stable across reruns.
    comp_id = {}
    for members in comps.values():
        cid = min(members)
        for m in members:
            comp_id[m] = cid

    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["seq_md5", "component"])
        for s in seqs:
            w.writerow([s["seq_md5"], comp_id[s["seq_md5"]]])

    sizes = sorted((len(v) for v in comps.values()), reverse=True)
    n30 = len({s["cluster_30"] for s in seqs})
    print("sequences %d" % len(seqs))
    print("cluster_30 groups %d -> components %d (merged %d)"
          % (n30, len(comps), n30 - len(comps)))
    print("component sizes: max %d, top10 %s" % (sizes[0], sizes[:10]))

    # Confirm the components really are clean at every level.
    for level in LEVELS:
        spans = collections.defaultdict(set)
        for s in seqs:
            spans[s[level]].add(comp_id[s["seq_md5"]])
        bad = sum(1 for v in spans.values() if len(v) > 1)
        print("  %s groups spanning components: %d" % (level, bad))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
