"""Can the 546-sequence component be broken by dropping a few sequences?

The component exists because the three clusterings are not nested. If only a
handful of sequences create the straddles, dropping them restores nesting and
lets the component split back into pieces small enough to assign.

Enforces nesting by, for each finer group that straddles coarser groups,
keeping the members in the largest coarser group and dropping the rest. Then
recomputes components on the survivors.
"""
import csv, collections, sys

D = "/project/6122149/purav/petai-labels/remapped/"
seqs = list(csv.DictReader(open(D + "sequences.csv", newline="")))
labels = list(csv.DictReader(open(D + "labels.csv", newline="")))
by = collections.defaultdict(list)
for l in labels:
    by[l["seq_md5"]].append(l)

LEVELS = ("cluster_90", "cluster_50", "cluster_30")

def components(rows):
    parent = {s["seq_md5"]: s["seq_md5"] for s in rows}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for level in LEVELS:
        g = collections.defaultdict(list)
        for s in rows:
            g[s[level]].append(s["seq_md5"])
        for m in g.values():
            r0 = find(m[0])
            for x in m[1:]:
                rx = find(x)
                if rx != r0:
                    parent[rx] = r0
    comps = collections.defaultdict(list)
    for s in rows:
        comps[find(s["seq_md5"])].append(s["seq_md5"])
    return comps

drop = set()
# Enforce nesting from the finest pair upward.
for fine, coarse in (("cluster_90", "cluster_50"), ("cluster_50", "cluster_30"),
                     ("cluster_90", "cluster_30")):
    g = collections.defaultdict(lambda: collections.defaultdict(list))
    for s in seqs:
        if s["seq_md5"] in drop:
            continue
        g[s[fine]][s[coarse]].append(s["seq_md5"])
    straddles = {k: v for k, v in g.items() if len(v) > 1}
    for k, v in straddles.items():
        keep = max(v.values(), key=len)
        for coarse_id, members in v.items():
            if members is not keep:
                drop.update(members)
    print("%s within %s: %d straddling groups" % (fine, coarse, len(straddles)))

print()
print("sequences dropped to restore nesting: %d of %d (%.2f%%)"
      % (len(drop), len(seqs), 100 * len(drop) / len(seqs)))
dl = sum(len(by[m]) for m in drop)
dn = sum(1 for m in drop for l in by[m] if l["is_negative"] == "True")
dh = sum(1 for m in drop for l in by[m] if l["source"] == "halo")
print("labels lost: %d of %d (%.2f%%) | negatives lost: %d of 200 | halo lost: %d of 197"
      % (dl, len(labels), 100 * dl / len(labels), dn, dh))

survivors = [s for s in seqs if s["seq_md5"] not in drop]
comps = components(survivors)
sizes = sorted((len(v) for v in comps.values()), reverse=True)
print()
print("after dropping: %d sequences, %d components, sizes top10 %s"
      % (len(survivors), len(comps), sizes[:10]))

cid = {m: min(v) for v in comps.values() for m in v}
agg = collections.defaultdict(lambda: collections.Counter())
for s in survivors:
    c = agg[cid[s["seq_md5"]]]
    L = by[s["seq_md5"]]
    c["nseq"] += 1
    c["nlab"] += len(L)
    n = sum(1 for l in L if l["is_negative"] == "True")
    c["neg"] += n
    c["negseq"] += n > 0
    c["halo"] += sum(1 for l in L if l["source"] == "halo")
tot = collections.Counter()
for c in agg.values():
    tot.update(c)
print("surviving negatives: %d on %d sequences | halo %d" % (tot["neg"], tot["negseq"], tot["halo"]))
print()
print("components holding negatives (nseq, nlab, neg, negseq, halo):")
for k, c in sorted(agg.items(), key=lambda x: -x[1]["neg"]):
    if c["neg"]:
        print("  %5d %6d %5d %6d %5d" % (c["nseq"], c["nlab"], c["neg"], c["negseq"], c["halo"]))
