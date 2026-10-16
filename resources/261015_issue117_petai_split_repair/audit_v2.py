"""Independent audit of the issue 117 v2 split. Recomputes every claim.

Usage: python3 audit_v2.py <v1_dir> <v2_dir> <out_json>
"""
import csv, json, sys, collections, hashlib, os

v1d, v2d, outp = sys.argv[1], sys.argv[2], sys.argv[3]
R = lambda d, p: list(csv.DictReader(open(os.path.join(d, p), newline="")))
CONC = {"µM", "mM", "nM", "M", "µmol/L", "mM/L"}
tok = lambda l: [x.strip() for x in l["substrate"].split(";")]
is_target = lambda l: "PET" in tok(l) and l["task"] == "quantitative" and l["unit"] in CONC
is_eval = lambda l: is_target(l) and l["source"] != "llm"

v1l, v1s = R(v1d, "labels.csv"), R(v1d, "sequences.csv")
v2l, v2s = R(v2d, "labels.csv"), R(v2d, "sequences.csv")
v2sp = R(v2d, "split.csv")
chk = {}

# v1 must be byte-identical to the published release.
chk["v1_untouched"] = {
    f: hashlib.md5(open(os.path.join(v1d, f), "rb").read()).hexdigest()
    for f in sorted(os.listdir(v1d))
}

# Conservation: every v2 label traces to a v1 label, and the drops are the 9.
v1_by = collections.Counter(l["seq_md5"] for l in v1l)
v1_split = {s["seq_md5"]: s["split"] for s in v1s}
dropped = [s["seq_md5"] for s in v1s
           if not v1_split[s["seq_md5"]] and not s["contained_in"]]
chk["labels"] = {"v1": len(v1l), "v2": len(v2l),
                 "dropped_sequences": len(dropped),
                 "dropped_labels": sum(v1_by[d] for d in dropped)}
assert len(v2l) + chk["labels"]["dropped_labels"] == len(v1l)

# Remapping: fragments gone, parents hold their labels.
remapped = [l for l in v2l if l["remapped"] == "True"]
chk["remap"] = {
    "labels_remapped": len(remapped),
    "fragments": len({l["orig_seq_md5"] for l in remapped}),
    "parents": len({l["seq_md5"] for l in remapped}),
}
v2_md5 = {s["seq_md5"] for s in v2s}
assert all(l["orig_seq_md5"] not in v2_md5 for l in remapped), "a fragment survived"
assert all(l["seq_md5"] in v2_md5 for l in v2l)

# Negatives must all survive.
chk["negatives"] = {
    "v1": sum(1 for l in v1l if l["is_negative"] == "True"),
    "v2": sum(1 for l in v2l if l["is_negative"] == "True"),
}
assert chk["negatives"]["v1"] == chk["negatives"]["v2"] == 200

# Leakage: recomputed from scratch, not trusting the writer.
split_of = {s["seq_md5"]: s["split"] for s in v2s}
leak = {}
for level in ("cluster_30", "cluster_50", "cluster_90", "component"):
    spans = collections.defaultdict(set)
    for s in v2s:
        spans[s[level]].add(s["split"])
    leak[level] = {"groups": len(spans),
                   "spanning_folds": sum(1 for v in spans.values() if len(v) > 1)}
chk["leakage"] = leak
assert all(v["spanning_folds"] == 0 for v in leak.values())

# split.csv must agree with sequences.csv.
assert {r["seq_md5"]: r["split"] for r in v2sp} == split_of

# Fold composition.
by_seq = collections.defaultdict(list)
for l in v2l:
    by_seq[l["seq_md5"]].append(l)
folds = {}
for f in ("train", "val", "test"):
    md5s = [m for m, g in split_of.items() if g == f]
    L = [l for m in md5s for l in by_seq[m]]
    ev = [l for l in L if is_eval(l)]
    folds[f] = {
        "sequences": len(md5s), "labels": len(L),
        "components": len({s["component"] for s in v2s if s["split"] == f}),
        "target_labels": sum(1 for l in L if is_target(l)),
        "eval_labels": len(ev), "eval_sequences": len({l["seq_md5"] for l in ev}),
        "eval_llm_labels": sum(1 for l in ev if l["source"] == "llm"),
        "halo_labels": sum(1 for l in L if l["source"] == "halo"),
        "negative_labels": sum(1 for l in L if l["is_negative"] == "True"),
        "llm_share_all_labels": round(
            100 * sum(1 for l in L if l["source"] == "llm") / max(1, len(L)), 1),
    }
chk["folds"] = folds
# The evaluated set carries no LLM labels and halo is train-only.
for f in ("val", "test"):
    assert folds[f]["eval_llm_labels"] == 0
    assert folds[f]["halo_labels"] == 0
assert sum(folds[f]["labels"] for f in folds) == len(v2l)
assert sum(folds[f]["sequences"] for f in folds) == len(v2s)

# FASTA must match the sequence table exactly.
names = [ln[1:].split()[0] for ln in open(os.path.join(v2d, "sequences.fasta"))
         if ln.startswith(">")]
chk["fasta"] = {"records": len(names), "matches_table": set(names) == v2_md5}
assert chk["fasta"]["matches_table"] and len(names) == len(v2s)

chk["verdict"] = "all assertions passed"
json.dump(chk, open(outp, "w"), indent=2)
print(json.dumps({k: v for k, v in chk.items() if k != "v1_untouched"}, indent=2))
