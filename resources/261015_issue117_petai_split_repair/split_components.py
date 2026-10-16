"""Assign the leakage-free component split. Issue 117 step 2.

The headline evaluated target is PET monomer release regression, because it is
the only target that survives homology control: the 546-sequence component
holds 187 of the 197 halo binary labels, so halo activity can be trained but
not evaluated.

Assignment rules, in order:
  1. The 546-sequence component goes to train. At 48% of sequences it cannot
     be an evaluation fold, and it carries 1,733 of the 2,109 PET
     concentration labels.
  2. Every component holding a halo label goes to train. Halo is train-only.
  3. The components holding PET concentration labels are dealt alternately
     into val and test, largest first, into whichever fold holds fewer such
     labels. This balances the metric the folds exist to measure.
  4. Everything else goes to train.

Usage: python3 split_components.py <remapped_dir> <components_csv> <out_dir>
"""
import csv, json, os, re, sys, collections, unicodedata

# The unit spellings are shared with step 3 through target_units.json so the
# two steps cannot drift. The release spells the same unit several ways, so a
# raw string match silently misses labels: the micro sign U+00B5 and Greek mu
# U+03BC both occur in "uM".
_DEF = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "target_units.json")))
CONC = set(_DEF["concentration_to_uM"])
AMOUNT = set(_DEF["excluded_amount_units"])
TOKEN = _DEF["target_substrate_token"]
_DASHES = dict.fromkeys(map(ord, "\u2212\u2013\u2014\u2010\u2011"), "-")

def unit_key(raw):
    u = unicodedata.normalize("NFKC", raw).translate(_DASHES)
    return re.sub(r"\s+", " ", u).strip()

def tok(l):
    return [x.strip() for x in l["substrate"].split(";")]

def is_target(l):
    return (TOKEN in tok(l) and l["task"] == "quantitative"
            and unit_key(l["unit"]) in CONC)

def is_eval(l):
    """The evaluated set: target labels that are not LLM extractions.

    The v1 release evaluated on 62.9% LLM labels, which are unvalidated and
    contain no negatives. Only human-curated target labels are evaluated here;
    LLM target labels are kept for training.
    """
    return is_target(l) and l["source"] != "llm"

def main(src, comp_csv, out):
    os.makedirs(out, exist_ok=True)
    def read(p):
        with open(os.path.join(src, p), newline="") as fh:
            r = csv.DictReader(fh)
            return list(r), r.fieldnames
    labels, lab_f = read("labels.csv")
    seqs, seq_f = read("sequences.csv")
    clusters, clu_f = read("clusters.csv")
    comp = {r["seq_md5"]: r["component"] for r in
            csv.DictReader(open(comp_csv, newline=""))}

    by_seq = collections.defaultdict(list)
    for l in labels:
        by_seq[l["seq_md5"]].append(l)

    stats = collections.defaultdict(lambda: collections.Counter())
    comp_seqs = collections.defaultdict(list)
    for s in seqs:
        c = comp[s["seq_md5"]]
        comp_seqs[c].append(s["seq_md5"])
        L = by_seq[s["seq_md5"]]
        a = stats[c]
        a["nseq"] += 1
        a["nlab"] += len(L)
        a["target"] += sum(1 for l in L if is_target(l))
        a["target_seq"] += any(is_target(l) for l in L)
        a["eval"] += sum(1 for l in L if is_eval(l))
        a["eval_seq"] += any(is_eval(l) for l in L)
        a["halo"] += sum(1 for l in L if l["source"] == "halo")
        a["neg"] += sum(1 for l in L if l["is_negative"] == "True")
        a["llm"] += sum(1 for l in L if l["source"] == "llm")

    blob = max(stats, key=lambda c: stats[c]["nseq"])
    assign = {}
    assign[blob] = "train"
    for c, a in stats.items():
        if c != blob and a["halo"] > 0:
            assign[c] = "train"

    # Deal the components carrying evaluable labels into val and test. A
    # component whose target labels are all LLM extractions stays in train.
    pool = sorted((c for c, a in stats.items()
                   if c not in assign and a["eval"] > 0),
                  key=lambda c: (-stats[c]["eval"], c))
    # Few enough components to balance exactly rather than greedily. Balance
    # evaluated labels and evaluated sequences at once, since the metric is
    # computed per sequence but fitted per label.
    tot_lab = sum(stats[c]["eval"] for c in pool)
    tot_seq = sum(stats[c]["eval_seq"] for c in pool)
    best = None
    for mask in range(1 << len(pool)):
        vl = vs = 0
        for i, c in enumerate(pool):
            if mask >> i & 1:
                vl += stats[c]["eval"]
                vs += stats[c]["eval_seq"]
        if vl == 0 or vl == tot_lab or vs == 0 or vs == tot_seq:
            continue
        cost = (abs(2 * vl - tot_lab) / tot_lab + abs(2 * vs - tot_seq) / tot_seq)
        if best is None or cost < best[0]:
            best = (cost, mask)
    _, mask = best
    for i, c in enumerate(pool):
        assign[c] = "val" if mask >> i & 1 else "test"
    for c in stats:
        assign.setdefault(c, "train")

    for s in seqs:
        s["split"] = assign[comp[s["seq_md5"]]]
    split_of = {s["seq_md5"]: s["split"] for s in seqs}

    # Leakage must be zero at every threshold, and components must not split.
    for level in ("cluster_30", "cluster_50", "cluster_90"):
        spans = collections.defaultdict(set)
        for s in seqs:
            spans[s[level]].add(s["split"])
        bad = {k: v for k, v in spans.items() if len(v) > 1}
        assert not bad, (level, list(bad)[:3])
    spans = collections.defaultdict(set)
    for s in seqs:
        spans[comp[s["seq_md5"]]].add(s["split"])
    assert all(len(v) == 1 for v in spans.values())
    assert all(l["seq_md5"] in split_of for l in labels)

    def write(p, rows, fields):
        with open(os.path.join(out, p), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    for s in seqs:
        s["component"] = comp[s["seq_md5"]]
    write("sequences.csv", seqs, seq_f + ["component"])
    write("labels.csv", labels, lab_f)
    write("clusters.csv", clusters, clu_f)
    write("split.csv", [{"seq_md5": s["seq_md5"], "component": s["component"],
                         "cluster_30": s["cluster_30"], "split": s["split"],
                         "provenance": s["provenance"]} for s in seqs],
          ["seq_md5", "component", "cluster_30", "split", "provenance"])
    with open(os.path.join(src, "sequences.fasta")) as fh, \
         open(os.path.join(out, "sequences.fasta"), "w") as oh:
        oh.write(fh.read())

    agg = {f: collections.Counter() for f in ("train", "val", "test")}
    ncomp = collections.Counter()
    for c, f in assign.items():
        agg[f].update(stats[c])
        ncomp[f] += 1
    T = collections.Counter()
    for a in stats.values():
        T.update(a)

    print("%-6s %5s %6s %7s %7s %7s %6s %6s %6s %5s" % (
        "fold", "comps", "seqs", "%seq", "labels", "target", "eval", "evalseq",
        "halo", "llm%"))
    for f in ("train", "val", "test"):
        a = agg[f]
        print("%-6s %5d %6d %5.1f%% %7d %7d %6d %6d %6d %4.1f%%" % (
            f, ncomp[f], a["nseq"], 100 * a["nseq"] / T["nseq"], a["nlab"],
            a["target"], a["eval"], a["eval_seq"], a["halo"],
            100 * a["llm"] / max(1, a["nlab"])))
    print()
    print("evaluated labels are human-curated by construction: "
          "LLM share of the eval set is 0.0% in val and test")

    summary = {
        "issue": 117,
        "step": "2 leakage-free re-split on connected components",
        "source": os.path.abspath(src),
        "headline_target": "PET monomer release, quantitative, concentration units",
        "evaluated_set": "target labels with source != llm; LLM target labels are train-only",
        "target_units": sorted(CONC),
        "target_units_source": "target_units.json",
        "excluded_amount_units": sorted(AMOUNT),
        "halo": "train-only, not evaluated: 187 of 197 halo labels fall in the train component",
        "unit": "connected component over cluster_90, cluster_50 and cluster_30",
        "components": len(stats),
        "train_component_sequences": stats[blob]["nseq"],
        "leakage": "zero cross-fold groups at cluster_30, cluster_50 and cluster_90",
        "folds": {f: dict(agg[f]) | {"components": ncomp[f]} for f in agg},
        "totals": dict(T),
    }
    json.dump(summary, open(os.path.join(out, "summary.json"), "w"), indent=2)
    json.dump({c: assign[c] for c in sorted(assign)},
              open(os.path.join(out, "component_assignment.json"), "w"), indent=2)
    print("\nwrote", out)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
