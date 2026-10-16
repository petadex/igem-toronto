"""Remap contained fragments onto their parent sequences.

Issue 117 step 1. Reads the v1 PETai label release and writes a new release
alongside it. The v1 files are never modified.

Of the 82 sequences with no split assignment, 73 have a contained_in parent
that already sits in a fold. Their labels are reattached to the parent, which
recovers 3,150 of the 3,188 stranded labels. The remaining 9 sequences have no
parent and are dropped with their 38 labels.

Lineage is preserved, not discarded: every remapped label keeps its original
md5 in orig_seq_md5 and carries remapped=True.

Usage: python3 remap_fragments.py <v1_dir> <out_dir>
"""
import csv, json, sys, collections, os

def read(path):
    with open(path, newline="") as fh:
        r = csv.DictReader(fh)
        return list(r), r.fieldnames

def write(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main(src, out):
    os.makedirs(out, exist_ok=True)
    labels, lab_f = read(os.path.join(src, "labels.csv"))
    seqs, seq_f = read(os.path.join(src, "sequences.csv"))
    split, split_f = read(os.path.join(src, "split.csv"))
    clusters, clu_f = read(os.path.join(src, "clusters.csv"))

    split_of = {r["seq_md5"]: r["split"] for r in split}
    by_md5 = {s["seq_md5"]: s for s in seqs}

    # Partition the unsplit sequences.
    remap, dropped = {}, []
    for s in seqs:
        if split_of.get(s["seq_md5"]):
            continue
        parent = s["contained_in"]
        if parent and split_of.get(parent):
            remap[s["seq_md5"]] = parent
        else:
            dropped.append(s["seq_md5"])

    # A parent must itself be a full sequence, not another fragment.
    for child, parent in remap.items():
        assert by_md5[parent].get("is_fragment") != "True", (child, parent)

    # Rewrite the labels.
    new_labels = []
    n_remapped = n_dropped_lab = 0
    for l in labels:
        md5 = l["seq_md5"]
        if md5 in dropped:
            n_dropped_lab += 1
            continue
        row = dict(l)
        if md5 in remap:
            row["seq_md5"] = remap[md5]
            row["orig_seq_md5"] = md5
            row["remapped"] = "True"
            n_remapped += 1
        else:
            row["orig_seq_md5"] = md5
            row["remapped"] = "False"
        new_labels.append(row)
    new_lab_f = lab_f + ["orig_seq_md5", "remapped"]

    # Recount per-sequence aggregates on the surviving sequences.
    counts = collections.Counter(l["seq_md5"] for l in new_labels)
    negs = collections.Counter(l["seq_md5"] for l in new_labels if l["is_negative"] == "True")
    refs = collections.defaultdict(set)
    for l in new_labels:
        refs[l["seq_md5"]].add(l["ref"])

    new_seqs = []
    for s in seqs:
        md5 = s["seq_md5"]
        if md5 in remap or md5 in dropped:
            continue
        row = dict(s)
        row["n_labels"] = str(counts[md5])
        row["n_negative"] = str(negs[md5])
        row["n_refs"] = str(len(refs[md5]))
        row["n_absorbed_fragments"] = str(sum(1 for c, p in remap.items() if p == md5))
        new_seqs.append(row)
    new_seq_f = seq_f + ["n_absorbed_fragments"]

    surviving = {s["seq_md5"] for s in new_seqs}
    new_split = [r for r in split if r["seq_md5"] in surviving]
    new_clusters = [r for r in clusters if r["seq_md5"] in surviving]

    write(os.path.join(out, "labels.csv"), new_labels, new_lab_f)
    write(os.path.join(out, "sequences.csv"), new_seqs, new_seq_f)
    write(os.path.join(out, "split.csv"), new_split, split_f)
    write(os.path.join(out, "clusters.csv"), new_clusters, clu_f)

    # The map is the audit trail for the remap.
    write(os.path.join(out, "fragment_map.csv"),
          [{"fragment_md5": c, "parent_md5": p,
            "fragment_labels": str(sum(1 for l in labels if l["seq_md5"] == c)),
            "parent_split": split_of[p]} for c, p in sorted(remap.items())],
          ["fragment_md5", "parent_md5", "fragment_labels", "parent_split"])
    write(os.path.join(out, "dropped_sequences.csv"),
          [{"seq_md5": d, "n_labels": str(sum(1 for l in labels if l["seq_md5"] == d)),
            "reason": "no contained_in parent"} for d in sorted(dropped)],
          ["seq_md5", "n_labels", "reason"])

    # Subset the FASTA to the surviving sequences.
    kept = skipped = 0
    with open(os.path.join(src, "sequences.fasta")) as fh, \
         open(os.path.join(out, "sequences.fasta"), "w") as oh:
        keep = False
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name in surviving
                kept += keep
                skipped += not keep
            if keep:
                oh.write(line)

    def tok(l):
        return [x.strip() for x in l["substrate"].split(";")]

    summary = {
        "source_release": os.path.abspath(src),
        "issue": 117,
        "step": "1 remap contained fragments onto parents",
        "v1": {"sequences": len(seqs), "labels": len(labels), "split_rows": len(split)},
        "v2": {"sequences": len(new_seqs), "labels": len(new_labels),
               "split_rows": len(new_split), "fasta_records": kept},
        "fragments_remapped": len(remap),
        "labels_remapped": n_remapped,
        "parents_receiving": len(set(remap.values())),
        "sequences_dropped": len(dropped),
        "labels_dropped": n_dropped_lab,
        "recovered": {
            "pet_substrate_token": sum(1 for l in new_labels
                                       if l["remapped"] == "True" and "PET" in tok(l)),
            "melting": sum(1 for l in new_labels if l["remapped"] == "True"
                           and "melt" in (l["variable"] + l["task"]).lower()),
            "negatives": sum(1 for l in new_labels
                             if l["remapped"] == "True" and l["is_negative"] == "True"),
        },
        "negatives_total": sum(1 for l in new_labels if l["is_negative"] == "True"),
        "labels_by_split": dict(collections.Counter(
            split_of[l["seq_md5"]] for l in new_labels)),
        "negatives_by_split": dict(collections.Counter(
            split_of[l["seq_md5"]] for l in new_labels if l["is_negative"] == "True")),
    }

    # Invariants. Every label is accounted for and every sequence resolves.
    assert len(new_labels) + n_dropped_lab == len(labels)
    assert all(l["seq_md5"] in surviving for l in new_labels)
    assert len(new_seqs) == len(seqs) - len(remap) - len(dropped)
    assert kept == len(new_seqs)

    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
