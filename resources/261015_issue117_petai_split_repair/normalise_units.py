"""Normalise units and name the regression targets explicitly. Issue 117 step 3.

The release spells units 621 different ways across 13,054 labels. 15 of those
are pure encoding duplicates: NFKC folds the micro sign U+00B5 into Greek mu
and the superscript minus into a hyphen, which takes 621 to 590. Lowercasing
would take it to 575 and is wrong, because nm and nM, mm and mM, and um and uM
all occur in the data and mean different things.

The rest is not a normalisation problem but a target-definition problem. This
script names seven targets, each with one canonical unit and an explicit map of
the spellings that convert into it. A label that does not map is left
unassigned rather than coerced, and anything suspicious is flagged rather than
silently rescaled:

  ambiguous_fraction  a percent target whose value is in (0, 1], so it may be
                      a fraction recorded as if it were a percentage
  out_of_range        a melting temperature outside 20 to 120 C, which the
                      degC column contains (it reaches 401.1) and the K column
                      does too (it contains an 8.8)
  mw_required         a Km in g/L or mg/mL, needing a molecular weight
  scaled_unit         a unit carrying its own scale factor, such as
                      "s-1 M-1 x 10-6", where guessing the direction is unsafe

Usage: python3 normalise_units.py <v2_dir> <out_dir>
"""
import csv, json, os, re, sys, collections, unicodedata

DASHES = dict.fromkeys(map(ord, "−–—‐‑"), "-")

def unit_key(raw):
    """NFKC fold, unify dashes, collapse whitespace. Case is left alone."""
    u = unicodedata.normalize("NFKC", raw).translate(DASHES)
    return re.sub(r"\s+", " ", u).strip()

# Concentration spellings, in micromolar. Shared with step 2 so the evaluated
# target cannot drift between the split and the normalisation.
_DEF = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "target_units.json")))
CONC = {k: float(v) for k, v in _DEF["concentration_to_uM"].items()}
TOKEN = _DEF["target_substrate_token"]
# Turnover spellings, in reciprocal seconds.
RATE = {"1/s": 1.0, "s-1": 1.0, "sec-1": 1.0, "1/min": 1 / 60, "min-1": 1 / 60}
# Catalytic efficiency, in reciprocal molar reciprocal seconds.
EFF = {"1/s/M": 1.0, "M-1s-1": 1.0, "M-1 s-1": 1.0,
       "mM-1s-1": 1e3, "mM-1 s-1": 1e3, "μM-1s-1": 1e6, "μM-1 s-1": 1e6}
PERCENT = {"%": 1.0, "percent": 1.0, "% relative activity": 1.0,
           "percent weight loss": 1.0, "% (w/w)": 1.0}
MW_REQUIRED = {"g/L", "mg/mL", "mg/L", "g/l", "mg/ml"}

SCALED = re.compile(r"10\s*-?\d|×\s*10|x\s*10")

def has(variable, *needles):
    v = variable.lower()
    return any(n.lower() in v for n in needles)

def classify(row, key):
    """Return (target, canonical_unit, factor, offset, flag). factor None = unassigned."""
    var, task = row["variable"], row["task"]
    substrates = [x.strip() for x in row["substrate"].split(";")]

    if task == "binary":
        return "halo_binary", "z vs empty vector", 1.0, 0.0, ""

    # Melting temperature. The CD "fraction of folded protein" rows are curve
    # points, not a Tm, and are deliberately excluded.
    if has(var, "melting temperature", "tm ") or var.strip() == "Tm":
        if key == "°C":
            return "melting_temperature", "°C", 1.0, 0.0, ""
        if key == "K":
            return "melting_temperature", "°C", 1.0, -273.15, ""
        return None, None, None, None, ""

    if SCALED.search(key):
        return None, None, None, None, "scaled_unit"

    # Kinetics, before the generic concentration rule, since Km is a
    # concentration and would otherwise be read as a release measurement.
    if has(var, "kcat/km", "catalytic efficiency"):
        if key in EFF:
            return "kcat_over_km", "M-1 s-1", EFF[key], 0.0, ""
        return None, None, None, None, ""
    if has(var, "kcat"):
        if key in RATE:
            return "kcat", "s-1", RATE[key], 0.0, ""
        return None, None, None, None, ""
    if var.strip() == "Km" or has(var, "km ("):
        if key in CONC:
            return "km", "μM", CONC[key], 0.0, ""
        if key in MW_REQUIRED:
            return None, None, None, None, "mw_required"
        return None, None, None, None, ""

    # Mass loss is a percentage but a different quantity from activity, so it
    # gets its own target instead of being folded into relative activity.
    if has(var, "weight loss", "weight/mass loss", "mass loss"):
        if key in PERCENT:
            return "weight_loss", "percent", PERCENT[key], 0.0, ""
        return None, None, None, None, ""

    # The unit is more reliable than the variable here. 369 labels carry the
    # unit "% relative activity" while naming an assay in the variable, such as
    # pNP release, so the unit decides. The substrate stays on the row, since
    # relative activity on pNP is not the same measurement as on PET.
    if key == "% relative activity" or has(var, "relative activity",
                                           "residual activity"):
        if key in PERCENT:
            return "relative_activity", "percent", PERCENT[key], 0.0, ""
        return None, None, None, None, ""

    # PET monomer release. Concentration only: mmol and umol are amounts and
    # need an assay volume, so they stay unassigned.
    if TOKEN in substrates and key in CONC:
        return "pet_release_conc", "μM", CONC[key], 0.0, ""

    return None, None, None, None, ""

def main(src, out):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(src, "labels.csv"), newline="") as fh:
        r = csv.DictReader(fh)
        labels, fields = list(r), r.fieldnames
    split_of = {s["seq_md5"]: s["split"] for s in
                csv.DictReader(open(os.path.join(src, "sequences.csv"), newline=""))}

    raw_units = {l["unit"] for l in labels}
    keys = {unit_key(u) for u in raw_units}

    rows = []
    for l in labels:
        key = unit_key(l["unit"])
        target, cunit, factor, offset, flag = classify(l, key)
        value = float(l["value"])
        canon = ""
        if factor is not None:
            canon = value * factor + offset
            if target == "melting_temperature" and not 20.0 <= canon <= 120.0:
                target, canon, flag = None, "", "out_of_range"
            elif cunit == "percent" and 0 < value <= 1:
                flag = "ambiguous_fraction"
        rows.append(dict(l, unit_key=key, target=target or "",
                         value_canonical="" if canon == "" else repr(round(canon, 6)),
                         unit_canonical=cunit or "", flag=flag,
                         split=split_of[l["seq_md5"]]))

    new_fields = fields + ["unit_key", "target", "value_canonical",
                           "unit_canonical", "flag", "split"]
    with open(os.path.join(out, "labels_normalised.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=new_fields)
        w.writeheader()
        w.writerows(rows)

    assigned = [r for r in rows if r["target"]]
    per_target = collections.defaultdict(lambda: collections.Counter())
    for r in assigned:
        t = per_target[r["target"]]
        t["labels"] += 1
        t[r["split"]] += 1
        t["sequences"] = 0
    seqs_per = collections.defaultdict(set)
    for r in assigned:
        seqs_per[r["target"]].add(r["seq_md5"])

    print("units: %d raw spellings, %d after NFKC and dash folding" % (len(raw_units), len(keys)))
    print("labels assigned to a target: %d of %d (%.1f%%)"
          % (len(assigned), len(rows), 100 * len(assigned) / len(rows)))
    print()
    print("%-22s %7s %6s %7s %6s %6s %-12s" %
          ("target", "labels", "seqs", "train", "val", "test", "canonical"))
    for t in sorted(per_target, key=lambda t: -per_target[t]["labels"]):
        c = per_target[t]
        unit = next(r["unit_canonical"] for r in assigned if r["target"] == t)
        print("%-22s %7d %6d %7d %6d %6d %-12s"
              % (t, c["labels"], len(seqs_per[t]), c["train"], c["val"], c["test"], unit))
    print()
    flags = collections.Counter(r["flag"] for r in rows if r["flag"])
    print("flagged:", dict(flags) or "none")
    print()
    unassigned = collections.Counter(r["unit_key"] for r in rows if not r["target"])
    print("largest unassigned unit groups:")
    for u, n in unassigned.most_common(12):
        print("  %-32s %5d" % (repr(u), n))

    summary = {
        "issue": 117, "step": "3 normalise units and name targets",
        "raw_unit_spellings": len(raw_units),
        "unit_keys_after_folding": len(keys),
        "lowercasing_rejected": "nm/nM, mm/mM and μm/μM all occur and differ",
        "labels": len(rows), "assigned": len(assigned),
        "targets": {t: {"labels": per_target[t]["labels"],
                        "sequences": len(seqs_per[t]),
                        "canonical_unit": next(r["unit_canonical"] for r in assigned
                                               if r["target"] == t),
                        "by_split": {f: per_target[t][f] for f in ("train", "val", "test")}}
                    for t in per_target},
        "flags": dict(flags),
        "excluded_by_design": {
            "fraction of folded protein": "CD curve points, not a Tm",
            "mmol and μmol": "amounts, not concentrations; need an assay volume",
        },
    }
    json.dump(summary, open(os.path.join(out, "targets_summary.json"), "w"), indent=2)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
