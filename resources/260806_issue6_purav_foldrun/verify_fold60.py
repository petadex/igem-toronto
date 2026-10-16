"""Per-ORF completeness check for the 60% fold: every input sequence of each tier has an ok row in the
progress CSVs, a structures/<NNN>/<id>.cif.gz and a metrics/<NNN>/<id>.json.
    python verify_fold60.py            # L40S tiers + H100 tail (whatever exists)"""
import csv, glob, os
TIERS = {"short_full": "short_full.fa", "mid_fast": "mid_fast.fa", "tail_h100": "tail_h100_le2500.fa"}
for t, fa in TIERS.items():
    if not os.path.isdir(f"data/fold60pid/{t}"):
        print(f"{t}: not started"); continue
    ids=[l[1:].strip() for l in open(f"data/fold60pid_input/{fa}") if l.startswith(">")]
    ok=set(); bad={}
    for p in glob.glob(f"data/fold60pid/{t}/progress*.csv"):
        for r in csv.DictReader(open(p)):
            if r["status"]=="ok": ok.add(r["id"])
            else: bad[r["id"]]=r["status"]
    missing_prog=[i for i in ids if i not in ok]
    missing_file=[]; missing_json=[]
    for i in ids:
        sub=f"{int(i[3:])%1000:03d}"
        if not os.path.exists(f"data/fold60pid/{t}/structures/{sub}/{i}.cif.gz"): missing_file.append(i)
        if not os.path.exists(f"data/fold60pid/{t}/metrics/{sub}/{i}.json"): missing_json.append(i)
    print(f"{t}: input {len(ids):,} | ok in progress {len(ok & set(ids)):,} | not ok {len(missing_prog)} | "
          f"missing .cif.gz {len(missing_file)} | missing .json {len(missing_json)} | non-ok rows {len(bad)} {list(bad.items())[:3]}")
    if missing_file[:5]: print("   e.g.", missing_file[:5])
