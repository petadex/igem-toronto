#!/usr/bin/env python
"""Run record for the 60% centroid fold: timings, GPU usage, throughput, memory and incidents,
from SLURM accounting (sacct) and the per-shard progress CSVs. Writes markdown tables to
data/fold60pid/qc/run_record.md (pasted into notebook Step 16) plus TSVs for plotting.

    python fold60_run_record.py
"""
import csv
import glob
import os
import statistics as st
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timedelta

OUT = "data/fold60pid/qc"
NAMES = "fold60-short_full,fold60-mid_fast,fold60-tail"
TIER_OF = {"fold60-short_full": "short_full", "fold60-mid_fast": "mid_fast", "fold60-tail": "tail_h100"}
HW = {"short_full": "L40S", "mid_fast": "L40S", "tail_h100": "H100"}
BAD_NODES = ["kn140", "kn161", "kn088", "kn117", "kn112"]
# from the shard logs: "CUDA is not available" vs "cudaErrorECCUncorrectable"
CAUSE = {"kn140": "CUDA unavailable (GPU visible to nvidia-smi)", "kn088": "CUDA unavailable (GPU visible to nvidia-smi)",
         "kn112": "CUDA unavailable (GPU visible to nvidia-smi)",
         "kn161": "uncorrectable ECC GPU memory error (also during the Phase 2 sweep)",
         "kn117": "uncorrectable ECC GPU memory error"}


def sacct():
    cols = "JobID,JobName,State,Submit,Start,End,ElapsedRaw,NodeList"
    out = subprocess.run(["sacct", "-X", "-n", "-P", "-u", os.environ.get("USER", "purav"),
                          "-S", "2026-10-15T22:00", "--name", NAMES, "-o", cols],
                         capture_output=True, text=True).stdout
    rows = []
    for line in out.strip().split("\n"):
        f = line.split("|")
        if len(f) != 8:
            continue
        r = dict(zip(cols.split(","), f))
        r["State"] = r["State"].split()[0]
        for k in ("Submit", "Start", "End"):
            try:
                r[k] = datetime.fromisoformat(r[k])
            except ValueError:
                r[k] = None
        r["ElapsedRaw"] = int(r["ElapsedRaw"] or 0)
        r["tier"] = TIER_OF.get(r["JobName"], "?")
        rows.append(r)
    return rows


def hm(sec):
    sec = int(sec)
    return f"{sec // 3600}h{(sec % 3600) // 60:02d}m"


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * (len(xs) - 1) + 0.5))]


def main():
    os.makedirs(OUT, exist_ok=True)
    J = sacct()
    md = []
    w = md.append

    # ---------------- timeline ----------------
    starts = [r["Start"] for r in J if r["Start"]]
    ends = [r["End"] for r in J if r["End"]]
    first_submit = min(r["Submit"] for r in J if r["Submit"])
    by_tier_end = defaultdict(list)
    for r in J:
        if r["End"] and r["State"] == "COMPLETED":
            by_tier_end[r["tier"]].append(r["End"])
    w("#### Timeline\n")
    w("| event | time |\n|---|---|")
    w(f"| first submission (plumbing test) | {first_submit:%Y-%m-%d %H:%M} |")
    w(f"| first shard started | {min(starts):%Y-%m-%d %H:%M} |")
    for t in ("short_full", "mid_fast", "tail_h100"):
        if by_tier_end[t]:
            w(f"| last {t} shard finished | {max(by_tier_end[t]):%Y-%m-%d %H:%M} |")
    w(f"| wall-clock, first start to last finish | {hm((max(ends) - min(starts)).total_seconds())} |\n")

    # ---------------- GPU-hours and outcomes per tier ----------------
    w("#### GPU-hours and job outcomes\n")
    w("| tier | hardware | state | jobs | GPU-h |\n|---|---|---|---|---|")
    tot = defaultdict(float)
    for (t, s), rs in sorted(defaultdict(list, {k: v for k, v in
            ((k, [r for r in J if (r['tier'], r['State']) == k]) for k in {(r['tier'], r['State']) for r in J})}).items()):
        h = sum(r["ElapsedRaw"] for r in rs) / 3600
        tot[t] += h
        w(f"| {t} | {HW.get(t, '?')} | {s} | {len(rs)} | {h:.1f} |")
    w(f"| **total** | | | **{len(J)}** | **{sum(tot.values()):.1f}** |")
    w(f"\nL40S {tot['short_full'] + tot['mid_fast']:.1f} GPU-h, H100 {tot['tail_h100']:.1f} GPU-h.\n")

    # ---------------- shard run times and queue waits (completed only) ----------------
    w("#### Shard run time and queue wait (completed jobs)\n")
    w("| tier | completed jobs | run time min / median / max | queue wait median / max |\n|---|---|---|---|")
    for t in ("short_full", "mid_fast", "tail_h100"):
        rs = [r for r in J if r["tier"] == t and r["State"] == "COMPLETED" and r["Start"] and r["Submit"]]
        if not rs:
            continue
        el = [r["ElapsedRaw"] for r in rs]
        wait = [(r["Start"] - r["Submit"]).total_seconds() for r in rs]
        w(f"| {t} | {len(rs)} | {hm(min(el))} / {hm(st.median(el))} / {hm(max(el))} | {hm(st.median(wait))} / {hm(max(wait))} |")
    w("")

    # ---------------- concurrency per hour ----------------
    t0 = min(starts).replace(minute=0, second=0, microsecond=0)
    t1 = max(ends)
    conc = []
    t = t0
    while t <= t1:
        mid = t + timedelta(minutes=30)
        run = defaultdict(int)
        for r in J:
            if r["Start"] and r["End"] and r["Start"] <= mid < r["End"] and r["ElapsedRaw"] > 120:
                run[HW.get(r["tier"], "?")] += 1
        conc.append((t, run["L40S"], run["H100"]))
        t += timedelta(hours=1)
    w("#### GPUs in use (sampled at the half hour)\n")
    w("| hour | L40S | H100 |\n|---|---|---|")
    for t, a, b in conc:
        w(f"| {t:%m-%d %H}:30 | {a} | {b} |")
    with open(f"{OUT}/run_concurrency.tsv", "w") as fh:
        fh.write("hour\tl40s\th100\n" + "".join(f"{t:%Y-%m-%dT%H:30}\t{a}\t{b}\n" for t, a, b in conc))
    w(f"\nPeak {max(c[1] for c in conc)} L40S; {max(c[2] for c in conc)} H100.\n")

    # ---------------- per-fold timings and memory from progress CSVs ----------------
    folds = defaultdict(list)  # tier -> (L, wall, vram, residues/s)
    for t in ("short_full", "mid_fast", "tail_h100"):
        for p in glob.glob(f"data/fold60pid/{t}/progress*.csv"):
            for r in csv.DictReader(open(p)):
                if r["status"] == "ok" and r["wall_s"]:
                    folds[t].append((int(r["seq_len"]), float(r["wall_s"]),
                                     float(r["peak_vram_gb"]) if r.get("peak_vram_gb") else None))
    w("#### Per-fold time and throughput\n")
    w("| tier | folds | fold time median (p10-p90) | residues folded | sum of fold time | folds per GPU-hour |\n|---|---|---|---|---|---|")
    eff = {}
    for t, fs in folds.items():
        wall = [x[1] for x in fs]
        fold_h = sum(wall) / 3600
        eff[t] = fold_h / tot[t] if tot[t] else None
        w(f"| {t} | {len(fs):,} | {st.median(wall):.2f} s ({q(wall, .1):.2f}-{q(wall, .9):.2f}) | "
          f"{sum(x[0] for x in fs) / 1e6:.1f} M | {fold_h:.1f} h | {len(fs) / tot[t]:,.0f} |")
    w("\n**End-to-end efficiency** (sum of fold time / GPU-hours billed, including model load, file writes, "
      "failed and timed-out jobs): " + ", ".join(f"{t} {v:.1%}" for t, v in eff.items() if v) + ".\n")

    bins = [(40, 100), (100, 150), (150, 250), (250, 400), (400, 600), (600, 900), (900, 1300),
            (1300, 1741), (1741, 2000), (2000, 2250), (2250, 2501)]
    w("#### Fold time and peak GPU memory by length\n")
    w("| length | tier | folds | fold time median | peak memory median / max |\n|---|---|---|---|---|")
    rows_tsv = []
    for lo, hi in bins:
        for t, fs in folds.items():
            s = [x for x in fs if lo <= x[0] < hi]
            if len(s) < 20:
                continue
            v = [x[2] for x in s if x[2] is not None]
            w(f"| {lo}-{hi - 1} aa | {t} | {len(s):,} | {st.median(x[1] for x in s):.2f} s | "
              f"{st.median(v):.1f} / {max(v):.1f} GB |")
            rows_tsv.append((lo, hi - 1, t, len(s), st.median(x[1] for x in s), st.median(v), max(v)))
    with open(f"{OUT}/run_length_profile.tsv", "w") as fh:
        fh.write("len_lo\tlen_hi\ttier\tfolds\tfold_s_median\tvram_gb_median\tvram_gb_max\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows_tsv)
    w("")

    # ---------------- H100 per shard ----------------
    w("#### H100 tail, per shard\n")
    w("| shard | job | node | sequences | length | run time | peak memory |\n|---|---|---|---|---|---|---|")
    for r in sorted((r for r in J if r["tier"] == "tail_h100"), key=lambda r: int(r["JobID"].split("_")[1])):
        k = int(r["JobID"].split("_")[1])
        p = f"data/fold60pid/tail_h100/progress.shard{k}of6.csv"
        ok = [x for x in csv.DictReader(open(p)) if x["status"] == "ok"] if os.path.exists(p) else []
        L = [int(x["seq_len"]) for x in ok]
        v = [float(x["peak_vram_gb"]) for x in ok]
        w(f"| {k} | {r['JobID']} | {r['NodeList']} | {len(ok)} | {min(L)}-{max(L)} aa | {hm(r['ElapsedRaw'])} | {max(v):.1f} GB |")
    w("")

    # ---------------- incidents ----------------
    w("#### Incidents\n")
    w("| node | failed jobs | first / last failure | median time to fail | cause |\n|---|---|---|---|---|")
    for n in BAD_NODES:
        rs = [r for r in J if r["NodeList"] == n and r["State"] == "FAILED"]
        if not rs:
            continue
        tf = [r["End"] for r in rs if r["End"]]
        w(f"| {n} | {len(rs)} | {min(tf):%m-%d %H:%M} / {max(tf):%m-%d %H:%M} | "
          f"{st.median(r['ElapsedRaw'] for r in rs):.0f} s | {CAUSE.get(n, '?')} |")
    other = [r for r in J if r["State"] == "FAILED" and r["NodeList"] not in BAD_NODES]
    to = [r for r in J if r["State"] == "TIMEOUT"]
    cx = [r for r in J if r["State"] == "CANCELLED"]
    w(f"\nOther failures: {len(other)}. Timeouts: {len(to)} ({', '.join(r['JobID'] for r in to)}), all resumed. "
      f"Cancelled: {len(cx)} (duplicate shard jobs, cancelled while pending; never ran).\n")

    open(f"{OUT}/run_record.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
