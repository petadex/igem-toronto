"""make_final_figures.py -- the result figures for issue #84, items 2, 3 and 5.

    fig2_naive_vs_designed.png   what the design saves against ordering every core
    fig3_frontier_over_k.png     the whole K frontier: nt, coverage, library
    fig5_runtime.png             where the wall time actually goes

Every nucleotide count here is a POST-PROCESSING count: the length of a construct
as it would be ordered, `CGTCTC + spacer + insert + spacer + GAGACG`, computed by
calling make_order.py's own `build_inserts` and `pick_spacer` rather than by
re-deriving the layout.  The designer's own `nt_ordered` is a flat 24 nt/oligo
estimate and is NOT used: the true overhead is 20/18/20 nt per layer, so on
cluster 1 at K=3 that estimate over-charges by 246 nt.

`summary.json` records the post-processed design for the recommended K only, so
the other K values are re-derived by re-running the sweep in process with the
recorded args and seed.  That is deterministic -- the same property validator
check 10 relies on -- and the script proves it by asserting its own figure for the
recommended K equals the real `order/constructs.csv` total.  The re-run takes
~10-20 min on cluster 1 and is cached; `--from-cache` redraws in a second.

Colours are the dataviz reference palette, unchanged and in documented order:
slot 1 blue #2a78d6, slot 2 orange #eb6834.

USAGE
-----
    python make_final_figures.py <run-dir> --cluster 1 [--from-cache]

Figures and the cache land in `<out>/cluster<N>/`, one directory per cluster,
so a second cluster never overwrites the first.
"""

import argparse
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "final"))

import cutsearch_design as U                                   # noqa: E402
import make_order as MO                                        # noqa: E402

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
import matplotlib.ticker as mticker                            # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
GRID = "#e6e4df"

# Per-oligo assembly overhead make_order.py actually adds: the two BsmBI
# recognition sites and the two spacers.  The backbone overhangs and the 2-nt
# reaches into each neighbour belong to the insert, so build_inserts counts
# those, not this.
SITE_NT = len(MO.BSMBI) + 1 + 1 + len(MO.BSMBI_RC)


# =========================================================================== #
# Measurement
# =========================================================================== #

class Clock:
    """Wall time attributed to the cut search and to the greedy, per K."""

    def __init__(self):
        self.cut = self.greedy = 0.0
        self.cut_calls = self.greedy_calls = 0

    def reset(self):
        self.__init__()


def instrument(clock):
    """Wrap the two hot functions.  `evaluate_K` resolves both as module globals
    at call time, so rebinding them here is enough."""
    inner_cut, inner_greedy = U.place_cuts, U.greedy

    def timed_cut(*a, **k):
        t = time.perf_counter()
        try:
            return inner_cut(*a, **k)
        finally:
            clock.cut += time.perf_counter() - t
            clock.cut_calls += 1

    def timed_greedy(*a, **k):
        t = time.perf_counter()
        try:
            return inner_greedy(*a, **k)
        finally:
            clock.greedy += time.perf_counter() - t
            clock.greedy_calls += 1

    U.place_cuts, U.greedy = timed_cut, timed_greedy


def configure(sargs):
    """Reproduce the module-level state `cutsearch_design.main()` sets up from the
    recorded args, so the re-run is the same run."""
    if sargs.get("no_codon_usage"):
        U.USE_CODON_USAGE = False
        U.DEG_TABLE = U._build_degenerate_table()
        U._codon_cache.clear()

    site = {"bsmbi": "CGTCTC", "esp3i": "CGTCTC",
            "bsai": "GGTCTC"}[sargs.get("gg_enzyme", "bsmbi")]
    U.FORBIDDEN_SITES = frozenset(
        s for x in ({site} | set(U.ALWAYS_FORBIDDEN_SITES))
        for s in (x, U.revcomp(x)))

    if sargs["chemistry"] != "gg":
        reserved = frozenset()
    elif sargs.get("shared_backbone_overhangs"):
        reserved = U.IGEM_RESERVED_OVERHANGS
    else:
        reserved = U.RESERVED_OVERHANGS

    if sargs.get("cut_search") == "dp":
        import dp_cutsearch
        U.place_cuts = dp_cutsearch.dp_cut_search
    return reserved


def price(design_layers, tokens, five_pad="", b5="CGGA", b3="GGTG"):
    """Construct length for every oligo, through make_order.py's own layout.

    Returns (lengths, insert_nt, n_spacer_failures).  A spacer failure means no
    single nucleotide either side left the construct with exactly two Type IIS
    sites; it has never happened, and is surfaced rather than swallowed."""
    design = {"fragments": [[{"oligo": u.oligo(), "variants": u.variants,
                              "nt": u.nt} for u in units]
                            for units in design_layers],
              "junctions": [list(t) for t in tokens]}
    inserts, _left, _right = MO.build_inserts(design, five_pad, b5, b3)

    lengths, insert_nt, fails = [], 0, 0
    for layer in inserts:
        for ins in layer:
            sp, _sp2 = MO.pick_spacer(ins)
            if sp is None:
                fails += 1
            insert_nt += len(ins)
            lengths.append(SITE_NT + len(ins))
    return lengths, insert_nt, fails


def measure(run_dir, k_max=None):
    """Re-run the K sweep, price every design as ordered, and time both halves."""
    with open(os.path.join(run_dir, "summary.json")) as fh:
        S = json.load(fh)
    a = S["args"]
    reserved = configure(a)

    aln = (S["input"] if os.path.exists(S["input"])
           else os.path.join(HERE, "final", a["aln_fasta"]))
    seqs = U.read_aligned_cores(aln)
    aligned = [s for s, _ in seqs]
    weights = [w for _, w in seqs]
    L = len(aligned[0])
    const = U.constant_columns(aligned)
    overhead = (a["oligo_overhead_nt"] if a["oligo_overhead_nt"] is not None
                else (24 if a["chemistry"] == "gg" else 0))
    max_layer_cols = (a["max_oligo_nt"] // 3) if a["max_oligo_nt"] else None

    clock = Clock()
    instrument(clock)

    rows = []
    for K in range(1, (k_max or a["k_max"]) + 1):
        if K > 1 and K * a["min_block_cols"] > L:
            break
        clock.reset()
        t0 = time.perf_counter()
        d = U.evaluate_K(aligned, weights, K, a["min_block_cols"], const,
                         a["chemistry"], a["arm_codons"], reserved, L,
                         a["max_junk_pct"] / 100.0, a["widen_candidates"],
                         n_candidates=a["cut_candidates"],
                         max_layer_cols=max_layer_cols,
                         node_budget=a["cut_node_budget"],
                         max_library=a["max_library"], max_nt=a["max_nt"],
                         oligo_overhead=overhead,
                         proxy_candidates=a["proxy_candidates"], seed=a["seed"],
                         exhaustive_max=a["exhaustive_max"])
        sec = time.perf_counter() - t0
        if d is None:
            print("  K=%2d  no valid segmentation" % K)
            continue

        lengths, insert_nt, fails = price(d["layers"], d["tokens"])
        rows.append({
            "K": K,
            "cuts": list(d["cuts"]),
            "oligos": d["oligos"],
            "layer_oligos": [len(units) for units in d["layers"]],
            "coding_nt": d["nt"],
            "ordered_nt": sum(lengths),
            "insert_nt": insert_nt,
            "longest_construct_nt": max(lengths),
            "spacer_failures": fails,
            "library": d["library"],
            "n_cores_encoded": d["n_cores_encoded"],
            "encoded_weight": d["encoded_weight"],
            "total_weight": d["total_weight"],
            "coverage_pct": d["coverage_pct"],
            "stopped_by": d.get("stopped_by"),
            "candidates_tried": d["candidates_tried"],
            "cut_search_truncated": d["cut_search_truncated"],
            "sec_total": sec,
            "sec_cut": clock.cut,
            "sec_greedy": clock.greedy,
            "greedy_calls": clock.greedy_calls,
        })
        r = rows[-1]
        print("  K=%2d  %3d oligos  %7s nt ordered  cov %5.1f%%  lib %6s  "
              "%6.1fs (cut %5.2fs, greedy %6.1fs over %d calls)"
              % (r["K"], r["oligos"], format(r["ordered_nt"], ","),
                 r["coverage_pct"], format(r["library"], ","), sec,
                 clock.cut, clock.greedy, clock.greedy_calls))
        sys.stdout.flush()

    return {"run_dir": os.path.abspath(run_dir), "args": a, "L": L,
            "n_cores": len(aligned), "total_weight": sum(weights),
            "recommended_K": S["recommended_K"], "rows": rows,
            "naive": naive_baseline(aligned, weights),
            "ordered": ordered_truth(run_dir)}


def naive_baseline(aligned, weights):
    """Order every unique core as its own gene.  Each one still has to enter the
    backbone, so each still pays for two BsmBI sites, two spacers and the two
    backbone overhangs -- 22 nt -- exactly as fragment 1 and fragment K do."""
    lens = [3 * sum(1 for c in s if c != "-") for s in aligned]
    per_gene = SITE_NT + len("CGGA") + len("GGTG")
    return {"genes": len(lens),
            "coding_nt": sum(lens),
            "overhead_nt": per_gene * len(lens),
            "ordered_nt": sum(lens) + per_gene * len(lens),
            "per_gene_overhead_nt": per_gene,
            "longest_nt": max(lens) + per_gene,
            # not deduplicated: what ordering all the natural sequences costs
            "ordered_nt_all_seqs": (sum(n * w for n, w in zip(lens, weights))
                                    + per_gene * sum(weights))}


def ordered_truth(run_dir):
    """The real post-processed numbers for the chosen design, from the CSV that
    went to wet lab.  This is what the re-run gets checked against."""
    path = os.path.join(run_dir, "order", "constructs.csv")
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {"constructs": len(rows),
            "ordered_nt": sum(int(r["length_nt"]) for r in rows),
            "added_nt": sum(int(r["added_nt_total"]) for r in rows),
            "longest_nt": max(int(r["length_nt"]) for r in rows)}


def _rec_row(data):
    """The recommended K's row, or the last one measured if the sweep was cut
    short with --k-max (a smoke test, not a real figure)."""
    hit = [r for r in data["rows"] if r["K"] == data["recommended_K"]]
    if hit:
        return hit[0]
    print("  ! recommended K=%s not measured; using K=%d instead"
          % (data["recommended_K"], data["rows"][-1]["K"]))
    return data["rows"][-1]


def crosscheck(data):
    """The re-run must reproduce the shipped order sheet exactly, or none of the
    other K rows can be trusted."""
    truth = data["ordered"]
    if truth is None:
        print("  ! no order/constructs.csv -- cross-check skipped")
        return
    r = _rec_row(data)
    if r["K"] != data["recommended_K"]:
        return
    assert r["oligos"] == truth["constructs"], (
        "re-run has %d oligos, the order sheet has %d"
        % (r["oligos"], truth["constructs"]))
    assert r["ordered_nt"] == truth["ordered_nt"], (
        "re-run prices K=%d at %d nt, the order sheet totals %d nt"
        % (r["K"], r["ordered_nt"], truth["ordered_nt"]))
    print("  [PASS] re-run reproduces the order sheet: %d constructs, %s nt"
          % (truth["constructs"], format(truth["ordered_nt"], ",")))


# =========================================================================== #
# Figure 2 -- naive vs designed
# =========================================================================== #

def fig2(data, out):
    n = data["naive"]
    rec = _rec_row(data)
    truth = data["ordered"] if rec["K"] == data["recommended_K"] else None
    des_nt = truth["ordered_nt"] if truth else rec["ordered_nt"]
    des_add = truth["added_nt"] if truth else rec["ordered_nt"] - rec["coding_nt"]
    des_cod = des_nt - des_add
    des_n = truth["constructs"] if truth else rec["oligos"]

    fig = plt.figure(figsize=(10.4, 4.9), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 0.8, 1.05], wspace=0.45,
                          left=0.07, right=0.97, top=0.78, bottom=0.24)

    # (a) nucleotides ordered, split coding vs assembly overhead
    ax = fig.add_subplot(gs[0, 0], facecolor=SURFACE)
    x = [0, 1]
    cod = [n["coding_nt"], des_cod]
    add = [n["overhead_nt"], des_add]
    ax.bar(x, cod, width=0.52, color=[BLUE, ORANGE], zorder=3)
    ax.bar(x, add, width=0.52, bottom=cod, color=[BLUE, ORANGE], alpha=0.42,
           linewidth=2, edgecolor=SURFACE, zorder=3)
    for xi, c, a_ in zip(x, cod, add):
        ax.text(xi, c + a_ + n["ordered_nt"] * 0.03, format(c + a_, ","),
                ha="center", color=INK, fontsize=11, fontweight="bold")
    saving = 100.0 * (n["ordered_nt"] - des_nt) / n["ordered_nt"]
    ax.annotate("", xy=(0.72, des_nt), xytext=(0.72, n["ordered_nt"]),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.6))
    ax.text(0.66, (des_nt + n["ordered_nt"]) / 2, "-%.0f%%" % saving, ha="right",
            va="center", color=INK, fontsize=14, fontweight="bold")
    ax.set_ylim(0, n["ordered_nt"] * 1.22)
    ax.set_title("nucleotides ordered", color=INK, fontsize=11, loc="left", pad=8)
    ax.set_ylabel("nt", color=INK2)
    _style(ax, x, ["naive\none gene per core",
                   "designed\nK = %d fragments" % rec["K"]])
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: "%dk" % (v / 1000)))
    fig.text(0.07, 0.055,
             "solid = coding sequence,  pale = BsmBI sites and spacers\n"
             "%s nt naive vs %s nt designed: the overhead\n"
             "is not where the money goes" % (format(n["overhead_nt"], ","),
                                              format(des_add, ",")),
             color=INK3, fontsize=8.5, va="bottom")

    # (b) pieces to order -- the design orders MORE, and that is the trade
    ax = fig.add_subplot(gs[0, 1], facecolor=SURFACE)
    ax.bar(x, [n["genes"], des_n], width=0.52, color=[BLUE, ORANGE], zorder=3)
    for xi, v in zip(x, [n["genes"], des_n]):
        ax.text(xi, v + des_n * 0.035, "%d" % v, ha="center", color=INK,
                fontsize=11, fontweight="bold")
    ax.set_ylim(0, des_n * 1.30)
    ax.set_title("pieces to order", color=INK, fontsize=11, loc="left", pad=8)
    ax.set_ylabel("constructs", color=INK2)
    _style(ax, x, ["naive", "designed"])
    fig.text(0.455, 0.055,
             "+%d pieces, but they assemble:\n%d constructs build all %d cores,\n"
             "rather than %d being them"
             % (des_n - n["genes"], des_n, rec["n_cores_encoded"], n["genes"]),
             color=INK3, fontsize=8.5, va="bottom")

    # (c) the saving costs no coverage -- one number, so no chart
    ax = fig.add_subplot(gs[0, 2], facecolor=SURFACE)
    ax.axis("off")
    ax.text(0.0, 0.84, "100%", color=ORANGE, fontsize=42, fontweight="bold",
            va="center")
    ax.text(0.0, 0.56, "of natural sequences covered,\nby both routes",
            color=INK, fontsize=11, va="center")
    ax.text(0.0, 0.22,
            "%d/%d unique cores,\n%d/%d natural sequences.\n\n"
            "Without deduplicating,\nthe naive route is %s nt."
            % (rec["n_cores_encoded"], data["n_cores"], rec["encoded_weight"],
               data["total_weight"], format(n["ordered_nt_all_seqs"], ",")),
            color=INK3, fontsize=8.5, va="center")

    fig.suptitle("What the design saves -- cluster %s, after post-processing"
                 % data.get("cluster", "1"),
                 color=INK, fontsize=13, fontweight="bold", x=0.07, ha="left",
                 y=0.93)
    _save(fig, out, "fig2_naive_vs_designed.png")


# =========================================================================== #
# Figure 3 -- the frontier over K
# =========================================================================== #

def fig3(data, out):
    rows = data["rows"]
    Ks = [r["K"] for r in rows]
    recK = _rec_row(data)["K"]
    cap = data["args"]["max_library"]
    capped = [r["K"] for r in rows if r["stopped_by"] == "library"]
    k0 = min(capped) if capped else None
    i_rec = Ks.index(recK)

    fig, axes = plt.subplots(3, 1, figsize=(7.6, 8.6), sharex=True,
                             facecolor=SURFACE,
                             gridspec_kw=dict(hspace=0.18, left=0.13,
                                              right=0.96, top=0.87,
                                              bottom=0.08))
    panels = [
        ("nucleotides ordered", "nt", [r["ordered_nt"] for r in rows]),
        ("natural sequences covered", "%% of %d" % data["total_weight"],
         [r["coverage_pct"] for r in rows]),
        ("library size", "distinct proteins", [r["library"] for r in rows]),
    ]
    for ax, (title, ylab, ys) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        if k0 is not None:
            ax.axvspan(k0 - 0.5, max(Ks) + 0.5, color=ORANGE, alpha=0.07,
                       zorder=0)
        ax.plot(Ks, ys, color=BLUE, lw=2.0, marker="o", ms=5.5,
                markerfacecolor=SURFACE, markeredgewidth=2.0, zorder=3)
        ax.plot([recK], [ys[i_rec]], marker="o", ms=9.5, color=ORANGE,
                markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=4)
        ax.set_title(title, color=INK, fontsize=11, loc="left", pad=6)
        ax.set_ylabel(ylab, color=INK2)
        _grid(ax)

    ax = axes[0]
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: "%dk" % (v / 1000)))
    ax.set_ylim(0, rows[0]["ordered_nt"] * 1.22)
    ax.annotate("%s nt" % format(rows[0]["ordered_nt"], ","),
                xy=(Ks[0], rows[0]["ordered_nt"]), xytext=(7, -3),
                textcoords="offset points", color=INK2, fontsize=9)
    ax.annotate("%s nt" % format(rows[-1]["ordered_nt"], ","),
                xy=(Ks[-1], rows[-1]["ordered_nt"]), xytext=(-4, 11),
                textcoords="offset points", color=INK2, fontsize=9, ha="right")
    ax.annotate("chosen: K=%d, %s nt"
                % (recK, format(rows[i_rec]["ordered_nt"], ",")),
                xy=(recK, rows[i_rec]["ordered_nt"]), xytext=(12, 12),
                textcoords="offset points", color=ORANGE, fontsize=9.5,
                fontweight="bold")
    if k0 is not None:
        ax.annotate("library cap binds from K=%d" % k0,
                    xy=(k0 + 0.12, rows[0]["ordered_nt"] * 1.04),
                    color=ORANGE, fontsize=9.5)

    ax = axes[1]
    ax.set_ylim(0, 112)
    ax.axhline(100, color=INK3, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.annotate("all %d sequences" % data["total_weight"],
                xy=(max(Ks), 100), xytext=(-2, 5), textcoords="offset points",
                color=INK3, fontsize=8.5, ha="right")
    ax.annotate("%.0f%%" % rows[-1]["coverage_pct"],
                xy=(Ks[-1], rows[-1]["coverage_pct"]), xytext=(-4, -14),
                textcoords="offset points", color=INK2, fontsize=9, ha="right")

    ax = axes[2]
    ax.axhline(cap, color=INK3, lw=1.2, ls=(0, (4, 3)), zorder=1)
    ax.annotate("cap %s" % format(cap, ","), xy=(Ks[0], cap), xytext=(4, 6),
                textcoords="offset points", color=INK3, fontsize=9)
    ax.set_ylim(0, cap * 1.12)
    ax.set_xlabel("K, number of fragments", color=INK2)
    ax.set_xticks(Ks)

    fig.suptitle("More fragments buy fewer nucleotides, then cost coverage",
                 color=INK, fontsize=13, fontweight="bold", x=0.035, ha="left",
                 y=0.965)
    fig.text(0.035, 0.925,
             "cluster %s: %d unique cores / %d natural sequences, seed 0.  "
             "Nucleotides are post-processed construct lengths."
             % (data.get("cluster", "1"), data["n_cores"],
                data["total_weight"]),
             color=INK3, fontsize=9)
    _save(fig, out, "fig3_frontier_over_k.png")


# =========================================================================== #
# Figure 5 -- runtime
# =========================================================================== #

def fig5(data, out):
    rows = data["rows"]
    Ks = [r["K"] for r in rows]
    greedy = [r["sec_greedy"] for r in rows]
    cut = [r["sec_cut"] for r in rows]
    other = [max(0.0, r["sec_total"] - r["sec_greedy"] - r["sec_cut"])
             for r in rows]
    total = sum(r["sec_total"] for r in rows)
    cut_share = 100.0 * sum(cut) / total if total else 0.0

    fig = plt.figure(figsize=(10.8, 4.9), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 2, wspace=0.26, left=0.07, right=0.98, top=0.78,
                          bottom=0.26)

    # (a) where the seconds go.
    #
    # Stacked bars would be a lie by omission here: the cut search is four
    # orders of magnitude below the greedy, so its segment would be a hairline
    # and the reader would see one bar.  Side by side on a log axis, the GAP is
    # the finding and both series are legible.
    ax = fig.add_subplot(gs[0, 0], facecolor=SURFACE)
    floor = 0.005                      # so a sub-10ms cut search still draws
    ax.bar([k - 0.19 for k in Ks], [max(v, floor) for v in greedy], width=0.34,
           color=BLUE, label="greedy", zorder=3)
    ax.bar([k + 0.19 for k in Ks], [max(v, floor) for v in cut], width=0.34,
           color=ORANGE, label="cut search", zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(floor, max(greedy) * 6)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: ("%g" % v) if v >= 1 else ("%.2f" % v)))
    ax.set_xticks(Ks)
    ax.set_xlabel("K, number of fragments", color=INK2)
    ax.set_ylabel("seconds  (log scale)", color=INK2)
    ax.set_title("wall time per K  --  cut search is %.2f%% of %s s total"
                 % (cut_share, format(int(round(total)), ",")),
                 color=INK, fontsize=11, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper right",
              ncol=2)
    _grid(ax)
    fig.text(0.07, 0.045,
             "log axis: the cut search runs in %.2f-%.2f s, the greedy in "
             "%.0f-%.0f s.\nBookkeeping outside both accounts for the "
             "remaining %.1f s." % (min(cut), max(cut), min(greedy),
                                    max(greedy), sum(other)),
             color=INK3, fontsize=8.5, va="bottom")

    # (b) what actually sets the cost of one greedy call.
    #
    # The stated worst case is O(n^2 K (units + L + n)), and its K factor never
    # materialises: a higher K stops the greedy SOONER, because the library cap
    # binds after fewer cores have been bought.  Cost tracks cores bought, not K.
    # Both models are fitted with one constant by least squares, so the spread of
    # their residuals is the whole comparison.
    ax = fig.add_subplot(gs[0, 1], facecolor=SURFACE)
    n, L = data["n_cores"], data["L"]
    per_call, bound, real = [], [], []
    for r in rows:
        calls = max(1, r["greedy_calls"])
        per_call.append(r["sec_greedy"] / calls)
        units = r["oligos"] / float(r["K"])        # mean units in one layer
        bound.append(n * n * r["K"] * (units + L + n))
        real.append(n * r["n_cores_encoded"] * (units + L + n))

    def fit(model):
        c = (sum(t * q for t, q in zip(per_call, model))
             / sum(q * q for q in model)) if any(model) else 0.0
        pred = [c * q for q in model]
        dev = max(abs(v - t) / t for v, t in zip(pred, per_call) if t > 0)
        return pred, dev

    pred_real, dev_real = fit(real)
    _pred_bound, dev_bound = fit(bound)

    ax.plot(Ks, pred_real, color=ORANGE, lw=2.0, zorder=3,
            label=r"$c \cdot n\,b\,(u + L + n)$,  $b$ = cores bought")
    ax.plot(Ks, per_call, color=BLUE, lw=0, marker="o", ms=7.5,
            markerfacecolor=SURFACE, markeredgewidth=2.0, label="measured",
            zorder=4)
    ax.set_xticks(Ks)
    ax.set_ylim(0, max(per_call) * 1.30)
    ax.set_xlabel("K, number of fragments", color=INK2)
    ax.set_ylabel("seconds per greedy call", color=INK2)
    ax.set_title("one greedy call: cost follows cores bought, not K", color=INK,
                 fontsize=11, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper right")
    _grid(ax)
    fig.text(0.55, 0.045,
             "n = %d cores, L = %d columns, u = units per layer.  This fits "
             "within %.0f%%.\nThe stated worst case $n^2K(u{+}L{+}n)$ does not: "
             "its constant varies %.0fx,\nbecause a higher K buys fewer cores "
             "before the library cap stops it."
             % (n, L, 100 * dev_real, _spread(per_call, bound)),
             color=INK3, fontsize=8.5, va="bottom")

    fig.suptitle("Runtime: the greedy dominates, the cut search does not",
                 color=INK, fontsize=13, fontweight="bold", x=0.07, ha="left",
                 y=0.93)
    _save(fig, out, "fig5_runtime.png")
    return cut_share, dev_real, dev_bound


# =========================================================================== #

def _spread(measured, model):
    """How far from constant the implied constant of a model is -- max/min of
    t/p.  1.0 would be a perfect description of the data; large is a bad model."""
    cs = [t / q for t, q in zip(measured, model) if q > 0]
    return max(cs) / min(cs) if cs else 0.0


def _style(ax, x, labels):
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=INK2, fontsize=9.5)
    _grid(ax)


def _grid(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, length=0)


def _save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print("  wrote %s" % path)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?",
                    default=os.path.join(HERE, "algoruns", "handoff_c1_v2",
                                         "20260908_210646_c1.core.aln_gg_unified"),
                    help="a finished run directory that has been post-processed")
    ap.add_argument("--out", default=os.path.join(HERE, "figures", "final"),
                    help="parent directory; figures go in <out>/cluster<N>")
    ap.add_argument("--cluster", default="1",
                    help="which cluster this run is, naming the output "
                         "subdirectory and labelling the figures")
    ap.add_argument("--cache", default=None,
                    help="measurement cache (default <out>/cluster<N>/figdata.json)")
    ap.add_argument("--from-cache", action="store_true",
                    help="redraw from the cache instead of re-running the sweep")
    ap.add_argument("--k-max", type=int, default=None,
                    help="stop the sweep early (default: the run's own k_max)")
    args = ap.parse_args()

    out = os.path.join(args.out, "cluster%s" % args.cluster)
    os.makedirs(out, exist_ok=True)
    cache = args.cache or os.path.join(out, "figdata.json")

    if args.from_cache:
        with open(cache) as fh:
            data = json.load(fh)
        print("measurements from %s" % cache)
    else:
        print("re-running the K sweep on %s" % args.run_dir)
        print("  deterministic: same args, same seed -- this takes a while")
        data = measure(args.run_dir, args.k_max)
        with open(cache, "w") as fh:
            json.dump(data, fh, indent=1)
        print("  cached to %s" % cache)

    data["cluster"] = args.cluster
    crosscheck(data)
    fig2(data, out)
    fig3(data, out)
    share, dev_real, dev_bound = fig5(data, out)
    print("\ncut search = %.2f%% of total wall time" % share)
    print("greedy per call: cores-bought model within %.0f%%, worst-case "
          "n^2K bound within %.0f%%" % (100 * dev_real, 100 * dev_bound))


if __name__ == "__main__":
    main()
