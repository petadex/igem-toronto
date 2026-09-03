"""make_dp_figures.py -- figures comparing the DFS and DP/CSP cut searches (issue #84).

    fig4_cutsearch.png     search cost and search work, DFS vs DP, by K
    fig5_design_ab.png     the designs those searches actually produce
    fig6_highk_frontier.png  coverage bought per nucleotide, K = 1..11

Measurements are embedded rather than re-read from run directories: the runs
live in untracked algoruns/ and temp dirs that are not committed, so a figure
script that read them would not regenerate for anyone else.  Every block below
names the exact command that produced it.

Colours are the dataviz reference palette, unchanged and in documented order:
slot 1 blue #2a78d6, slot 2 orange #eb6834.

USAGE
-----
    python make_dp_figures.py [--out figures]
"""

import argparse
import os
from math import comb

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"

N_SITES = 101          # legal Golden Gate sites on cluster 1 at min_block 20
TIME_FLOOR = 0.01      # sub-10ms timings are reported as 0.00; shown at the floor


# --------------------------------------------------------------------------- #
# Measurements.
#
# (A) cut search alone -- `python dp_cutsearch.py --k-max 6 <c1>` for the paired
#     rows, `--k-max 8 --no-compare` for the DP-only tail.  Cost is
#     log(library at full coverage); seconds are wall time for place_cuts /
#     dp_cut_search only, not the designer.
# --------------------------------------------------------------------------- #
CUT = {                     # K: (dfs_cost, dfs_sec, dfs_truncated, dp_cost, dp_sec)
    2: (6.7405, 0.00, False, 6.7405, 0.01),
    3: (7.6573, 0.08, False, 7.6573, 0.05),
    4: (9.5122, 5.56, True,  9.5122, 0.03),
    5: (11.2834, 7.67, True, 11.2834, 0.11),
    6: (13.9004, 10.35, True, 13.8441, 0.14),
    7: (None, None, None, 15.6152, 0.15),
    8: (None, None, None, 18.0588, 0.14),
}

# (B) whole designer, library cap 2,500.
#     DFS  -- the k-max 6 probe (exhaustive_max 150, no nt cap).
#     DP   -- `ab_dp_vs_dfs.py --arms dp --k-max 12 --max-library 2500
#              --max-nt 45000 --exhaustive-max 50`.
#     K=2 is NOT comparable between the arms: exhaustive_max differs, so the DFS
#     saw 331 candidates there and the DP saw 50.  Marked in the figure.
DESIGN_DFS = {              # K: (natural sequences encoded, nt ordered)
    1: (65, 40263), 2: (65, 28089), 3: (65, 26652),
    4: (57, 21303), 5: (54, 18390), 6: (49, 11931),
}
DESIGN_DP = {
    1: (65, 40263), 2: (65, 28152), 3: (65, 26706),
    4: (57, 21303), 5: (54, 18558), 6: (49, 12798),
    7: (43, 7143), 8: (38, 4746), 9: (32, 3744),
    10: (20, 2058), 11: (20, 2082),
}

# (C) the DP frontier in full, same run as (B).
#     K: (cores, natseq, library, oligos, nt_ordered, longest_nt, stopped_by)
FRONTIER = {
    1: (54, 65, 54, 50, 40263, 786, "no improving move"),
    2: (54, 65, 1150, 61, 28152, 528, "no improving move"),
    3: (54, 65, 2116, 61, 26706, 528, "no improving move"),
    4: (46, 57, 1976, 48, 21303, 570, "library"),
    5: (43, 54, 1820, 46, 18558, 537, "library"),
    6: (38, 49, 2156, 41, 12798, 447, "library"),
    7: (32, 43, 1496, 33, 7143, 384, "library"),
    8: (27, 38, 2310, 30, 4746, 327, "library"),
    9: (21, 32, 2352, 29, 3744, 234, "library"),
    10: (15, 20, 1728, 22, 2058, 180, "library"),
    11: (15, 20, 1728, 23, 2082, 114, "library"),
}
TOTAL_W = 65

DFS_BUDGET = 2_000_000
DP_EDGE_BASE, DP_EDGE_STEP = 202, 4187     # measured: edges = base + (K-2)*step


def style(ax):
    """Recessive axes and grid; the data carries the ink."""
    ax.set_facecolor(SURFACE)
    ax.grid(True, color="#e6e5e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d4d2cb")
    ax.tick_params(colors=INK2, labelsize=9, length=0)


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, name)
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", p)


def fig_cutsearch(outdir):
    """Two panels, because seconds and path counts are different scales and a
    dual axis would be a lie.  Left: what the search costs to run.  Right: why."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # ---- left: wall time -------------------------------------------------- #
    ks = sorted(CUT)
    d_k = [k for k in ks if CUT[k][1] is not None]
    d_t = [max(CUT[k][1], TIME_FLOOR) for k in d_k]
    p_k = ks
    p_t = [max(CUT[k][4], TIME_FLOOR) for k in ks]

    ax1.plot(d_k, d_t, color=ORANGE, linewidth=2, marker="o", markersize=8,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=3, label="DFS")
    ax1.plot(p_k, p_t, color=BLUE, linewidth=2, marker="o", markersize=8,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=3, label="DP + CSP")

    # hollow markers where the DFS returned a truncated (unproven) answer
    t_k = [k for k in d_k if CUT[k][2]]
    t_t = [max(CUT[k][1], TIME_FLOOR) for k in t_k]
    ax1.scatter(t_k, t_t, s=150, facecolors=SURFACE, edgecolors=ORANGE,
                linewidths=2, zorder=4)
    ax1.annotate("hollow = truncated,\nanswer not proven",
                 xy=(4, max(CUT[4][1], TIME_FLOOR)), xytext=(4.15, 0.55),
                 fontsize=8.5, color=ORANGE,
                 arrowprops=dict(arrowstyle="-", color=ORANGE, linewidth=1))
    ax1.annotate("DFS cannot\nreach K >= 7", xy=(6.55, 3.2), fontsize=8.5,
                 color=INK2, ha="center")
    ax1.axvspan(6.5, 8.4, color="#f0efe9", zorder=0)

    ax1.set_yscale("log")
    ax1.set_xlabel("fragments K", color=INK2, fontsize=9.5)
    ax1.set_ylabel("cut-search wall time (s, log)", color=INK2, fontsize=9.5)
    ax1.set_title("The search itself", color=INK, fontsize=11, loc="left", pad=10)
    ax1.set_xlim(1.6, 8.4)
    ax1.set_xticks(ks)
    style(ax1)
    leg = ax1.legend(frameon=False, fontsize=9, loc="upper left")
    for t in leg.get_texts():
        t.set_color(INK2)

    # ---- right: why ------------------------------------------------------- #
    place = [comb(N_SITES, k - 1) for k in ks]
    edges = [DP_EDGE_BASE + (k - 2) * DP_EDGE_STEP for k in ks]
    ax2.plot(ks, place, color=ORANGE, linewidth=2, marker="o", markersize=8,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=3,
             label="paths the DFS must enumerate")
    ax2.plot(ks, edges, color=BLUE, linewidth=2, marker="o", markersize=8,
             markeredgecolor=SURFACE, markeredgewidth=2, zorder=3,
             label="edges the DP relaxes")
    ax2.axhline(DFS_BUDGET, color=INK3, linewidth=1.2, linestyle=(0, (4, 3)),
                zorder=2)
    ax2.annotate("DFS node budget (2M)", xy=(2.1, DFS_BUDGET * 1.4), fontsize=8.5,
                 color=INK2)
    ax2.set_yscale("log")
    ax2.set_xlabel("fragments K", color=INK2, fontsize=9.5)
    ax2.set_ylabel("search work (log)", color=INK2, fontsize=9.5)
    ax2.set_title("Why: exponential vs linear in K", color=INK, fontsize=11,
                  loc="left", pad=10)
    ax2.set_xlim(1.6, 8.4)
    ax2.set_xticks(ks)
    style(ax2)
    leg = ax2.legend(frameon=False, fontsize=9, loc="lower right")
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.suptitle("Cut search on cluster 1 (54 cores, 262 columns, 101 legal sites)",
                 color=INK, fontsize=12.5, x=0.5, y=1.02)
    save(fig, outdir, "fig4_cutsearch.png")


def fig_design_ab(outdir):
    """The honest null: an exact search does not move the design at this cap."""
    ks = sorted(DESIGN_DFS)
    x = range(len(ks))
    w = 0.38
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    for ax, idx, lab, fmt in ((ax1, 0, "natural sequences encoded", "{:.0f}"),
                              (ax2, 1, "nucleotides ordered", "{:,.0f}")):
        dfs = [DESIGN_DFS[k][idx] for k in ks]
        dp = [DESIGN_DP[k][idx] for k in ks]
        ax.bar([i - w / 2 for i in x], dfs, w - 0.02, color=ORANGE, zorder=3,
               label="DFS")
        ax.bar([i + w / 2 for i in x], dp, w - 0.02, color=BLUE, zorder=3,
               label="DP + CSP")
        ax.set_xticks(list(x))
        ax.set_xticklabels(["K=%d" % k for k in ks])
        ax.set_ylabel(lab, color=INK2, fontsize=9.5)
        style(ax)
        if idx == 1:
            ax.yaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(lambda v, _: format(int(v), ",")))

    ax1.axhline(TOTAL_W, color=INK3, linewidth=1.2, linestyle=(0, (4, 3)), zorder=4)
    ax1.annotate("all 65 natural sequences", xy=(0.05, TOTAL_W + 1.2), fontsize=8.5,
                 color=INK2)
    ax1.set_ylim(0, TOTAL_W + 9)
    ax1.set_title("Coverage: identical at every K", color=INK, fontsize=11,
                  loc="left", pad=10)
    ax2.set_title("Cost: DP marginally worse at K=3, 5, 6", color=INK, fontsize=11,
                  loc="left", pad=10)

    # K=2 is a config difference, not a result -- say so on the figure itself.
    for ax in (ax1, ax2):
        ax.annotate("*", xy=(1, 0), xytext=(1, -0.13), textcoords=("data",
                    "axes fraction"), ha="center", fontsize=14, color=INK3)
    # The left panel has no free space -- its bars reach the top of the axis at
    # K=1..3 and the floor elsewhere -- and above the axis collides with the
    # panel title.  The right panel's bars fall away to the right, so the legend
    # goes there and serves both panels.
    leg = ax2.legend(frameon=False, fontsize=9, loc="upper right")
    for t in leg.get_texts():
        t.set_color(INK2)

    fig.suptitle("Designs produced, library cap 2,500 (cluster 1)",
                 color=INK, fontsize=12.5, x=0.5, y=1.02)
    fig.text(0.5, -0.05, "* K=2 not comparable: the DFS run used exhaustive_max "
             "150 (331 candidates), the DP run used 50",
             ha="center", fontsize=8.5, color=INK3)
    save(fig, outdir, "fig5_design_ab.png")


def fig_highk(outdir):
    """What high K actually buys: a coverage/cost frontier, one point per K."""
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ks = sorted(FRONTIER)
    nt = [FRONTIER[k][4] for k in ks]
    nat = [FRONTIER[k][1] for k in ks]
    capped = [FRONTIER[k][6] == "library" for k in ks]

    ax.plot(nt, nat, color=BLUE, linewidth=2, zorder=2)
    full = [i for i, c in enumerate(capped) if not c]
    cap = [i for i, c in enumerate(capped) if c]
    ax.scatter([nt[i] for i in full], [nat[i] for i in full], s=90, color=BLUE,
               edgecolors=SURFACE, linewidths=2, zorder=3,
               label="stopped: no improving move")
    ax.scatter([nt[i] for i in cap], [nat[i] for i in cap], s=90,
               facecolors=SURFACE, edgecolors=BLUE, linewidths=2, zorder=3,
               label="stopped: library cap bound")

    # Hand-placed label offsets.  Several points sit almost on top of each other
    # -- K=2/K=3 differ by 1,446 nt at the same coverage, and K=10/K=11 by 24 nt
    # at the same coverage -- so a uniform offset collides.
    OFF = {1: (0, 10), 2: (20, 8), 3: (-20, 8), 4: (2, 10), 5: (0, 10),
           6: (0, 10), 7: (0, 10), 8: (-4, 10), 9: (-20, 2), 10: (0, -18),
           11: (22, 4)}
    for k, x, y in zip(ks, nt, nat):
        ax.annotate("K=%d" % k, xy=(x, y), xytext=OFF[k],
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color=INK2)

    ax.axhline(TOTAL_W, color=INK3, linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
    ax.annotate("all 65 natural sequences", xy=(2600, TOTAL_W + 1.4), fontsize=8.5,
                color=INK2)
    ax.set_xlabel("nucleotides ordered", color=INK2, fontsize=9.5)
    ax.set_ylabel("natural sequences encoded", color=INK2, fontsize=9.5)
    ax.set_ylim(0, TOTAL_W + 9)
    ax.set_xlim(0, 44000)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: format(int(v), ",")))
    style(ax)
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right")
    for t in leg.get_texts():
        t.set_color(INK2)
    ax.set_title("What fragmenting buys: coverage against DNA ordered\n"
                 "cluster 1, library cap 2,500, K = 1..11",
                 color=INK, fontsize=12, loc="left", pad=12)
    save(fig, outdir, "fig6_highk_frontier.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "figures"))
    args = ap.parse_args()
    fig_cutsearch(args.out)
    fig_design_ab(args.out)
    fig_highk(args.out)


if __name__ == "__main__":
    main()
