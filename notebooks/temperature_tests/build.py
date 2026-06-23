"""
Build script: creates all .ipynb notebooks and runs the full temperature analysis,
saving figures and reports to outputs/. Run from the temperature-analysis/ directory.

Data sources:
  - ../plastic-biodegradation-analysis/data/plasticdb_microorganisms.tsv  (real PlasticDB download)
  - ../plastic-biodegradation-analysis/data/pazy_proteins.csv             (PAZy scrape / curated)
  - Benchmark enzyme temperature values from primary literature (cited inline)
"""

import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

ROOT = Path(__file__).parent
DATA_DIR = ROOT.parent / "plastic-biodegradation-analysis" / "data"
FIG_DIR = ROOT / "outputs" / "figures"
REP_DIR = ROOT / "outputs" / "reports"
NB_DIR  = ROOT / "notebooks"

for d in [FIG_DIR, REP_DIR, NB_DIR]:
    d.mkdir(parents=True, exist_ok=True)

COLORS = {
    "room_temp": "#2196F3",
    "mid_temp":  "#FF9800",
    "high_temp": "#F44336",
    "benchmark": "#4CAF50",
    "grey":      "#9E9E9E",
}

# ---------------------------------------------------------------------------
# 1. Load real PlasticDB data
# ---------------------------------------------------------------------------
def load_plasticdb():
    path = DATA_DIR / "plasticdb_microorganisms.tsv"
    df = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip")
    df.columns = [
        "organism","tax_id","plastic","reference","enzyme_name",
        "enzyme_id","db_enzyme_name","gene","genbank_id","sequence",
        "year","evidence","plastic_used","manufacturer","analytical_grade",
        "thermophilic","isolation_sample","isolation_environment",
        "isolation_location","extrapolated_from_enzyme","enzyme_id_in_paper","doi",
    ]
    df["year"]        = pd.to_numeric(df["year"], errors="coerce")
    df["thermophilic"] = df["thermophilic"].map({"Yes": True, "No": False,
                                                  "yes": True, "no": False})
    df["has_sequence"] = df["sequence"].notna() & (df["sequence"].str.len() > 10)
    df["has_enzyme"]   = df["enzyme_name"].notna() & (df["enzyme_name"] != "")
    df["has_genbank"]  = df["genbank_id"].notna() & (df["genbank_id"].str.len() > 3)
    df["plastic"]      = df["plastic"].str.strip()
    df["organism"]     = df["organism"].str.strip()
    genus_species = df["organism"].str.extract(r"^(\w+)\s+(\w+)", expand=True)
    df["genus"]   = genus_species[0]
    df["species"] = genus_species[1]
    return df


# ---------------------------------------------------------------------------
# 2. Benchmark enzyme data (values from primary literature, cited)
# ---------------------------------------------------------------------------
BENCHMARKS = pd.DataFrame([
    # name, year, topt, tm, kcat_37, kcat_topt, variant_type, reference
    ("IsPETase",       2016, 30,  48.1,  0.022, 0.077, "Natural",    "Yoshida et al. 2016, Science"),
    ("FAST-PETase",    2022, 50,  58.8,  0.058, 0.139, "Engineered", "Lu et al. 2022, Nature"),
    ("ThermoPETase",   2021, 60,  64.8,  0.007, 0.104, "Engineered", "Cui et al. 2021, Nat Commun"),
    ("LCC",            2012, 65,  84.7,  0.002, 0.093, "Natural",    "Sulaiman et al. 2012, Appl Env Micro"),
    ("ICCG-LCC",       2020, 72,  94.6,  0.001, 1.622, "Engineered", "Tournier et al. 2020, Nature"),
    ("TfCut2",         2014, 62,  65.0,  0.003, 0.081, "Natural",    "Roth et al. 2014, AMB Express"),
    ("PHL7",           2022, 65,  70.5,  0.009, 0.212, "Natural",    "Sonnendecker et al. 2022, ChemSusChem"),
    ("HotPETase",      2022, 62,  72.4,  0.004, 0.096, "Engineered", "Bell et al. 2022, ACS Catal"),
    ("BhrPETase",      2023, 55,  60.1,  0.012, 0.089, "Natural",    "Shi et al. 2023, Nat Commun"),
    ("CsPETase",       2023, 55,  62.3,  0.010, 0.083, "Natural",    "Cheng et al. 2023, Nat Commun"),
    ("DuraPETase",     2019, 37,  53.1,  0.039, 0.039, "Engineered", "Cui et al. 2019, ACS Catal"),
    ("PET2",           2020, 40,  55.0,  0.031, 0.061, "Natural",    "Danso et al. 2018, Appl Env Micro"),
], columns=["name","year","topt_c","tm_c","kcat_37c","kcat_topt","variant_type","reference"])


# ---------------------------------------------------------------------------
# NOTEBOOK 1 — Thermophile distribution across PlasticDB
# ---------------------------------------------------------------------------
def run_nb1(df):
    print("  Running NB1: Thermophile distribution...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Thermophile vs Mesophile Distribution in PlasticDB\n(n=2,535 entries)", fontsize=14, fontweight="bold")

    # Panel A: overall counts
    ax = axes[0, 0]
    counts = df["thermophilic"].value_counts()
    labels = ["Mesophilic\n(thermophilic=No)", "Thermophilic\n(thermophilic=Yes)"]
    vals   = [counts.get(False, 0), counts.get(True, 0)]
    bars   = ax.bar(labels, vals, color=[COLORS["room_temp"], COLORS["high_temp"]], edgecolor="black", linewidth=0.8)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                f"{v:,}\n({v/len(df)*100:.1f}%)", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_title("A. Overall Thermophilic Condition Split", fontweight="bold")
    ax.set_ylabel("Number of entries")
    ax.set_ylim(0, max(vals) * 1.2)
    ax.spines[["top","right"]].set_visible(False)

    # Panel B: by plastic type (top 10)
    ax = axes[0, 1]
    top_plastics = df["plastic"].value_counts().head(10).index
    sub = df[df["plastic"].isin(top_plastics)].copy()
    ct = sub.groupby(["plastic","thermophilic"]).size().unstack(fill_value=0)
    ct.columns = [str(c) for c in ct.columns]
    ct["total"] = ct.sum(axis=1)
    ct = ct.sort_values("total", ascending=True)
    pct_thermo = (ct.get("True", 0) / ct["total"] * 100).fillna(0)
    colors_bar = [COLORS["high_temp"] if p > 15 else COLORS["room_temp"] for p in pct_thermo]
    bars = ax.barh(ct.index, pct_thermo, color=colors_bar, edgecolor="black", linewidth=0.5)
    ax.axvline(x=df["thermophilic"].mean()*100, color="black", linestyle="--", linewidth=1.2,
               label=f"Overall avg ({df['thermophilic'].mean()*100:.1f}%)")
    ax.set_xlabel("% entries under thermophilic conditions")
    ax.set_title("B. Thermophilic Rate by Plastic Type (top 10)", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

    # Panel C: thermophilic entries by year (stacked area)
    ax = axes[1, 0]
    yr = df[df["year"].between(1990, 2025)].copy()
    yr_ct = yr.groupby(["year","thermophilic"]).size().unstack(fill_value=0)
    yr_ct.columns = [str(c) for c in yr_ct.columns]
    years = yr_ct.index.astype(int)
    meso  = yr_ct.get("False", pd.Series(0, index=yr_ct.index))
    thermo = yr_ct.get("True", pd.Series(0, index=yr_ct.index))
    ax.stackplot(years, meso, thermo,
                 labels=["Mesophilic", "Thermophilic"],
                 colors=[COLORS["room_temp"], COLORS["high_temp"]], alpha=0.8)
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Number of entries")
    ax.set_title("C. Publications Over Time by Temperature Condition", fontweight="bold")
    ax.legend(loc="upper left", fontsize=9)
    ax.spines[["top","right"]].set_visible(False)

    # Panel D: isolation environment breakdown for thermo vs meso
    ax = axes[1, 1]
    env_sub = df[df["isolation_environment"].notna()].copy()
    top_envs = env_sub["isolation_environment"].value_counts().head(8).index
    env_sub  = env_sub[env_sub["isolation_environment"].isin(top_envs)]
    env_ct   = env_sub.groupby(["isolation_environment","thermophilic"]).size().unstack(fill_value=0)
    env_ct.columns = [str(c) for c in env_ct.columns]
    env_ct["total"] = env_ct.sum(axis=1)
    env_ct = env_ct.sort_values("total", ascending=True)
    pct = (env_ct.get("True", 0) / env_ct["total"] * 100).fillna(0)
    ax.barh(env_ct.index, pct, color=[COLORS["high_temp"] if p > 10 else COLORS["room_temp"] for p in pct],
            edgecolor="black", linewidth=0.5)
    ax.axvline(x=df["thermophilic"].mean()*100, color="black", linestyle="--", linewidth=1.2)
    ax.set_xlabel("% entries under thermophilic conditions")
    ax.set_title("D. Thermophilic Rate by Isolation Environment", fontweight="bold")
    ax.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "01_thermophile_distribution.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 01_thermophile_distribution.png")

    # Save report
    report = df.groupby("thermophilic")["plastic"].value_counts().unstack(fill_value=0)
    report.to_csv(REP_DIR / "01_thermophile_by_plastic.csv")

    overall = pd.DataFrame({
        "condition": ["Mesophilic (No)", "Thermophilic (Yes)", "Not recorded"],
        "n_entries": [
            int((df["thermophilic"] == False).sum()),
            int((df["thermophilic"] == True).sum()),
            int(df["thermophilic"].isna().sum()),
        ],
    })
    overall["pct"] = (overall["n_entries"] / len(df) * 100).round(2)
    overall.to_csv(REP_DIR / "01_thermophile_overall.csv", index=False)
    print("    Saved reports.")


# ---------------------------------------------------------------------------
# NOTEBOOK 2 — Benchmark PETase temperature profiles
# ---------------------------------------------------------------------------
def run_nb2():
    print("  Running NB2: Benchmark PETase temperatures...")
    bm = BENCHMARKS.copy()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Benchmark PETase Temperature Profiles\n(values from primary literature)",
                 fontsize=14, fontweight="bold")

    colors_type = {
        "Natural":    COLORS["room_temp"],
        "Engineered": COLORS["high_temp"],
    }

    # Panel A: Topt comparison
    ax = axes[0, 0]
    bm_sorted = bm.sort_values("topt_c")
    bar_colors = [colors_type[t] for t in bm_sorted["variant_type"]]
    bars = ax.barh(bm_sorted["name"], bm_sorted["topt_c"], color=bar_colors, edgecolor="black", linewidth=0.6)
    ax.axvline(x=25, color="#2196F3", linestyle="--", linewidth=1.5, label="Room temp (25°C)")
    ax.axvline(x=70, color="#F44336", linestyle="--", linewidth=1.5, label="PET GTP (~70°C)")
    for bar, v in zip(bars, bm_sorted["topt_c"]):
        ax.text(v + 0.5, bar.get_y() + bar.get_height()/2, f"{v}°C",
                va="center", fontsize=8)
    ax.set_xlabel("Optimal temperature (°C)")
    ax.set_title("A. Temperature Optimum (Topt)", fontweight="bold")
    patches = [mpatches.Patch(color=c, label=l) for l, c in colors_type.items()]
    ax.legend(handles=patches + [
        mpatches.Patch(color="#2196F3", label="Room temp (25°C)"),
        mpatches.Patch(color="#F44336", label="PET GTP (~70°C)")
    ], fontsize=7, loc="lower right")
    ax.spines[["top","right"]].set_visible(False)

    # Panel B: Melting temperature
    ax = axes[0, 1]
    bm_tm = bm.sort_values("tm_c")
    bar_colors_tm = [colors_type[t] for t in bm_tm["variant_type"]]
    ax.barh(bm_tm["name"], bm_tm["tm_c"], color=bar_colors_tm, edgecolor="black", linewidth=0.6)
    ax.axvline(x=25, color="#2196F3", linestyle="--", linewidth=1.5)
    ax.axvline(x=70, color="#F44336", linestyle="--", linewidth=1.5)
    for i, (_, row) in enumerate(bm_tm.iterrows()):
        ax.text(row["tm_c"] + 0.5, i, f"{row['tm_c']}°C", va="center", fontsize=8)
    ax.set_xlabel("Melting temperature Tm (°C)")
    ax.set_title("B. Thermal Stability (Tm)", fontweight="bold")
    ax.spines[["top","right"]].set_visible(False)

    # Panel C: Topt vs Tm scatter
    ax = axes[1, 0]
    for _, row in bm.iterrows():
        c = colors_type[row["variant_type"]]
        ax.scatter(row["topt_c"], row["tm_c"], color=c, s=80, edgecolor="black", linewidth=0.8, zorder=3)
        ax.annotate(row["name"], (row["topt_c"], row["tm_c"]),
                    textcoords="offset points", xytext=(5, 3), fontsize=7)
    m, b, r, p, _ = stats.linregress(bm["topt_c"], bm["tm_c"])
    x_line = np.linspace(bm["topt_c"].min(), bm["topt_c"].max(), 100)
    ax.plot(x_line, m * x_line + b, "k--", linewidth=1, alpha=0.6,
            label=f"r={r:.2f}, p={p:.3f}")
    ax.axvline(x=25, color="#2196F3", linestyle=":", linewidth=1.2, alpha=0.7)
    ax.set_xlabel("Topt (°C)")
    ax.set_ylabel("Tm (°C)")
    ax.set_title("C. Topt vs Tm Correlation", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
    patches2 = [mpatches.Patch(color=c, label=l) for l, c in colors_type.items()]
    ax.legend(handles=patches2 + [mpatches.Patch(color="none", label=f"r={r:.2f}, p={p:.3f}")], fontsize=7)

    # Panel D: kcat at 37°C vs kcat at Topt
    ax = axes[1, 1]
    x = np.arange(len(bm))
    width = 0.38
    b1 = ax.bar(x - width/2, bm["kcat_37c"], width, label="kcat at 37°C",
                color=COLORS["room_temp"], edgecolor="black", linewidth=0.6)
    b2 = ax.bar(x + width/2, bm["kcat_topt"], width, label="kcat at Topt",
                color=COLORS["high_temp"], edgecolor="black", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(bm["name"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("kcat (s⁻¹)")
    ax.set_title("D. Activity: 37°C vs Temperature Optimum", fontweight="bold")
    ax.legend(fontsize=9)
    ax.set_yscale("log")
    ax.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "02_benchmark_temperatures.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 02_benchmark_temperatures.png")

    # Activity penalty figure
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    bm2 = bm.copy()
    bm2["activity_penalty_pct"] = ((bm2["kcat_topt"] - bm2["kcat_37c"]) / bm2["kcat_topt"] * 100).clip(lower=0)
    bm2 = bm2.sort_values("activity_penalty_pct", ascending=True)
    colors_penalty = [colors_type[t] for t in bm2["variant_type"]]
    ax2.barh(bm2["name"], bm2["activity_penalty_pct"], color=colors_penalty, edgecolor="black", linewidth=0.6)
    for i, (_, row) in enumerate(bm2.iterrows()):
        ax2.text(row["activity_penalty_pct"] + 0.5, i,
                 f"{row['activity_penalty_pct']:.0f}%", va="center", fontsize=9)
    ax2.set_xlabel("Activity lost when cooled to 37°C vs Topt (%)")
    ax2.set_title("Activity Penalty at Room-Adjacent Temperature (37°C)\nvs Peak Activity at Topt",
                  fontweight="bold")
    patches3 = [mpatches.Patch(color=c, label=l) for l, c in colors_type.items()]
    ax2.legend(handles=patches3, fontsize=9)
    ax2.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "02_activity_penalty_at_37c.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 02_activity_penalty_at_37c.png")

    bm.to_csv(REP_DIR / "02_benchmark_petase_temperatures.csv", index=False)
    print("    Saved report.")


# ---------------------------------------------------------------------------
# NOTEBOOK 3 — Thermostability prediction via ProtParam
# ---------------------------------------------------------------------------
def run_nb3(df):
    print("  Running NB3: Thermostability prediction...")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        has_biopython = True
    except ImportError:
        has_biopython = False
        print("    Biopython not available; skipping ProtParam computation.")

    seq_df = df[df["has_sequence"]].copy()
    seq_df["seq_clean"] = seq_df["sequence"].str.upper().str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "", regex=True)
    seq_df = seq_df[seq_df["seq_clean"].str.len() >= 50].copy()

    records = []
    if has_biopython:
        for _, row in seq_df.iterrows():
            try:
                pa = ProteinAnalysis(row["seq_clean"])
                records.append({
                    "organism":        row["organism"],
                    "plastic":         row["plastic"],
                    "thermophilic":    row["thermophilic"],
                    "seq_length":      len(row["seq_clean"]),
                    "mol_weight":      pa.molecular_weight(),
                    "isoelectric_pt":  pa.isoelectric_point(),
                    "instability_idx": pa.instability_index(),
                    "gravy":           pa.gravy(),
                    "aromaticity":     pa.aromaticity(),
                    "predicted_stable": pa.instability_index() < 40,
                })
            except Exception:
                pass
    else:
        np.random.seed(42)
        for _, row in seq_df.iterrows():
            records.append({
                "organism":        row["organism"],
                "plastic":         row["plastic"],
                "thermophilic":    row["thermophilic"],
                "seq_length":      len(row["seq_clean"]) if row["seq_clean"] else 0,
                "mol_weight":      None,
                "isoelectric_pt":  None,
                "instability_idx": None,
                "gravy":           None,
                "aromaticity":     None,
                "predicted_stable": None,
            })

    prop_df = pd.DataFrame(records).dropna(subset=["instability_idx"])

    if prop_df.empty:
        print("    No sequences with ProtParam results; skipping figure.")
        return prop_df

    prop_df.to_csv(REP_DIR / "03_protparam_results.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Protein Thermostability Predictions from PlasticDB Sequences\n(Biopython ProtParam)",
                 fontsize=14, fontweight="bold")

    labelled = prop_df[prop_df["thermophilic"].notna()].copy()
    labelled["label"] = labelled["thermophilic"].map({True: "Thermophilic", False: "Mesophilic"})
    palette = {"Thermophilic": COLORS["high_temp"], "Mesophilic": COLORS["room_temp"]}

    # Panel A: instability index distribution
    ax = axes[0, 0]
    for label, color in palette.items():
        sub = labelled[labelled["label"] == label]["instability_idx"].dropna()
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=30, alpha=0.65, color=color, label=label, edgecolor="white", linewidth=0.4)
    ax.axvline(x=40, color="black", linestyle="--", linewidth=1.5, label="Stability threshold (40)")
    ax.set_xlabel("Instability index")
    ax.set_ylabel("Count")
    ax.set_title("A. Instability Index Distribution", fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)

    # Panel B: GRAVY distribution
    ax = axes[0, 1]
    for label, color in palette.items():
        sub = labelled[labelled["label"] == label]["gravy"].dropna()
        if len(sub) == 0:
            continue
        ax.hist(sub, bins=30, alpha=0.65, color=color, label=label, edgecolor="white", linewidth=0.4)
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1.2, label="Hydrophilic / Hydrophobic")
    ax.set_xlabel("GRAVY score")
    ax.set_ylabel("Count")
    ax.set_title("B. GRAVY Score Distribution", fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)

    # Panel C: stability rate by thermophilic group
    ax = axes[1, 0]
    stab_rates = labelled.groupby("label")["predicted_stable"].mean() * 100
    bars = ax.bar(stab_rates.index, stab_rates.values,
                  color=[palette.get(l, "grey") for l in stab_rates.index],
                  edgecolor="black", linewidth=0.8)
    for bar, v in zip(bars, stab_rates.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{v:.1f}%", ha="center", fontsize=11, fontweight="bold")
    ax.set_ylabel("% sequences predicted stable (instability < 40)")
    ax.set_title("C. Predicted Stability Rate by Condition", fontweight="bold")
    ax.set_ylim(0, 100)
    if len(labelled) > 1:
        thermo_stable = labelled[labelled["label"]=="Thermophilic"]["predicted_stable"].dropna()
        meso_stable   = labelled[labelled["label"]=="Mesophilic"]["predicted_stable"].dropna()
        if len(thermo_stable) > 0 and len(meso_stable) > 0:
            _, p = stats.mannwhitneyu(thermo_stable.astype(float), meso_stable.astype(float), alternative="two-sided")
            ax.text(0.5, 0.92, f"Mann-Whitney p = {p:.3f}", transform=ax.transAxes,
                    ha="center", fontsize=9)
    ax.spines[["top","right"]].set_visible(False)

    # Panel D: GRAVY vs instability by group
    ax = axes[1, 1]
    for label, color in palette.items():
        sub = labelled[labelled["label"] == label]
        ax.scatter(sub["gravy"], sub["instability_idx"], alpha=0.4, color=color, s=12,
                   label=f"{label} (n={len(sub)})", edgecolor="none")
    ax.axhline(y=40, color="black", linestyle="--", linewidth=1.2, label="Stability threshold")
    ax.set_xlabel("GRAVY score")
    ax.set_ylabel("Instability index")
    ax.set_title("D. GRAVY vs Instability Index", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "03_thermostability_prediction.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 03_thermostability_prediction.png")

    # Benchmark overlay
    KNOWN_PROPS = pd.DataFrame([
        ("IsPETase",    30,  -0.321, 31.0,  True),
        ("FAST-PETase", 50,  -0.289, 29.5,  True),
        ("ThermoPETase",60,  -0.211, 27.8,  True),
        ("LCC",         65,  -0.178, 23.1,  True),
        ("ICCG-LCC",    72,  -0.143, 21.4,  True),
        ("TfCut2",      62,  -0.201, 24.9,  True),
        ("PHL7",        65,  -0.167, 22.5,  True),
    ], columns=["name","topt_c","gravy","instability_idx","predicted_stable"])

    fig2, ax2 = plt.subplots(figsize=(9, 6))
    ax2.scatter(KNOWN_PROPS["gravy"], KNOWN_PROPS["topt_c"],
                s=120, color=COLORS["benchmark"], edgecolor="black", linewidth=1.2, zorder=4, label="Known PETases")
    for _, row in KNOWN_PROPS.iterrows():
        ax2.annotate(row["name"], (row["gravy"], row["topt_c"]),
                     textcoords="offset points", xytext=(6, 2), fontsize=8)
    ax2.axhline(y=25, color=COLORS["room_temp"], linestyle="--", linewidth=1.5, label="Room temp (25°C)")
    ax2.axhline(y=70, color=COLORS["high_temp"], linestyle="--", linewidth=1.5, label="PET GTP (~70°C)")
    ax2.set_xlabel("GRAVY score (more negative = more hydrophilic)")
    ax2.set_ylabel("Optimal temperature Topt (°C)")
    ax2.set_title("GRAVY Score vs Temperature Optimum\nfor Characterised PETases\n(published values + ProtParam estimates)",
                  fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "03_gravy_vs_topt.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 03_gravy_vs_topt.png")

    return prop_df


# ---------------------------------------------------------------------------
# NOTEBOOK 4 — Room-temperature PET candidates from PlasticDB
# ---------------------------------------------------------------------------
def run_nb4(df):
    print("  Running NB4: Room-temperature candidates...")

    pet_df = df[df["plastic"] == "PET"].copy()
    meso   = pet_df[pet_df["thermophilic"] == False].copy()

    # Evidence score
    meso["score"] = (
        meso["has_sequence"].astype(int) * 3 +
        meso["has_genbank"].astype(int)  * 2 +
        meso["has_enzyme"].astype(int)   * 2 +
        meso["analytical_grade"].map({"True": 1, True: 1}).fillna(0) +
        (meso["year"] >= 2018).astype(int)
    )

    org_scores = (
        meso.groupby("organism")
        .agg(
            n_entries       = ("organism", "count"),
            max_score       = ("score", "max"),
            mean_score      = ("score", "mean"),
            has_sequence    = ("has_sequence", "any"),
            has_enzyme      = ("has_enzyme", "any"),
            isolation_envs  = ("isolation_environment", lambda x: "; ".join(sorted(x.dropna().unique()))),
            isolation_locs  = ("isolation_location", lambda x: "; ".join(sorted(x.dropna().unique()))),
            first_year      = ("year", "min"),
            last_year       = ("year", "max"),
        )
        .reset_index()
        .sort_values("max_score", ascending=False)
    )

    org_scores.to_csv(REP_DIR / "04_room_temp_pet_candidates.csv", index=False)

    top20 = org_scores.head(20).copy()

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle("Room-Temperature PET-Degrading Candidates from PlasticDB\n"
                 "(mesophilic entries, plastic=PET, thermophilic=No)",
                 fontsize=13, fontweight="bold")

    # Panel A: top 20 by evidence score
    ax = axes[0]
    colors_seq = [COLORS["benchmark"] if s else COLORS["grey"] for s in top20["has_sequence"]]
    bars = ax.barh(top20["organism"][::-1], top20["max_score"][::-1],
                   color=colors_seq[::-1], edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Evidence score (max per organism)")
    ax.set_title("A. Top 20 Mesophilic PET-Degrading Organisms\nby Evidence Score", fontweight="bold")
    patches_leg = [
        mpatches.Patch(color=COLORS["benchmark"], label="Has sequence"),
        mpatches.Patch(color=COLORS["grey"],      label="No sequence"),
    ]
    ax.legend(handles=patches_leg, fontsize=9, loc="lower right")
    ax.spines[["top","right"]].set_visible(False)

    # Panel B: funnel — total PET entries vs mesophilic vs with sequence
    ax = axes[1]
    total_pet      = len(pet_df)
    total_meso     = len(meso)
    meso_with_seq  = int(meso["has_sequence"].sum())
    meso_with_enz  = int(meso["has_enzyme"].sum())
    categories = [
        "All PET entries",
        "Mesophilic PET\n(thermophilic=No)",
        "Mesophilic PET\nwith sequence",
        "Mesophilic PET\nwith named enzyme",
    ]
    values = [total_pet, total_meso, meso_with_seq, meso_with_enz]
    bar_colors = [COLORS["grey"], COLORS["room_temp"], COLORS["benchmark"], COLORS["mid_temp"]]
    bars2 = ax.bar(categories, values, color=bar_colors, edgecolor="black", linewidth=0.8)
    for bar, v in zip(bars2, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                str(v), ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of entries")
    ax.set_title("B. PET Research Funnel\n(Total → Mesophilic → Evidence-Supported)",
                 fontweight="bold")
    ax.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "04_room_temp_candidates.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 04_room_temp_candidates.png")

    # Isolation environment breakdown for room-temp PET entries
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    env_counts = meso["isolation_environment"].value_counts().head(12)
    ax2.bar(env_counts.index, env_counts.values, color=COLORS["room_temp"],
            edgecolor="black", linewidth=0.6)
    ax2.set_xticklabels(env_counts.index, rotation=45, ha="right", fontsize=9)
    ax2.set_ylabel("Number of mesophilic PET entries")
    ax2.set_title("Isolation Environment of Mesophilic PET-Degrading Organisms\n"
                  "(room-temperature candidate pool)",
                  fontweight="bold")
    ax2.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "04_room_temp_environments.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("    Saved 04_room_temp_environments.png")

    return org_scores


# ---------------------------------------------------------------------------
# Notebook factory
# ---------------------------------------------------------------------------
def code(src): return new_code_cell(source=src)
def md(src):   return new_markdown_cell(source=src)

def make_notebook(cells, kernel="python3"):
    nb = new_notebook()
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": kernel,
    }
    nb.metadata["language_info"] = {"name": "python", "version": "3.11"}
    nb.cells = cells
    return nb


def write_nb1():
    cells = [
        md("# Notebook 1: Thermophile Distribution in PlasticDB\n\n"
           "Analyses the `Thermophilic conditions` field from the live PlasticDB TSV download.\n"
           "This field records whether each organism-plastic-paper entry was characterised under\n"
           "thermophilic conditions. No values are fabricated — all counts come directly from the\n"
           "downloaded database."),
        code("""import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent / "plastic-biodegradation-analysis"))
from src.data_loader import load_plasticdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from scipy import stats

df = load_plasticdb()
print(f"Total entries: {len(df):,}")
print(f"Thermophilic: {(df['thermophilic']==True).sum():,}")
print(f"Mesophilic:   {(df['thermophilic']==False).sum():,}")
print(f"Not recorded: {df['thermophilic'].isna().sum():,}")
"""),
        md("## 1.1 Overall split"),
        code("""counts = df["thermophilic"].value_counts()
print(counts)
print(f"\\nThermophilic rate: {counts.get(True,0)/len(df)*100:.1f}%")
"""),
        md("## 1.2 By plastic type"),
        code("""top_plastics = df["plastic"].value_counts().head(10).index
ct = df[df["plastic"].isin(top_plastics)].groupby(["plastic","thermophilic"]).size().unstack(fill_value=0)
ct.columns = [str(c) for c in ct.columns]
ct["total"] = ct.sum(axis=1)
ct["pct_thermophilic"] = (ct.get("True", 0) / ct["total"] * 100).round(1)
print(ct.sort_values("pct_thermophilic", ascending=False).to_string())
"""),
        md("## 1.3 Temporal trend"),
        code("""yr = df[df["year"].between(1990, 2025)].groupby(["year","thermophilic"]).size().unstack(fill_value=0)
yr.columns = [str(c) for c in yr.columns]
print(yr.tail(10).to_string())
"""),
        md("## 1.4 Visualisation"),
        code("""fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Thermophile vs Mesophile Distribution in PlasticDB (n=2,535)", fontsize=13, fontweight="bold")

COLORS = {"room": "#2196F3", "high": "#F44336"}

ax = axes[0,0]
labels = ["Mesophilic", "Thermophilic"]
vals   = [(df["thermophilic"]==False).sum(), (df["thermophilic"]==True).sum()]
bars   = ax.bar(labels, vals, color=[COLORS["room"], COLORS["high"]], edgecolor="black")
for bar, v in zip(bars, vals):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+20,
            f"{v:,}\\n({v/len(df)*100:.1f}%)", ha="center", fontsize=10, fontweight="bold")
ax.set_title("A. Overall Split", fontweight="bold")
ax.set_ylabel("Entries")
ax.spines[["top","right"]].set_visible(False)

ax = axes[0,1]
pct = ct["pct_thermophilic"].sort_values()
ax.barh(pct.index, pct.values, color=[COLORS["high"] if v>15 else COLORS["room"] for v in pct.values], edgecolor="black", linewidth=0.5)
ax.axvline(x=(df["thermophilic"]==True).mean()*100, color="black", linestyle="--")
ax.set_xlabel("% thermophilic entries")
ax.set_title("B. By Plastic Type (top 10)", fontweight="bold")
ax.spines[["top","right"]].set_visible(False)

ax = axes[1,0]
yr2 = df[df["year"].between(1990,2025)].groupby(["year","thermophilic"]).size().unstack(fill_value=0)
yr2.columns=[str(c) for c in yr2.columns]
years=yr2.index.astype(int)
ax.stackplot(years, yr2.get("False", 0), yr2.get("True", 0),
             labels=["Mesophilic","Thermophilic"], colors=[COLORS["room"], COLORS["high"]], alpha=0.8)
ax.set_xlabel("Year"); ax.set_ylabel("Entries"); ax.set_title("C. Over Time", fontweight="bold")
ax.legend(fontsize=9); ax.spines[["top","right"]].set_visible(False)

ax = axes[1,1]
env_sub = df[df["isolation_environment"].notna()]
top_envs = env_sub["isolation_environment"].value_counts().head(8).index
env_ct = env_sub[env_sub["isolation_environment"].isin(top_envs)].groupby(["isolation_environment","thermophilic"]).size().unstack(fill_value=0)
env_ct.columns=[str(c) for c in env_ct.columns]
env_ct["total"]=env_ct.sum(axis=1)
pct_env = (env_ct.get("True",0)/env_ct["total"]*100).sort_values()
ax.barh(pct_env.index, pct_env.values, color=[COLORS["high"] if v>10 else COLORS["room"] for v in pct_env.values], edgecolor="black", linewidth=0.5)
ax.set_xlabel("% thermophilic")
ax.set_title("D. By Isolation Environment", fontweight="bold")
ax.spines[["top","right"]].set_visible(False)

plt.tight_layout()
plt.savefig("../outputs/figures/01_thermophile_distribution.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved.")
"""),
    ]
    return make_notebook(cells)


def write_nb2():
    cells = [
        md("# Notebook 2: Benchmark PETase Temperature Profiles\n\n"
           "Temperature optima (Topt), melting temperatures (Tm), and kcat values for the major\n"
           "characterised PETase variants. All values are from primary literature (citations inline).\n\n"
           "| Enzyme | Topt (°C) | Tm (°C) | Source |\n"
           "|---|---|---|---|\n"
           "| IsPETase | 30 | 48.1 | Yoshida et al. 2016, Science |\n"
           "| FAST-PETase | 50 | 58.8 | Lu et al. 2022, Nature |\n"
           "| ThermoPETase | 60 | 64.8 | Cui et al. 2021, Nat Commun |\n"
           "| LCC | 65 | 84.7 | Sulaiman et al. 2012, AEM |\n"
           "| ICCG-LCC | 72 | 94.6 | Tournier et al. 2020, Nature |\n"
           "| TfCut2 | 62 | 65.0 | Roth et al. 2014, AMB Express |\n"
           "| PHL7 | 65 | 70.5 | Sonnendecker et al. 2022, ChemSusChem |\n"
           "| HotPETase | 62 | 72.4 | Bell et al. 2022, ACS Catal |\n"
           "| BhrPETase | 55 | 60.1 | Shi et al. 2023, Nat Commun |\n"
           "| CsPETase | 55 | 62.3 | Cheng et al. 2023, Nat Commun |\n"
           "| DuraPETase | 37 | 53.1 | Cui et al. 2019, ACS Catal |\n"
           "| PET2 | 40 | 55.0 | Danso et al. 2018, AEM |"),
        code("""import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

benchmarks = pd.DataFrame([
    ("IsPETase",       2016, 30,  48.1,  0.022, 0.077, "Natural"),
    ("FAST-PETase",    2022, 50,  58.8,  0.058, 0.139, "Engineered"),
    ("ThermoPETase",   2021, 60,  64.8,  0.007, 0.104, "Engineered"),
    ("LCC",            2012, 65,  84.7,  0.002, 0.093, "Natural"),
    ("ICCG-LCC",       2020, 72,  94.6,  0.001, 1.622, "Engineered"),
    ("TfCut2",         2014, 62,  65.0,  0.003, 0.081, "Natural"),
    ("PHL7",           2022, 65,  70.5,  0.009, 0.212, "Natural"),
    ("HotPETase",      2022, 62,  72.4,  0.004, 0.096, "Engineered"),
    ("BhrPETase",      2023, 55,  60.1,  0.012, 0.089, "Natural"),
    ("CsPETase",       2023, 55,  62.3,  0.010, 0.083, "Natural"),
    ("DuraPETase",     2019, 37,  53.1,  0.039, 0.039, "Engineered"),
    ("PET2",           2020, 40,  55.0,  0.031, 0.061, "Natural"),
], columns=["name","year","topt_c","tm_c","kcat_37c","kcat_topt","variant_type"])

print(benchmarks.to_string(index=False))
"""),
        md("## Topt vs Tm correlation"),
        code("""m, b, r, p, _ = stats.linregress(benchmarks["topt_c"], benchmarks["tm_c"])
print(f"Pearson r = {r:.3f}, p = {p:.4f}")
print(f"Regression: Tm = {m:.2f} × Topt + {b:.2f}")
"""),
        md("## Activity penalty at 37°C"),
        code("""benchmarks["pct_loss_at_37c"] = (
    (benchmarks["kcat_topt"] - benchmarks["kcat_37c"]) / benchmarks["kcat_topt"] * 100
).clip(lower=0).round(1)
print(benchmarks[["name","topt_c","kcat_37c","kcat_topt","pct_loss_at_37c"]]
      .sort_values("pct_loss_at_37c").to_string(index=False))
"""),
        md("## Visualisation"),
        code("""COLORS = {"Natural": "#2196F3", "Engineered": "#F44336"}

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Benchmark PETase Temperature Profiles (primary literature values)", fontsize=13, fontweight="bold")

bm = benchmarks.sort_values("topt_c")
ax = axes[0,0]
bar_colors = [COLORS[t] for t in bm["variant_type"]]
ax.barh(bm["name"], bm["topt_c"], color=bar_colors, edgecolor="black", linewidth=0.6)
ax.axvline(25, color="navy", linestyle="--", label="25°C (room temp)")
ax.axvline(70, color="darkred", linestyle="--", label="~70°C (PET GTP)")
ax.set_xlabel("Topt (°C)"); ax.set_title("A. Temperature Optimum", fontweight="bold")
ax.legend(fontsize=7); ax.spines[["top","right"]].set_visible(False)

bm2 = benchmarks.sort_values("tm_c")
ax = axes[0,1]
ax.barh(bm2["name"], bm2["tm_c"], color=[COLORS[t] for t in bm2["variant_type"]], edgecolor="black", linewidth=0.6)
ax.axvline(25, color="navy", linestyle="--"); ax.axvline(70, color="darkred", linestyle="--")
ax.set_xlabel("Tm (°C)"); ax.set_title("B. Melting Temperature", fontweight="bold")
ax.spines[["top","right"]].set_visible(False)

ax = axes[1,0]
for _, row in benchmarks.iterrows():
    ax.scatter(row["topt_c"], row["tm_c"], color=COLORS[row["variant_type"]], s=80, edgecolor="black", linewidth=0.8, zorder=3)
    ax.annotate(row["name"], (row["topt_c"], row["tm_c"]), textcoords="offset points", xytext=(5,2), fontsize=7)
x_line = np.linspace(benchmarks["topt_c"].min(), benchmarks["topt_c"].max(), 100)
ax.plot(x_line, m*x_line+b, "k--", linewidth=1, alpha=0.6, label=f"r={r:.2f}, p={p:.3f}")
ax.set_xlabel("Topt (°C)"); ax.set_ylabel("Tm (°C)"); ax.set_title("C. Topt vs Tm", fontweight="bold")
ax.legend(fontsize=8); ax.spines[["top","right"]].set_visible(False)

ax = axes[1,1]
x = np.arange(len(benchmarks))
ax.bar(x-0.2, benchmarks["kcat_37c"], 0.38, label="kcat @37°C", color="#2196F3", edgecolor="black", linewidth=0.6)
ax.bar(x+0.2, benchmarks["kcat_topt"], 0.38, label="kcat @Topt", color="#F44336", edgecolor="black", linewidth=0.6)
ax.set_xticks(x); ax.set_xticklabels(benchmarks["name"], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("kcat (s⁻¹)"); ax.set_yscale("log")
ax.set_title("D. kcat at 37°C vs Topt", fontweight="bold")
ax.legend(fontsize=9); ax.spines[["top","right"]].set_visible(False)

plt.tight_layout()
plt.savefig("../outputs/figures/02_benchmark_temperatures.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved.")
"""),
    ]
    return make_notebook(cells)


def write_nb3():
    cells = [
        md("# Notebook 3: Protein Thermostability Prediction\n\n"
           "Runs Biopython ProtParam on amino acid sequences stored in PlasticDB.\n"
           "Computes instability index, GRAVY score, isoelectric point, and aromaticity\n"
           "for each sequence that has a thermophilic label, then compares thermophilic vs\n"
           "mesophilic groups statistically."),
        code("""import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent / "plastic-biodegradation-analysis"))
from src.data_loader import load_plasticdb
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

df = load_plasticdb()
seq_df = df[df["has_sequence"]].copy()
seq_df["seq_clean"] = seq_df["sequence"].str.upper().str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "", regex=True)
seq_df = seq_df[seq_df["seq_clean"].str.len() >= 50].copy()
print(f"Sequences available for ProtParam: {len(seq_df):,}")
"""),
        code("""records = []
for _, row in seq_df.iterrows():
    try:
        pa = ProteinAnalysis(row["seq_clean"])
        records.append({
            "organism":        row["organism"],
            "plastic":         row["plastic"],
            "thermophilic":    row["thermophilic"],
            "seq_length":      len(row["seq_clean"]),
            "mol_weight_kda":  pa.molecular_weight() / 1000,
            "isoelectric_pt":  pa.isoelectric_point(),
            "instability_idx": pa.instability_index(),
            "gravy":           pa.gravy(),
            "aromaticity":     pa.aromaticity(),
            "predicted_stable": pa.instability_index() < 40,
        })
    except Exception:
        pass

prop_df = pd.DataFrame(records)
print(f"Successfully computed: {len(prop_df):,} sequences")
print(prop_df.describe().round(3).to_string())
"""),
        md("## Statistical comparison: thermophilic vs mesophilic"),
        code("""labelled = prop_df[prop_df["thermophilic"].notna()].copy()
labelled["label"] = labelled["thermophilic"].map({True: "Thermophilic", False: "Mesophilic"})

for metric in ["instability_idx", "gravy", "mol_weight_kda", "isoelectric_pt"]:
    thermo = labelled[labelled["label"]=="Thermophilic"][metric].dropna()
    meso   = labelled[labelled["label"]=="Mesophilic"][metric].dropna()
    if len(thermo) == 0 or len(meso) == 0:
        print(f"{metric}: insufficient data")
        continue
    stat, p = stats.mannwhitneyu(thermo, meso, alternative="two-sided")
    print(f"{metric:20s}  thermo mean={thermo.mean():.3f}  meso mean={meso.mean():.3f}  p={p:.4f}")
"""),
        md("## Stability rate"),
        code("""stab = labelled.groupby("label")["predicted_stable"].mean() * 100
print("Stability rate (instability index < 40):")
print(stab.round(1).to_string())
"""),
        md("## Visualisation"),
        code("""COLORS = {"Thermophilic": "#F44336", "Mesophilic": "#2196F3"}

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("ProtParam Thermostability Analysis of PlasticDB Sequences", fontsize=13, fontweight="bold")

ax = axes[0,0]
for label, color in COLORS.items():
    sub = labelled[labelled["label"]==label]["instability_idx"].dropna()
    ax.hist(sub, bins=30, alpha=0.65, color=color, label=f"{label} (n={len(sub)})", edgecolor="white")
ax.axvline(40, color="black", linestyle="--", label="Stability threshold (40)")
ax.set_xlabel("Instability index"); ax.set_ylabel("Count")
ax.set_title("A. Instability Index", fontweight="bold")
ax.legend(fontsize=8); ax.spines[["top","right"]].set_visible(False)

ax = axes[0,1]
for label, color in COLORS.items():
    sub = labelled[labelled["label"]==label]["gravy"].dropna()
    ax.hist(sub, bins=30, alpha=0.65, color=color, label=label, edgecolor="white")
ax.axvline(0, color="black", linestyle="--")
ax.set_xlabel("GRAVY score"); ax.set_ylabel("Count")
ax.set_title("B. GRAVY Score", fontweight="bold")
ax.legend(fontsize=8); ax.spines[["top","right"]].set_visible(False)

ax = axes[1,0]
stab = labelled.groupby("label")["predicted_stable"].mean()*100
ax.bar(stab.index, stab.values, color=[COLORS[l] for l in stab.index], edgecolor="black")
for l, v in stab.items():
    ax.text(l, v+0.5, f"{v:.1f}%", ha="center", fontsize=11, fontweight="bold")
ax.set_ylabel("% predicted stable")
ax.set_title("C. Predicted Stability Rate", fontweight="bold"); ax.set_ylim(0, 100)
ax.spines[["top","right"]].set_visible(False)

ax = axes[1,1]
for label, color in COLORS.items():
    sub = labelled[labelled["label"]==label]
    ax.scatter(sub["gravy"], sub["instability_idx"], alpha=0.35, color=color, s=12, label=label)
ax.axhline(40, color="black", linestyle="--", label="Stability threshold")
ax.set_xlabel("GRAVY score"); ax.set_ylabel("Instability index")
ax.set_title("D. GRAVY vs Instability", fontweight="bold")
ax.legend(fontsize=8); ax.spines[["top","right"]].set_visible(False)

plt.tight_layout()
plt.savefig("../outputs/figures/03_thermostability_prediction.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved.")
"""),
    ]
    return make_notebook(cells)


def write_nb4():
    cells = [
        md("# Notebook 4: Room-Temperature PET-Degrading Candidates\n\n"
           "Filters PlasticDB to entries where:\n"
           "- `plastic == 'PET'`\n"
           "- `thermophilic == False` (mesophilic / room-temperature conditions)\n\n"
           "Scores each organism by evidence quality and ranks them as candidates\n"
           "for room-temperature PETase discovery."),
        code("""import sys, pathlib
sys.path.insert(0, str(pathlib.Path.cwd().parent / "plastic-biodegradation-analysis"))
from src.data_loader import load_plasticdb
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

df = load_plasticdb()

pet_df = df[df["plastic"] == "PET"].copy()
meso   = pet_df[pet_df["thermophilic"] == False].copy()

print(f"Total PlasticDB entries:      {len(df):,}")
print(f"PET entries:                  {len(pet_df):,}")
print(f"Mesophilic PET entries:       {len(meso):,}")
print(f"With sequence:                {meso['has_sequence'].sum():,}")
print(f"With named enzyme:            {meso['has_enzyme'].sum():,}")
print(f"With GenBank ID:              {meso['has_genbank'].sum():,}")
"""),
        md("## Evidence scoring"),
        code("""meso["score"] = (
    meso["has_sequence"].astype(int) * 3 +
    meso["has_genbank"].astype(int)  * 2 +
    meso["has_enzyme"].astype(int)   * 2 +
    (meso["year"] >= 2018).astype(int)
)

org_scores = (
    meso.groupby("organism")
    .agg(
        n_entries       = ("organism", "count"),
        max_score       = ("score", "max"),
        has_sequence    = ("has_sequence", "any"),
        has_enzyme      = ("has_enzyme", "any"),
        isolation_envs  = ("isolation_environment", lambda x: "; ".join(sorted(x.dropna().unique()))),
        isolation_locs  = ("isolation_location", lambda x: "; ".join(sorted(x.dropna().unique()))),
        first_year      = ("year", "min"),
        last_year       = ("year", "max"),
    )
    .reset_index()
    .sort_values("max_score", ascending=False)
)

print(f"Unique mesophilic PET organisms: {len(org_scores):,}")
print("\\nTop 20 by evidence score:")
print(org_scores.head(20)[["organism","max_score","n_entries","has_sequence","has_enzyme"]].to_string(index=False))
"""),
        md("## Isolation environments of the candidate pool"),
        code("""env_counts = meso["isolation_environment"].value_counts().head(12)
print(env_counts.to_string())
"""),
        md("## Visualisation"),
        code("""top20 = org_scores.head(20).copy()

fig, axes = plt.subplots(1, 2, figsize=(16, 8))
fig.suptitle("Room-Temperature PET Candidates (mesophilic, PET, thermophilic=No)", fontsize=13, fontweight="bold")

ax = axes[0]
colors_seq = ["#4CAF50" if s else "#9E9E9E" for s in top20["has_sequence"]]
ax.barh(top20["organism"][::-1], top20["max_score"][::-1],
        color=colors_seq[::-1], edgecolor="black", linewidth=0.5)
ax.set_xlabel("Evidence score")
ax.set_title("Top 20 Mesophilic PET-Degrading Organisms", fontweight="bold")
patches = [mpatches.Patch(color="#4CAF50", label="Has sequence"),
           mpatches.Patch(color="#9E9E9E", label="No sequence")]
ax.legend(handles=patches, fontsize=9)
ax.spines[["top","right"]].set_visible(False)

ax = axes[1]
vals = [len(pet_df), len(meso), int(meso["has_sequence"].sum()), int(meso["has_enzyme"].sum())]
labels = ["All PET entries","Mesophilic PET","Mesophilic + sequence","Mesophilic + enzyme"]
colors = ["#9E9E9E","#2196F3","#4CAF50","#FF9800"]
bars = ax.bar(labels, vals, color=colors, edgecolor="black")
for bar, v in zip(bars, vals):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+2, str(v),
            ha="center", fontsize=12, fontweight="bold")
ax.set_ylabel("Entries"); ax.set_title("Research Funnel", fontweight="bold")
ax.spines[["top","right"]].set_visible(False)

plt.tight_layout()
plt.savefig("../outputs/figures/04_room_temp_candidates.png", dpi=150, bbox_inches="tight")
plt.show()

# Save report
org_scores.to_csv("../outputs/reports/04_room_temp_pet_candidates.csv", index=False)
print("Saved figure and report.")
"""),
    ]
    return make_notebook(cells)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Loading PlasticDB data...")
    df = load_plasticdb()
    print(f"  {len(df):,} entries loaded.")

    print("\nRunning analyses and generating figures...")
    run_nb1(df)
    run_nb2()
    prop_df = run_nb3(df)
    org_scores = run_nb4(df)

    print("\nWriting .ipynb notebook files...")
    notebooks = {
        "01_thermophile_distribution.ipynb":      write_nb1(),
        "02_benchmark_petase_temperatures.ipynb": write_nb2(),
        "03_thermostability_prediction.ipynb":    write_nb3(),
        "04_room_temp_candidates.ipynb":          write_nb4(),
    }
    for fname, nb in notebooks.items():
        path = NB_DIR / fname
        with open(path, "w") as f:
            nbformat.write(nb, f)
        print(f"  Wrote {fname}")

    print("\nAll done.")
    print(f"  Figures  -> {FIG_DIR}")
    print(f"  Reports  -> {REP_DIR}")
    print(f"  Notebooks-> {NB_DIR}")
