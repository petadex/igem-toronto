"""
Visualization functions for plastic biodegradation meta-analysis.

All functions return a matplotlib Figure or a plotly Figure.
Set `backend='plotly'` for interactive charts, `backend='matplotlib'` for publication-ready.
"""

from __future__ import annotations
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

FIGURE_DIR = Path(__file__).parent.parent / "outputs" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "primary":    "#2E86AB",
    "secondary":  "#A23B72",
    "accent":     "#F18F01",
    "green":      "#44BBA4",
    "red":        "#E94F37",
    "background": "#F7F9FC",
}

sns.set_theme(style="whitegrid", palette="deep", font_scale=1.1)


def _save(fig, name: str, backend: str):
    if backend == "matplotlib":
        path = FIGURE_DIR / f"{name}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
        return path
    else:
        path = FIGURE_DIR / f"{name}.html"
        fig.write_html(str(path))
        return path


def plot_plastic_distribution(df: pd.DataFrame, top_n: int = 20, backend: str = "plotly"):
    """Bar chart of the most-studied plastic types."""
    counts = df["plastic"].value_counts().head(top_n).reset_index()
    counts.columns = ["plastic", "count"]

    if backend == "plotly":
        fig = px.bar(
            counts, x="plastic", y="count",
            color="count", color_continuous_scale="Blues",
            title=f"Top {top_n} Plastic Types by Number of Reported Degraders",
            labels={"plastic": "Plastic Type", "count": "Number of Entries"},
            text="count",
        )
        fig.update_traces(textposition="outside")
        fig.update_layout(
            coloraxis_showscale=False,
            plot_bgcolor=PALETTE["background"],
            height=480,
        )
        _save(fig, "plastic_distribution", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(12, 5))
        bars = ax.bar(counts["plastic"], counts["count"],
                      color=sns.color_palette("Blues_r", len(counts)))
        ax.set_xlabel("Plastic Type", fontsize=12)
        ax.set_ylabel("Number of Entries", fontsize=12)
        ax.set_title(f"Top {top_n} Plastic Types by Number of Reported Degraders", fontsize=14)
        ax.tick_params(axis="x", rotation=45)
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                    str(int(bar.get_height())), ha="center", va="bottom", fontsize=8)
        plt.tight_layout()
        _save(fig, "plastic_distribution", backend)
        return fig


def plot_temporal_trends(temporal_df: pd.DataFrame, backend: str = "plotly"):
    """Dual-axis chart: entries per year + cumulative unique species."""
    if backend == "plotly":
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(
            go.Bar(x=temporal_df["year"], y=temporal_df["n_entries"],
                   name="Entries per Year", marker_color=PALETTE["primary"], opacity=0.7),
            secondary_y=False,
        )
        fig.add_trace(
            go.Scatter(x=temporal_df["year"], y=temporal_df["rolling_3yr"],
                       name="3-yr Rolling Avg", line=dict(color=PALETTE["accent"], width=2)),
            secondary_y=False,
        )
        fig.add_trace(
            go.Scatter(x=temporal_df["year"], y=temporal_df["cumulative_species"],
                       name="Cumulative Species", line=dict(color=PALETTE["secondary"], dash="dot", width=2)),
            secondary_y=True,
        )
        fig.update_layout(
            title="Publication Trends in Plastic Biodegradation Research",
            xaxis_title="Year",
            plot_bgcolor=PALETTE["background"],
            height=480,
            legend=dict(x=0.01, y=0.99),
        )
        fig.update_yaxes(title_text="Entries per Year", secondary_y=False)
        fig.update_yaxes(title_text="Cumulative Unique Species", secondary_y=True)
        _save(fig, "temporal_trends", backend)
        return fig
    else:
        fig, ax1 = plt.subplots(figsize=(13, 5))
        ax2 = ax1.twinx()
        ax1.bar(temporal_df["year"], temporal_df["n_entries"],
                color=PALETTE["primary"], alpha=0.6, label="Entries/yr")
        ax1.plot(temporal_df["year"], temporal_df["rolling_3yr"],
                 color=PALETTE["accent"], lw=2, label="3-yr avg")
        ax2.plot(temporal_df["year"], temporal_df["cumulative_species"],
                 color=PALETTE["secondary"], lw=2, ls="--", label="Cumulative species")
        ax1.set_xlabel("Year", fontsize=12)
        ax1.set_ylabel("Entries per Year", fontsize=12)
        ax2.set_ylabel("Cumulative Unique Species", fontsize=12)
        ax1.set_title("Publication Trends in Plastic Biodegradation Research", fontsize=14)
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
        plt.tight_layout()
        _save(fig, "temporal_trends", backend)
        return fig


def plot_geographic_heatmap(geo_df: pd.DataFrame, backend: str = "plotly"):
    """Choropleth / bubble map of research activity by country."""
    geo_df = geo_df.copy()
    geo_df["location_clean"] = geo_df["isolation_location"].str.strip()
    top_geo = geo_df.nlargest(40, "n_entries")

    if backend == "plotly":
        fig = px.bar(
            top_geo, x="n_entries", y="location_clean",
            orientation="h",
            color="n_species",
            color_continuous_scale="Teal",
            title="Research Activity by Country/Region (Top 40)",
            labels={"n_entries": "Number of Entries", "location_clean": "Country/Region",
                    "n_species": "Unique Species"},
            text="n_entries",
        )
        fig.update_traces(textposition="outside")
        fig.update_layout(
            plot_bgcolor=PALETTE["background"],
            height=900,
            yaxis={"categoryorder": "total ascending"},
        )
        _save(fig, "geographic_heatmap", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(10, 12))
        top_geo_sorted = top_geo.sort_values("n_entries")
        bars = ax.barh(top_geo_sorted["location_clean"], top_geo_sorted["n_entries"],
                       color=sns.color_palette("Blues", len(top_geo_sorted)))
        ax.set_xlabel("Number of Entries", fontsize=12)
        ax.set_title("Research Activity by Country/Region", fontsize=14)
        plt.tight_layout()
        _save(fig, "geographic_heatmap", backend)
        return fig


def plot_co_occurrence_heatmap(co_matrix: pd.DataFrame, top_n: int = 15, backend: str = "plotly"):
    """Heatmap of plastic co-occurrence by shared degrading organisms."""
    plastics = co_matrix.sum().nlargest(top_n).index
    sub_arr = co_matrix.loc[plastics, plastics].values.astype(float)
    np.fill_diagonal(sub_arr, 0)
    sub = pd.DataFrame(sub_arr, index=plastics, columns=plastics)

    if backend == "plotly":
        fig = px.imshow(
            sub,
            color_continuous_scale="Blues",
            title="Plastic Co-occurrence: Organisms that Degrade Multiple Plastics",
            labels={"color": "Shared Organisms"},
            text_auto=True,
        )
        fig.update_layout(height=600, plot_bgcolor=PALETTE["background"])
        _save(fig, "co_occurrence_heatmap", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(sub.astype(int), annot=True, fmt="d", cmap="Blues", ax=ax,
                    linewidths=0.5, cbar_kws={"label": "Shared Organisms"})
        ax.set_title("Plastic Co-occurrence Heatmap", fontsize=14)
        plt.tight_layout()
        _save(fig, "co_occurrence_heatmap", backend)
        return fig


def plot_novelty_scatter(novelty_df: pd.DataFrame, top_n: int = 30, backend: str = "plotly"):
    """
    Scatter: breadth_score vs rarity_score, bubble size = evidence_gap_score.
    Highlights top candidates for novel discovery.
    """
    top = novelty_df.head(top_n).copy()
    top["label"] = top["organism"].apply(lambda s: " ".join(s.split()[:2]))

    if backend == "plotly":
        fig = px.scatter(
            top,
            x="breadth_score",
            y="rarity_score",
            size="evidence_gap_score",
            color="novelty_potential",
            color_continuous_scale="Viridis",
            hover_name="organism",
            text="label",
            size_max=40,
            title=f"Top {top_n} Organisms by Novelty Potential",
            labels={
                "breadth_score": "Plastic Breadth Score",
                "rarity_score": "Taxonomic Rarity Score",
                "novelty_potential": "Novelty Potential",
            },
        )
        fig.update_traces(textposition="top center")
        fig.update_layout(height=580, plot_bgcolor=PALETTE["background"])
        _save(fig, "novelty_scatter", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(11, 7))
        scatter = ax.scatter(
            top["breadth_score"], top["rarity_score"],
            s=top["evidence_gap_score"] * 5 + 20,
            c=top["novelty_potential"],
            cmap="viridis", alpha=0.8,
        )
        plt.colorbar(scatter, ax=ax, label="Novelty Potential")
        for _, row in top.iterrows():
            ax.annotate(row["label"], (row["breadth_score"], row["rarity_score"]),
                        textcoords="offset points", xytext=(4, 4), fontsize=6)
        ax.set_xlabel("Plastic Breadth Score", fontsize=12)
        ax.set_ylabel("Taxonomic Rarity Score", fontsize=12)
        ax.set_title(f"Top {top_n} Organisms by Novelty Potential", fontsize=14)
        plt.tight_layout()
        _save(fig, "novelty_scatter", backend)
        return fig


def plot_evidence_quality_dist(df_scored: pd.DataFrame, backend: str = "plotly"):
    """Donut chart of evidence quality tier distribution."""
    tier_counts = df_scored["evidence_tier"].value_counts().reset_index()
    tier_counts.columns = ["tier", "count"]
    order = ["Excellent", "High", "Medium", "Low"]
    tier_counts["tier"] = pd.Categorical(tier_counts["tier"], categories=order, ordered=True)
    tier_counts = tier_counts.sort_values("tier")

    colors = ["#2E86AB", "#44BBA4", "#F18F01", "#E94F37"]

    if backend == "plotly":
        fig = px.pie(
            tier_counts, names="tier", values="count",
            color="tier",
            color_discrete_map=dict(zip(order, colors)),
            title="Evidence Quality Tier Distribution (PlasticDB)",
            hole=0.45,
        )
        fig.update_traces(textinfo="percent+label")
        fig.update_layout(height=420)
        _save(fig, "evidence_quality", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(7, 6))
        wedges, texts, autotexts = ax.pie(
            tier_counts["count"], labels=tier_counts["tier"],
            autopct="%1.1f%%", colors=colors, startangle=90,
            wedgeprops={"width": 0.55},
        )
        ax.set_title("Evidence Quality Tier Distribution", fontsize=14)
        plt.tight_layout()
        _save(fig, "evidence_quality", backend)
        return fig


def plot_genus_top20(df: pd.DataFrame, backend: str = "plotly"):
    """Horizontal bar — top 20 genera by number of entries."""
    top_genera = (
        df.groupby("genus")["organism"]
        .count()
        .nlargest(20)
        .reset_index()
        .rename(columns={"organism": "entries"})
        .sort_values("entries")
    )
    if backend == "plotly":
        fig = px.bar(
            top_genera, x="entries", y="genus", orientation="h",
            title="Top 20 Genera with Most Plastic-Degradation Entries",
            labels={"entries": "Number of Entries", "genus": "Genus"},
            color="entries", color_continuous_scale="Teal", text="entries",
        )
        fig.update_traces(textposition="outside")
        fig.update_layout(coloraxis_showscale=False,
                          plot_bgcolor=PALETTE["background"], height=600)
        _save(fig, "genus_top20", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(9, 8))
        ax.barh(top_genera["genus"], top_genera["entries"],
                color=sns.color_palette("Blues", len(top_genera)))
        ax.set_xlabel("Number of Entries")
        ax.set_title("Top 20 Genera — Plastic Biodegradation Entries", fontsize=14)
        plt.tight_layout()
        _save(fig, "genus_top20", backend)
        return fig


def plot_plastic_category_sunburst(df: pd.DataFrame, backend: str = "plotly"):
    """Sunburst of plastic category → specific plastic → entry count."""
    agg = (
        df.groupby(["plastic_category", "plastic"])
        .size()
        .reset_index(name="count")
    )
    if backend == "plotly":
        fig = px.sunburst(
            agg, path=["plastic_category", "plastic"], values="count",
            color="count", color_continuous_scale="Blues",
            title="Plastic Biodegradation Research Landscape",
        )
        fig.update_layout(height=540)
        _save(fig, "plastic_category_sunburst", backend)
        return fig
    else:
        return plot_plastic_distribution(df, backend=backend)


def plot_research_gaps(gaps_df: pd.DataFrame, backend: str = "plotly"):
    """Bar chart of plastics ranked by research gap score."""
    top = gaps_df.head(15).copy()
    if backend == "plotly":
        fig = px.bar(
            top, x="gap_score", y="plastic", orientation="h",
            color="pct_with_sequence",
            color_continuous_scale="RdYlGn",
            title="Research Gap Priority Score by Plastic Type",
            labels={"gap_score": "Gap Score (higher = more underexplored)",
                    "plastic": "Plastic Type",
                    "pct_with_sequence": "Frac. with Sequence"},
            text=top["n_organisms"].astype(str) + " spp.",
        )
        fig.update_traces(textposition="outside")
        fig.update_layout(
            plot_bgcolor=PALETTE["background"], height=580,
            yaxis={"categoryorder": "total ascending"},
        )
        _save(fig, "research_gaps", backend)
        return fig
    else:
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.barh(top["plastic"], top["gap_score"],
                color=sns.color_palette("RdYlGn_r", len(top)))
        ax.set_xlabel("Gap Score")
        ax.set_title("Research Gap Priority by Plastic Type", fontsize=14)
        plt.tight_layout()
        _save(fig, "research_gaps", backend)
        return fig


def plot_pazy_vs_plasticdb(cross_db: dict, backend: str = "plotly"):
    """
    Side-by-side bar showing PlasticDB vs PAZy coverage per plastic type.
    """
    plastics = sorted(set(cross_db["shared_plastics"]) |
                      set(cross_db["plasticdb_only_plastics"]) |
                      set(cross_db["pazy_only_plastics"]))
    in_pdb = [1 if p in cross_db["shared_plastics"] + cross_db["plasticdb_only_plastics"] else 0 for p in plastics]
    in_pazy = [1 if p in cross_db["shared_plastics"] + cross_db["pazy_only_plastics"] else 0 for p in plastics]

    if backend == "plotly":
        fig = go.Figure()
        fig.add_trace(go.Bar(name="PlasticDB", x=plastics, y=in_pdb,
                             marker_color=PALETTE["primary"]))
        fig.add_trace(go.Bar(name="PAZy", x=plastics, y=in_pazy,
                             marker_color=PALETTE["secondary"]))
        fig.update_layout(
            barmode="group",
            title="Plastic Type Coverage: PlasticDB vs PAZy",
            xaxis_title="Plastic Type",
            yaxis_title="Covered (1=Yes, 0=No)",
            plot_bgcolor=PALETTE["background"],
            height=400,
        )
        _save(fig, "pazy_vs_plasticdb", backend)
        return fig
    else:
        x = np.arange(len(plastics))
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.bar(x - 0.2, in_pdb, 0.4, label="PlasticDB", color=PALETTE["primary"])
        ax.bar(x + 0.2, in_pazy, 0.4, label="PAZy", color=PALETTE["secondary"])
        ax.set_xticks(x)
        ax.set_xticklabels(plastics, rotation=45, ha="right")
        ax.set_title("Plastic Type Coverage: PlasticDB vs PAZy", fontsize=14)
        ax.legend()
        plt.tight_layout()
        _save(fig, "pazy_vs_plasticdb", backend)
        return fig


def save_all_figures(results: dict, backend: str = "matplotlib"):
    """Generate and save every standard figure. Returns list of file paths."""
    df = results["df_scored"]
    paths = []
    paths.append(_save(plot_plastic_distribution(df, backend=backend), "plastic_distribution", backend))
    paths.append(_save(plot_temporal_trends(results["temporal_trends"], backend=backend), "temporal_trends", backend))
    paths.append(_save(plot_geographic_heatmap(results["geographic_distribution"], backend=backend), "geographic_heatmap", backend))
    paths.append(_save(plot_co_occurrence_heatmap(results["co_occurrence"], backend=backend), "co_occurrence_heatmap", backend))
    paths.append(_save(plot_novelty_scatter(results["novelty_scores"], backend=backend), "novelty_scatter", backend))
    paths.append(_save(plot_evidence_quality_dist(df, backend=backend), "evidence_quality", backend))
    paths.append(_save(plot_genus_top20(df, backend=backend), "genus_top20", backend))
    paths.append(_save(plot_plastic_category_sunburst(df, backend=backend), "plastic_category_sunburst", backend))
    paths.append(_save(plot_research_gaps(results["research_gaps"]["plastic_gaps"], backend=backend), "research_gaps", backend))
    paths.append(_save(plot_pazy_vs_plasticdb(results["cross_db"], backend=backend), "pazy_vs_plasticdb", backend))
    return paths
