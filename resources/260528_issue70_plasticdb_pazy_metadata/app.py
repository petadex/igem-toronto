"""
Plastic Biodegradation Research Explorer
Streamlit application for interactive meta-analysis.

Run: streamlit run app.py --server.port 5000
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.data_loader import load_all, PLASTIC_FULL_NAMES
from src.analysis import (
    evidence_quality_score,
    taxonomic_diversity,
    temporal_trend_analysis,
    geographic_distribution,
    isolation_environment_profile,
    research_gap_analysis,
    cross_database_comparison,
    compute_novelty_potential,
    plastic_co_occurrence,
)
from src.visualization import (
    plot_plastic_distribution,
    plot_temporal_trends,
    plot_geographic_heatmap,
    plot_co_occurrence_heatmap,
    plot_novelty_scatter,
    plot_evidence_quality_dist,
    plot_genus_top20,
    plot_plastic_category_sunburst,
    plot_research_gaps,
    plot_pazy_vs_plasticdb,
)
from src.novel_discovery import (
    identify_phylogenetic_gaps,
    underexplored_environments,
    generate_discovery_report,
)

st.set_page_config(
    page_title="Plastic Biodegradation Explorer",
    page_icon="♻️",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_data(show_spinner="Loading databases...")
def get_data():
    return load_all()


@st.cache_data(show_spinner="Running analysis...")
def get_results(_data):
    df = _data["plasticdb"]
    organisms = _data["organisms"]
    pazy = _data["pazy"]
    df_scored = evidence_quality_score(df)
    return {
        "df_scored": df_scored,
        "taxonomic_diversity": taxonomic_diversity(df),
        "temporal_trends": temporal_trend_analysis(df),
        "geographic_distribution": geographic_distribution(df),
        "isolation_environments": isolation_environment_profile(df),
        "research_gaps": research_gap_analysis(df),
        "co_occurrence": plastic_co_occurrence(organisms),
        "cross_db": cross_database_comparison(df, pazy),
        "novelty_scores": compute_novelty_potential(organisms, df),
    }


data = get_data()
results = get_results(data)
df = data["plasticdb"]
organisms = data["organisms"]
pazy = data["pazy"]
df_scored = results["df_scored"]


with st.sidebar:
    st.title("♻️ Plastic Biodegradation")
    st.caption("Meta-analysis | PlasticDB + PAZy")
    st.divider()
    page = st.radio(
        "Navigation",
        [
            "📊 Overview",
            "🧬 Taxonomy",
            "🧪 Plastic Substrates",
            "🗺️ Geography & Time",
            "🔍 Research Gaps",
            "💡 Novel Discovery",
            "🔗 Cross-Database",
            "📋 Data Explorer",
        ],
    )
    st.divider()
    st.caption(
        "Sources: [PlasticDB](https://plasticdb.org) · [PAZy](https://www.pazy.eu)\n\n"
        "Gambarini et al. 2022; Buchholz et al. 2022"
    )


if page == "📊 Overview":
    st.title("Plastic Biodegradation Research — Meta-Analysis")
    st.markdown(
        "Comprehensive analysis of **PlasticDB** (875 species, 329 proteins) and "
        "**PAZy** (thoroughly characterised plastic-active enzymes). "
        "Explore taxonomic diversity, substrate coverage, geographic trends, "
        "evidence quality, and candidate organisms for novel discovery."
    )

    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Total Entries", f"{len(df):,}")
    col2.metric("Unique Species", f"{df['organism'].nunique():,}")
    col3.metric("Unique Genera", f"{df['genus'].nunique():,}")
    col4.metric("Plastic Types", f"{df['plastic'].nunique():,}")
    col5.metric("Years Covered", f"{int(df['year'].min()):.0f}–{int(df['year'].max()):.0f}")

    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        st.plotly_chart(plot_plastic_distribution(df, backend="plotly"),
                        use_container_width=True)
    with col_b:
        st.plotly_chart(plot_evidence_quality_dist(df_scored, backend="plotly"),
                        use_container_width=True)

    st.plotly_chart(plot_temporal_trends(results["temporal_trends"], backend="plotly"),
                    use_container_width=True)

    td = results["taxonomic_diversity"]
    st.subheader("Key Metrics")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Unique Genera", td["n_unique_genera"])
    c2.metric("Shannon Diversity (genus)", f"{td['genus_shannon_diversity']:.3f}")
    c3.metric("Singleton Genera", td["singleton_genera"])
    c4.metric("% Singleton", f"{td['pct_singleton_genera']}%")


elif page == "🧬 Taxonomy":
    st.title("Taxonomic Analysis")

    tab1, tab2, tab3 = st.tabs(["Top Genera", "Species Distribution", "Plastic Breadth"])

    with tab1:
        st.plotly_chart(plot_genus_top20(df, backend="plotly"), use_container_width=True)
        td = results["taxonomic_diversity"]
        st.subheader("Top 10 Genera")
        top_g = pd.DataFrame(list(td["top_10_genera"].items()), columns=["Genus", "Entries"])
        st.dataframe(top_g, use_container_width=True, hide_index=True)

    with tab2:
        st.subheader("Top 30 Species by Number of Entries")
        top_sp = df["organism"].value_counts().head(30).reset_index()
        top_sp.columns = ["Organism", "Entries"]
        fig = px.bar(top_sp, x="Entries", y="Organism", orientation="h",
                     color="Entries", color_continuous_scale="Teal",
                     title="Top 30 Plastic-Degrading Species")
        fig.update_layout(yaxis={"categoryorder": "total ascending"},
                          height=700, coloraxis_showscale=False)
        st.plotly_chart(fig, use_container_width=True)

    with tab3:
        st.subheader("Organisms Degrading the Most Plastic Types")
        multi = (
            organisms[organisms["n_plastics"] >= 3]
            .sort_values("n_plastics", ascending=False)
            .head(40)
            [["organism", "genus", "n_plastics", "n_entries", "plastics_degraded",
              "first_year", "last_year", "has_sequence", "has_enzyme"]]
        )
        multi["plastics_degraded"] = multi["plastics_degraded"].apply(
            lambda x: ", ".join(x) if isinstance(x, list) else x
        )
        st.dataframe(multi, use_container_width=True, hide_index=True)


elif page == "🧪 Plastic Substrates":
    st.title("Plastic Substrate Analysis")

    col_a, col_b = st.columns(2)
    with col_a:
        st.plotly_chart(plot_plastic_category_sunburst(df, backend="plotly"),
                        use_container_width=True)
    with col_b:
        st.subheader("Plastic Summary Table")
        ps = data["plastics"][["plastic", "full_name", "category",
                                "n_entries", "n_unique_organisms",
                                "n_unique_genera", "pct_with_sequence",
                                "first_year", "last_year"]]
        st.dataframe(ps, use_container_width=True, hide_index=True)

    st.subheader("Plastic Co-occurrence (Shared Degrading Organisms)")
    st.caption("How often are pairs of plastics degraded by the same organism?")
    st.plotly_chart(plot_co_occurrence_heatmap(results["co_occurrence"], backend="plotly"),
                    use_container_width=True)

    st.subheader("Filter by Plastic Type")
    selected_plastic = st.selectbox("Select plastic", sorted(df["plastic"].dropna().unique()))
    plastic_df = df[df["plastic"] == selected_plastic]
    c1, c2, c3 = st.columns(3)
    c1.metric("Entries", len(plastic_df))
    c2.metric("Unique Species", plastic_df["organism"].nunique())
    c3.metric("% with Sequence", f"{plastic_df['has_sequence'].mean()*100:.1f}%")
    top_sp = plastic_df["organism"].value_counts().head(15).reset_index()
    top_sp.columns = ["Organism", "Entries"]
    fig2 = px.bar(top_sp, x="Organism", y="Entries",
                  title=f"Top 15 Organisms Degrading {selected_plastic}",
                  color="Entries", color_continuous_scale="Blues")
    fig2.update_layout(coloraxis_showscale=False)
    st.plotly_chart(fig2, use_container_width=True)


elif page == "🗺️ Geography & Time":
    st.title("Geographic & Temporal Trends")

    st.plotly_chart(plot_temporal_trends(results["temporal_trends"], backend="plotly"),
                    use_container_width=True)

    col_a, col_b = st.columns(2)
    with col_a:
        st.plotly_chart(plot_geographic_heatmap(results["geographic_distribution"],
                                                backend="plotly"),
                        use_container_width=True)
    with col_b:
        st.subheader("Isolation Environment Breakdown")
        env_df = results["isolation_environments"]
        fig_env = px.bar(
            env_df.head(20), x="n_entries", y="isolation_environment",
            orientation="h", color="n_species", color_continuous_scale="Teal",
            title="Entries by Isolation Environment",
            labels={"n_entries": "Entries", "isolation_environment": "Environment",
                    "n_species": "Species"},
        )
        fig_env.update_layout(yaxis={"categoryorder": "total ascending"},
                               height=600)
        st.plotly_chart(fig_env, use_container_width=True)

    st.subheader("Year-by-Year Data Table")
    tt = results["temporal_trends"]
    tt_display = tt[["year", "n_entries", "n_unique_species", "n_unique_plastics",
                      "rolling_3yr", "yoy_growth_pct", "cumulative_entries",
                      "cumulative_species"]].copy()
    tt_display["year"] = tt_display["year"].astype(int)
    tt_display["rolling_3yr"] = tt_display["rolling_3yr"].round(1)
    tt_display["yoy_growth_pct"] = tt_display["yoy_growth_pct"].round(1)
    st.dataframe(tt_display, use_container_width=True, hide_index=True)


elif page == "🔍 Research Gaps":
    st.title("Research Gap Analysis")
    st.markdown(
        "Gap scores identify plastics where degradation research is sparse, "
        "lacks molecular evidence, or relies on extrapolation. "
        "Higher score = greater research opportunity."
    )

    st.plotly_chart(plot_research_gaps(results["research_gaps"]["plastic_gaps"],
                                       backend="plotly"),
                    use_container_width=True)

    col_a, col_b = st.columns(2)
    with col_a:
        st.subheader("Plastics with No Recent Data (post-2020)")
        no_recent = results["research_gaps"]["no_recent_data_plastics"]
        if no_recent:
            st.write(", ".join(sorted(no_recent)))
        else:
            st.write("All plastics have recent entries.")

    with col_b:
        st.subheader("Under-studied Regions (<5 entries)")
        understudied = results["research_gaps"]["understudied_regions"]
        st.write(", ".join(understudied[:30]) if understudied else "None identified.")

    st.subheader("Gap Score Detail Table")
    gap_tbl = results["research_gaps"]["plastic_gaps"][
        ["plastic", "n_organisms", "n_entries",
         "pct_with_sequence", "pct_extrapolated", "last_year", "gap_score"]
    ].copy()
    gap_tbl["pct_with_sequence"] = (gap_tbl["pct_with_sequence"] * 100).round(1)
    gap_tbl["pct_extrapolated"] = (gap_tbl["pct_extrapolated"].fillna(0) * 100).round(1)
    gap_tbl["gap_score"] = gap_tbl["gap_score"].round(2)
    st.dataframe(gap_tbl, use_container_width=True, hide_index=True)


elif page == "💡 Novel Discovery":
    st.title("Novel Species Discovery Pipeline")
    st.markdown(
        "This pipeline scores each organism for **novelty potential** — combining "
        "plastic breadth, taxonomic rarity, recency, and molecular evidence gaps. "
        "It also identifies phylogenetic gaps and under-characterised environments."
    )

    st.plotly_chart(plot_novelty_scatter(results["novelty_scores"], backend="plotly"),
                    use_container_width=True)

    col_a, col_b = st.columns(2)
    with col_a:
        st.subheader("Top 20 Organisms by Novelty Potential")
        top_novel = results["novelty_scores"].head(20)[
            ["organism", "n_plastics", "novelty_potential",
             "breadth_score", "rarity_score", "recency_score",
             "evidence_gap_score", "has_sequence", "has_enzyme"]
        ].copy()
        top_novel["plastics"] = results["novelty_scores"].head(20)["plastics_degraded"].apply(
            lambda x: ", ".join(x) if isinstance(x, list) else x
        )
        top_novel["novelty_potential"] = top_novel["novelty_potential"].round(2)
        st.dataframe(top_novel, use_container_width=True, hide_index=True)

    with col_b:
        st.subheader("Phylogenetic Gaps — Singleton Genera")
        phylo = identify_phylogenetic_gaps(organisms, df)
        st.dataframe(
            phylo.head(20)[["genus", "n_species", "n_plastics",
                            "is_singleton", "known_genus", "discovery_priority"]],
            use_container_width=True, hide_index=True,
        )

    st.subheader("Underexplored Isolation Environments")
    env_gaps = underexplored_environments(df)
    st.dataframe(
        env_gaps.head(15)[["isolation_environment", "n_species", "n_entries",
                           "pct_with_sequence", "pct_with_enzyme",
                           "characterisation_gap", "exploration_score"]].round(3),
        use_container_width=True, hide_index=True,
    )

    st.subheader("Priority Candidates for Hard-to-Degrade Plastics")
    from src.novel_discovery import PRIORITY_PLASTICS_FOR_DISCOVERY, plastic_specific_candidates
    selected_hard = st.selectbox("Select plastic", PRIORITY_PLASTICS_FOR_DISCOVERY)
    candidates = plastic_specific_candidates(df, selected_hard)
    if not candidates.empty:
        st.dataframe(candidates, use_container_width=True, hide_index=True)
    else:
        st.info(f"No entries found for {selected_hard}.")


elif page == "🔗 Cross-Database":
    st.title("Cross-Database Comparison: PlasticDB vs PAZy")
    st.markdown(
        "**PlasticDB** captures all reported plastic-degrading microorganisms. "
        "**PAZy** focuses only on *thoroughly biochemically characterised* enzymes. "
        "The gap between them represents the research-to-characterisation frontier."
    )

    cdb = results["cross_db"]
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("PlasticDB Species", cdb["plasticdb_n_organisms"])
    col2.metric("PAZy Species", cdb["pazy_n_organisms"])
    col3.metric("PlasticDB Plastic Types", cdb["plasticdb_n_plastics"])
    col4.metric("Shared Plastic Types", len(cdb["shared_plastics"]))

    st.plotly_chart(plot_pazy_vs_plasticdb(cdb, backend="plotly"),
                    use_container_width=True)

    col_a, col_b, col_c = st.columns(3)
    with col_a:
        st.subheader("Shared Plastics")
        st.write(", ".join(sorted(cdb["shared_plastics"])) or "None")
    with col_b:
        st.subheader("PlasticDB Only")
        st.write(", ".join(sorted(cdb["plasticdb_only_plastics"])) or "None")
    with col_c:
        st.subheader("PAZy Only")
        st.write(", ".join(sorted(cdb["pazy_only_plastics"])) or "None")

    st.subheader("PAZy Characterised Enzymes")
    if not pazy.empty:
        st.dataframe(pazy, use_container_width=True, hide_index=True)


elif page == "📋 Data Explorer":
    st.title("Data Explorer")

    st.sidebar.subheader("Filters")
    all_plastics = sorted(df["plastic"].dropna().unique())
    sel_plastics = st.sidebar.multiselect("Plastic type", all_plastics, default=all_plastics[:5])
    year_range = st.sidebar.slider(
        "Year range",
        int(df["year"].min()), int(df["year"].max()),
        (2010, int(df["year"].max())),
    )
    only_with_seq = st.sidebar.checkbox("Only entries with protein sequence", False)
    only_with_enzyme = st.sidebar.checkbox("Only entries with enzyme info", False)

    filtered = df_scored[
        (df_scored["plastic"].isin(sel_plastics)) &
        (df_scored["year"].fillna(0).between(year_range[0], year_range[1]))
    ]
    if only_with_seq:
        filtered = filtered[filtered["has_sequence"]]
    if only_with_enzyme:
        filtered = filtered[filtered["has_enzyme"]]

    st.metric("Filtered rows", f"{len(filtered):,}")

    display_cols = [
        "organism", "genus", "plastic", "year", "evidence_tier",
        "evidence_score", "has_sequence", "has_enzyme",
        "isolation_environment", "isolation_location", "doi",
    ]
    st.dataframe(filtered[display_cols].reset_index(drop=True),
                 use_container_width=True, height=500)

    csv = filtered.drop(columns=["sequence"], errors="ignore").to_csv(index=False)
    st.download_button(
        "Download filtered data (CSV)",
        data=csv,
        file_name="plasticdb_filtered.csv",
        mime="text/csv",
    )
