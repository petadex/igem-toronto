"""
Meta-analysis functions for plastic biodegradation research.

Covers:
- Taxonomic diversity
- Plastic substrate coverage
- Geographic & temporal trends
- Evidence quality scoring
- Research gap identification
"""

from __future__ import annotations
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import MultiLabelBinarizer
from collections import Counter


def evidence_quality_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign a composite evidence quality score (0-100) to each entry.

    Components:
      - has_sequence (30 pts)
      - has_genbank  (20 pts)
      - has_enzyme   (20 pts)
      - analytical_grade plastic used (15 pts)
      - not extrapolated from enzyme (15 pts)
    """
    df = df.copy()
    score = pd.Series(0.0, index=df.index)
    score += df["has_sequence"].fillna(False).astype(float) * 30
    score += df["has_genbank"].fillna(False).astype(float) * 20
    score += df["has_enzyme"].fillna(False).astype(float) * 20
    score += df["analytical_grade"].fillna(False).astype(float) * 15
    score += (~df["extrapolated_from_enzyme"].fillna(True)).astype(float) * 15
    df["evidence_score"] = score
    df["evidence_tier"] = pd.cut(
        score,
        bins=[-1, 20, 50, 80, 100],
        labels=["Low", "Medium", "High", "Excellent"],
    )
    return df


def taxonomic_diversity(df: pd.DataFrame) -> dict:
    """
    Compute taxonomic diversity metrics at genus and species level.

    Returns a dict with counts, Shannon diversity index, and top genera.
    """
    genus_counts = df.dropna(subset=["genus"])["genus"].value_counts()
    species_counts = df.dropna(subset=["organism"])["organism"].value_counts()

    def shannon(counts):
        total = counts.sum()
        p = counts / total
        return -(p * np.log(p)).sum()

    return {
        "n_unique_genera": int(genus_counts.shape[0]),
        "n_unique_species": int(species_counts.shape[0]),
        "genus_shannon_diversity": round(float(shannon(genus_counts)), 4),
        "top_10_genera": genus_counts.head(10).to_dict(),
        "top_10_species": species_counts.head(10).to_dict(),
        "singleton_genera": int((genus_counts == 1).sum()),
        "pct_singleton_genera": round(float((genus_counts == 1).mean() * 100), 1),
    }


def plastic_coverage_matrix(organisms_df: pd.DataFrame, top_n_plastics: int = 20) -> pd.DataFrame:
    """
    Build a binary organism × plastic coverage matrix for the top-N plastics.
    Useful for co-occurrence and novelty analysis.
    """
    exploded = organisms_df[["organism", "plastics_degraded"]].explode("plastics_degraded")
    exploded = exploded.rename(columns={"plastics_degraded": "plastic"})
    top_plastics = exploded["plastic"].value_counts().head(top_n_plastics).index.tolist()
    exploded = exploded[exploded["plastic"].isin(top_plastics)]
    matrix = exploded.pivot_table(index="organism", columns="plastic", aggfunc=len, fill_value=0)
    matrix = (matrix > 0).astype(int)
    return matrix


def plastic_co_occurrence(organisms_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute how often pairs of plastics are degraded by the same organism.
    Returns a symmetric co-occurrence count matrix.
    """
    mlb = MultiLabelBinarizer()
    binary = pd.DataFrame(
        mlb.fit_transform(organisms_df["plastics_degraded"]),
        columns=mlb.classes_,
    )
    co = binary.T @ binary
    return co


def temporal_trend_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """
    Year-by-year publication counts with rolling average and growth rate.
    """
    yearly = (
        df.dropna(subset=["year"])
        .groupby("year")
        .agg(
            n_entries=("organism", "count"),
            n_unique_species=("organism", "nunique"),
            n_unique_plastics=("plastic", "nunique"),
            n_with_sequence=("has_sequence", "sum"),
        )
        .reset_index()
        .sort_values("year")
    )
    yearly["rolling_3yr"] = yearly["n_entries"].rolling(3, center=True).mean()
    yearly["yoy_growth_pct"] = yearly["n_entries"].pct_change() * 100
    yearly["cumulative_entries"] = yearly["n_entries"].cumsum()
    yearly["cumulative_species"] = (
        df.dropna(subset=["year"])
        .sort_values("year")
        .groupby("year")["organism"]
        .apply(lambda x: x)
        .groupby(level=0)
        .agg(lambda _: None)  # placeholder; computed below
    )
    cumulative_species = (
        df.dropna(subset=["year"])
        .sort_values("year")
        .drop_duplicates(subset=["organism"])
        .groupby("year")
        .size()
        .cumsum()
        .rename("cumulative_species")
        .reset_index()
    )
    yearly = yearly.drop(columns=["cumulative_species"]).merge(cumulative_species, on="year", how="left")
    yearly["cumulative_species"] = yearly["cumulative_species"].ffill()
    return yearly


def geographic_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise research activity by isolation location.
    """
    geo = (
        df.dropna(subset=["isolation_location"])
        .query("isolation_location != ''")
        .groupby("isolation_location")
        .agg(
            n_entries=("organism", "count"),
            n_species=("organism", "nunique"),
            n_plastics=("plastic", "nunique"),
            years=("year", lambda x: sorted(x.dropna().unique().tolist())),
            most_common_plastic=("plastic", lambda x: x.mode().iloc[0] if len(x) > 0 else None),
        )
        .reset_index()
        .sort_values("n_entries", ascending=False)
    )
    return geo


def isolation_environment_profile(df: pd.DataFrame) -> pd.DataFrame:
    """
    Break down discoveries by isolation environment type (soil, marine, etc.).
    """
    env = (
        df.dropna(subset=["isolation_environment"])
        .query("isolation_environment != ''")
        .groupby("isolation_environment")
        .agg(
            n_entries=("organism", "count"),
            n_species=("organism", "nunique"),
            n_plastics=("plastic", "nunique"),
        )
        .reset_index()
        .sort_values("n_entries", ascending=False)
    )
    return env


def research_gap_analysis(df: pd.DataFrame) -> dict:
    """
    Identify the biggest research gaps:
    - Plastics with few unique degrading organisms
    - Plastics with low sequence/enzyme evidence
    - Plastics with no recent (post-2020) entries
    - Under-studied geographic regions
    """
    plastic_summary = (
        df.groupby("plastic")
        .agg(
            n_organisms=("organism", "nunique"),
            n_entries=("plastic", "count"),
            pct_with_sequence=("has_sequence", "mean"),
            pct_extrapolated=("extrapolated_from_enzyme", "mean"),
            last_year=("year", "max"),
        )
        .reset_index()
    )
    plastic_summary["gap_score"] = (
        (1 / np.log1p(plastic_summary["n_organisms"])) * 40
        + (1 - plastic_summary["pct_with_sequence"]) * 35
        + plastic_summary["pct_extrapolated"].fillna(0) * 25
    )
    plastic_summary = plastic_summary.sort_values("gap_score", ascending=False)

    recent = df[df["year"] >= 2020]
    plastics_with_recent = set(recent["plastic"].dropna().unique())
    all_plastics = set(df["plastic"].dropna().unique())
    no_recent_data = all_plastics - plastics_with_recent

    geo_counts = df["isolation_location"].value_counts()
    understudied_regions = geo_counts[geo_counts < 5].index.tolist()

    top_organisms_per_plastic = (
        df.groupby("plastic")["organism"]
        .apply(lambda x: x.value_counts().head(3).index.tolist())
        .to_dict()
    )

    return {
        "plastic_gaps": plastic_summary,
        "no_recent_data_plastics": sorted(no_recent_data),
        "understudied_regions": understudied_regions[:20],
        "dominant_organisms_per_plastic": top_organisms_per_plastic,
    }


def cross_database_comparison(plasticdb_df: pd.DataFrame, pazy_df: pd.DataFrame) -> dict:
    """
    Compare PlasticDB (broad coverage) vs PAZy (thoroughly characterised) organisms
    and plastics.
    """
    plasticdb_plastics = set(plasticdb_df["plastic"].dropna().unique())
    pazy_plastics = set(pazy_df["plastic"].dropna().unique()) if "plastic" in pazy_df.columns else set()
    plasticdb_only = plasticdb_plastics - pazy_plastics
    pazy_only = pazy_plastics - plasticdb_plastics
    shared = plasticdb_plastics & pazy_plastics

    pdb_orgs = set(plasticdb_df["organism"].dropna().unique())
    pazy_orgs = set(pazy_df["organism"].dropna().unique()) if "organism" in pazy_df.columns else set()

    return {
        "plasticdb_n_plastics": len(plasticdb_plastics),
        "pazy_n_plastics": len(pazy_plastics),
        "shared_plastics": sorted(shared),
        "plasticdb_only_plastics": sorted(plasticdb_only),
        "pazy_only_plastics": sorted(pazy_only),
        "plasticdb_n_organisms": len(pdb_orgs),
        "pazy_n_organisms": len(pazy_orgs),
        "overlap_organisms": sorted(pdb_orgs & pazy_orgs),
        "coverage_ratio": round(len(shared) / len(plasticdb_plastics | pazy_plastics), 3),
    }


def compute_novelty_potential(organisms_df: pd.DataFrame, plasticdb_df: pd.DataFrame) -> pd.DataFrame:
    """
    Score each organism for potential novelty / research interest:
    - breadth: number of different plastics it degrades
    - rarity: how rare is the genus compared to others in DB
    - recency: how recently it was reported
    - environment diversity: different environments sampled from
    - evidence gap: high breadth but low molecular evidence = priority target
    """
    genus_freq = plasticdb_df["genus"].value_counts().to_dict()
    total_entries = len(plasticdb_df)
    max_plastics = organisms_df["n_plastics"].max()
    max_year = plasticdb_df["year"].max()

    records = []
    for _, row in organisms_df.iterrows():
        g = row.get("genus", "")
        freq = genus_freq.get(g, 1)
        rarity_score = max(0, 100 - (freq / total_entries * 1000))
        breadth_score = (row["n_plastics"] / max_plastics) * 100
        recency = row.get("last_year", 2000) or 2000
        recency_score = max(0, 100 - (max_year - recency) * 5)
        has_seq = row.get("has_sequence", False)
        has_enz = row.get("has_enzyme", False)
        evidence_gap = (100 - (50 * has_seq + 50 * has_enz)) if (breadth_score > 30) else 0

        novelty = (
            breadth_score * 0.30
            + rarity_score * 0.25
            + recency_score * 0.20
            + evidence_gap * 0.25
        )
        records.append({
            "organism": row["organism"],
            "genus": g,
            "n_plastics": row["n_plastics"],
            "plastics_degraded": row["plastics_degraded"],
            "breadth_score": round(breadth_score, 2),
            "rarity_score": round(rarity_score, 2),
            "recency_score": round(recency_score, 2),
            "evidence_gap_score": round(evidence_gap, 2),
            "novelty_potential": round(novelty, 2),
            "last_year": recency,
            "has_sequence": has_seq,
            "has_enzyme": has_enz,
        })

    result = pd.DataFrame(records).sort_values("novelty_potential", ascending=False)
    return result.reset_index(drop=True)


def run_full_analysis(data: dict) -> dict:
    """Run all analyses and return a results dict."""
    df = data["plasticdb"]
    organisms = data["organisms"]
    pazy = data["pazy"]

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
