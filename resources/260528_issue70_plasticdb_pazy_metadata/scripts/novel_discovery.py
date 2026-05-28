"""
Novel species discovery pipeline.

Combines novelty potential scoring with phylogenetic gap detection and
environment-based prioritization to suggest the most promising candidates
for experimental follow-up.
"""

from __future__ import annotations
import pandas as pd
import numpy as np
from collections import defaultdict


WELL_STUDIED_GENERA = {
    "Pseudomonas", "Bacillus", "Aspergillus", "Trichoderma",
    "Streptomyces", "Penicillium", "Ralstonia", "Ideonella",
    "Thermobifida", "Fusarium", "Rhodotorula", "Alcaligenes",
}

PRIORITY_PLASTICS_FOR_DISCOVERY = [
    "PE", "LDPE", "HDPE", "PP", "PS", "PVC",
]

POORLY_CHARACTERISED_PLASTICS = ["LDPE", "PE", "HDPE", "PP", "PS", "PVC"]


def identify_phylogenetic_gaps(organisms_df: pd.DataFrame, plasticdb_df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify genera represented by a single species (singletons) in the DB —
    these are likely vastly under-sampled clades worth deeper exploration.
    """
    genus_species = (
        plasticdb_df.dropna(subset=["genus", "organism"])
        .groupby("genus")["organism"]
        .nunique()
        .reset_index()
        .rename(columns={"organism": "n_species"})
    )
    genus_plastics = (
        plasticdb_df.dropna(subset=["genus"])
        .groupby("genus")["plastic"]
        .nunique()
        .reset_index()
        .rename(columns={"plastic": "n_plastics"})
    )
    merged = genus_species.merge(genus_plastics, on="genus")
    merged["is_singleton"] = merged["n_species"] == 1
    merged["known_genus"] = merged["genus"].isin(WELL_STUDIED_GENERA)
    merged["discovery_priority"] = (
        merged["is_singleton"].astype(int) * 40
        + (~merged["known_genus"]).astype(int) * 35
        + np.log1p(merged["n_plastics"]) * 10
    )
    return merged.sort_values("discovery_priority", ascending=False).reset_index(drop=True)


def underexplored_environments(plasticdb_df: pd.DataFrame) -> pd.DataFrame:
    """
    Find isolation environments with high species diversity but low molecular
    characterisation — these are sampling hotspots for novel degraders.
    """
    env = (
        plasticdb_df.dropna(subset=["isolation_environment"])
        .query("isolation_environment != ''")
        .groupby("isolation_environment")
        .agg(
            n_species=("organism", "nunique"),
            n_entries=("isolation_environment", "count"),
            pct_with_sequence=("has_sequence", "mean"),
            pct_with_enzyme=("has_enzyme", "mean"),
        )
        .reset_index()
    )
    env["characterisation_gap"] = 1 - (env["pct_with_sequence"] + env["pct_with_enzyme"]) / 2
    env["exploration_score"] = (
        np.log1p(env["n_species"]) * 30
        + env["characterisation_gap"] * 50
        + (env["n_entries"] < 10).astype(float) * 20
    )
    return env.sort_values("exploration_score", ascending=False).reset_index(drop=True)


def plastic_specific_candidates(plasticdb_df: pd.DataFrame, plastic: str) -> pd.DataFrame:
    """
    For a given plastic type, rank organisms by their research priority.
    Prioritises: recently reported, not well studied, have isolation metadata.
    """
    sub = plasticdb_df[plasticdb_df["plastic"] == plastic].copy()
    if sub.empty:
        return pd.DataFrame()

    summary = (
        sub.groupby("organism")
        .agg(
            n_entries=("organism", "count"),
            has_sequence=("has_sequence", "any"),
            has_enzyme=("has_enzyme", "any"),
            last_year=("year", "max"),
            isolation_environments=("isolation_environment",
                                    lambda x: "; ".join(x.dropna().unique()[:3])),
            isolation_locations=("isolation_location",
                                 lambda x: "; ".join(x.dropna().unique()[:3])),
        )
        .reset_index()
    )
    max_year = sub["year"].max()
    summary["priority_score"] = (
        (~summary["has_sequence"]).astype(float) * 40
        + (~summary["has_enzyme"]).astype(float) * 30
        + (max_year - summary["last_year"].fillna(max_year - 5)).clip(0) * (-2)
        + summary["n_entries"].apply(lambda n: 20 if n == 1 else 10 if n < 5 else 0)
    )
    return summary.sort_values("priority_score", ascending=False).reset_index(drop=True)


def generate_discovery_report(
    plasticdb_df: pd.DataFrame,
    organisms_df: pd.DataFrame,
    novelty_scores_df: pd.DataFrame,
) -> dict:
    """
    Compile a comprehensive discovery report with top candidates, gap analysis,
    and actionable research directions.
    """
    phylo_gaps = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
    env_gaps = underexplored_environments(plasticdb_df)

    top_novel = novelty_scores_df.head(20)

    hard_plastic_candidates = {}
    for plastic in PRIORITY_PLASTICS_FOR_DISCOVERY:
        candidates = plastic_specific_candidates(plasticdb_df, plastic)
        if not candidates.empty:
            hard_plastic_candidates[plastic] = candidates.head(10)

    undersampled_regions = (
        plasticdb_df.groupby("isolation_location")["organism"]
        .nunique()
        .reset_index()
        .rename(columns={"organism": "n_species"})
        .query("n_species < 5")
        .sort_values("n_species", ascending=False)
    )

    multi_plastic_rare_genera = (
        organisms_df[
            (organisms_df["n_plastics"] >= 3)
            & (~organisms_df["genus"].isin(WELL_STUDIED_GENERA))
        ]
        .sort_values("n_plastics", ascending=False)
        .head(20)
    )

    no_sequence_multi_plastic = (
        organisms_df[
            (organisms_df["n_plastics"] >= 2)
            & (~organisms_df["has_sequence"])
        ]
        .sort_values("n_plastics", ascending=False)
        .head(20)
    )

    return {
        "top_novel_organisms": top_novel,
        "phylogenetic_gaps": phylo_gaps.head(20),
        "underexplored_environments": env_gaps.head(15),
        "hard_plastic_candidates": hard_plastic_candidates,
        "undersampled_regions": undersampled_regions.head(20),
        "multi_plastic_rare_genera": multi_plastic_rare_genera,
        "no_sequence_but_broad_degraders": no_sequence_multi_plastic,
        "summary_stats": {
            "total_singleton_genera": int((phylo_gaps["n_species"] == 1).sum()),
            "n_novel_top20": len(top_novel),
            "n_understudied_environments": len(env_gaps),
            "hard_plastics_covered": len(hard_plastic_candidates),
        },
    }


def format_discovery_report_text(report: dict) -> str:
    """Render discovery report as a human-readable text summary."""
    lines = [
        "=" * 70,
        "  PLASTIC BIODEGRADATION — NOVEL SPECIES DISCOVERY REPORT",
        "=" * 70,
        "",
        "--- TOP 10 ORGANISMS BY NOVELTY POTENTIAL ---",
    ]
    for i, row in report["top_novel_organisms"].head(10).iterrows():
        plastics = ", ".join(row["plastics_degraded"]) if isinstance(row["plastics_degraded"], list) else str(row["plastics_degraded"])
        lines.append(
            f"  {i+1:2}. {row['organism']:<40s}  "
            f"Novelty={row['novelty_potential']:.1f}  "
            f"Plastics={row['n_plastics']}  [{plastics}]"
        )

    lines += [
        "",
        "--- PHYLOGENETIC GAPS (Singleton Genera Worth Exploring) ---",
    ]
    for _, row in report["phylogenetic_gaps"].head(10).iterrows():
        lines.append(
            f"  • {row['genus']:<30s}  species_in_db={row['n_species']}  "
            f"plastics={row['n_plastics']}  "
            f"known={'Yes' if row['known_genus'] else 'No'}"
        )

    lines += [
        "",
        "--- UNDEREXPLORED ISOLATION ENVIRONMENTS ---",
    ]
    for _, row in report["underexplored_environments"].head(8).iterrows():
        lines.append(
            f"  • {row['isolation_environment']:<35s}  "
            f"species={row['n_species']}  "
            f"char_gap={row['characterisation_gap']:.0%}"
        )

    lines += [
        "",
        "--- PRIORITY HARD PLASTICS (PE/PP/PS/PVC/etc.) ---",
    ]
    for plastic, df in report["hard_plastic_candidates"].items():
        lines.append(f"  [{plastic}] Top candidate: {df.iloc[0]['organism']}  "
                     f"(score={df.iloc[0]['priority_score']:.0f})")

    lines += [
        "",
        "--- MULTI-PLASTIC DEGRADERS LACKING SEQUENCE DATA ---",
    ]
    for _, row in report["no_sequence_but_broad_degraders"].head(8).iterrows():
        lines.append(
            f"  • {row['organism']:<40s}  n_plastics={row['n_plastics']}"
        )

    lines += [
        "",
        f"SUMMARY: {report['summary_stats']['total_singleton_genera']} singleton genera, "
        f"{report['summary_stats']['n_understudied_environments']} underexplored environments.",
        "=" * 70,
    ]
    return "\n".join(lines)
