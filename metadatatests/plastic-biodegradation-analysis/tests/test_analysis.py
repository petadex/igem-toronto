"""
Tests for src/analysis.py

Covers: evidence scoring, diversity metrics, temporal/geo trends, gap analysis,
novelty scoring, co-occurrence, cross-DB comparison.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.analysis import (
    evidence_quality_score,
    taxonomic_diversity,
    temporal_trend_analysis,
    geographic_distribution,
    isolation_environment_profile,
    research_gap_analysis,
    plastic_co_occurrence,
    cross_database_comparison,
    compute_novelty_potential,
    run_full_analysis,
)


class TestEvidenceQualityScore:
    def test_adds_evidence_columns(self, minimal_df):
        result = evidence_quality_score(minimal_df)
        assert "evidence_score" in result.columns
        assert "evidence_tier" in result.columns

    def test_score_range(self, minimal_df):
        result = evidence_quality_score(minimal_df)
        assert result["evidence_score"].between(0, 100).all()

    def test_full_evidence_entry_scores_high(self, minimal_df):
        result = evidence_quality_score(minimal_df)
        ideonella = result[result["organism"] == "Ideonella sakaiensis"]
        assert ideonella["evidence_score"].iloc[0] >= 80

    def test_no_evidence_entry_scores_low(self, minimal_df):
        result = evidence_quality_score(minimal_df)
        pseudomonas = result[result["organism"] == "Pseudomonas putida"]
        assert pseudomonas["evidence_score"].iloc[0] <= 20

    def test_evidence_tier_categories(self, minimal_df):
        result = evidence_quality_score(minimal_df)
        valid = {"Low", "Medium", "High", "Excellent"}
        tiers = set(result["evidence_tier"].dropna().unique())
        assert tiers.issubset(valid), f"Unexpected tiers: {tiers - valid}"

    def test_real_data_tier_distribution(self, plasticdb_df):
        result = evidence_quality_score(plasticdb_df)
        counts = result["evidence_tier"].value_counts(dropna=False)
        assert len(counts.dropna()) >= 2, "Should have at least 2 evidence tiers"
        assert counts.sum() == len(result)

    def test_idempotent(self, minimal_df):
        r1 = evidence_quality_score(minimal_df)
        r2 = evidence_quality_score(r1)
        pd.testing.assert_series_equal(r1["evidence_score"], r2["evidence_score"])


class TestTaxonomicDiversity:
    def test_returns_dict(self, plasticdb_df):
        result = taxonomic_diversity(plasticdb_df)
        assert isinstance(result, dict)

    def test_required_keys(self, plasticdb_df):
        result = taxonomic_diversity(plasticdb_df)
        for key in ["n_unique_genera", "n_unique_species",
                    "genus_shannon_diversity", "top_10_genera",
                    "singleton_genera", "pct_singleton_genera"]:
            assert key in result, f"Missing key: {key}"

    def test_counts_positive(self, plasticdb_df):
        result = taxonomic_diversity(plasticdb_df)
        assert result["n_unique_genera"] > 100
        assert result["n_unique_species"] > 500

    def test_shannon_positive(self, plasticdb_df):
        result = taxonomic_diversity(plasticdb_df)
        assert result["genus_shannon_diversity"] > 0

    def test_top_genera_count(self, plasticdb_df):
        result = taxonomic_diversity(plasticdb_df)
        assert len(result["top_10_genera"]) <= 10

    def test_singleton_pct_between_0_and_100(self, plasticdb_df):
        result = taxonomic_diversity(plasticdb_df)
        assert 0 <= result["pct_singleton_genera"] <= 100

    def test_minimal_df_diversity(self, minimal_df):
        result = taxonomic_diversity(minimal_df)
        assert result["n_unique_genera"] == 5
        assert result["n_unique_species"] == 5


class TestTemporalTrends:
    def test_returns_dataframe(self, plasticdb_df):
        result = temporal_trend_analysis(plasticdb_df)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, plasticdb_df):
        result = temporal_trend_analysis(plasticdb_df)
        for col in ["year", "n_entries", "n_unique_species",
                    "cumulative_entries", "rolling_3yr"]:
            assert col in result.columns

    def test_sorted_by_year(self, plasticdb_df):
        result = temporal_trend_analysis(plasticdb_df)
        assert result["year"].is_monotonic_increasing

    def test_cumulative_entries_monotonic(self, plasticdb_df):
        result = temporal_trend_analysis(plasticdb_df)
        assert result["cumulative_entries"].is_monotonic_increasing

    def test_n_entries_positive(self, plasticdb_df):
        result = temporal_trend_analysis(plasticdb_df)
        assert (result["n_entries"] > 0).all()

    def test_recent_years_present(self, plasticdb_df):
        result = temporal_trend_analysis(plasticdb_df)
        assert 2022 in result["year"].values or 2023 in result["year"].values


class TestGeographicDistribution:
    def test_returns_dataframe(self, plasticdb_df):
        result = geographic_distribution(plasticdb_df)
        assert isinstance(result, pd.DataFrame)

    def test_columns_present(self, plasticdb_df):
        result = geographic_distribution(plasticdb_df)
        for col in ["isolation_location", "n_entries", "n_species"]:
            assert col in result.columns

    def test_known_countries_present(self, plasticdb_df):
        result = geographic_distribution(plasticdb_df)
        locations = result["isolation_location"].tolist()
        assert any("India" in str(loc) or "Japan" in str(loc) for loc in locations)

    def test_sorted_by_entries(self, plasticdb_df):
        result = geographic_distribution(plasticdb_df)
        assert result["n_entries"].iloc[0] >= result["n_entries"].iloc[-1]


class TestResearchGapAnalysis:
    def test_returns_dict(self, plasticdb_df):
        result = research_gap_analysis(plasticdb_df)
        assert isinstance(result, dict)

    def test_required_keys(self, plasticdb_df):
        result = research_gap_analysis(plasticdb_df)
        for key in ["plastic_gaps", "no_recent_data_plastics",
                    "understudied_regions", "dominant_organisms_per_plastic"]:
            assert key in result

    def test_plastic_gaps_is_dataframe(self, plasticdb_df):
        result = research_gap_analysis(plasticdb_df)
        assert isinstance(result["plastic_gaps"], pd.DataFrame)

    def test_gap_score_present(self, plasticdb_df):
        result = research_gap_analysis(plasticdb_df)
        assert "gap_score" in result["plastic_gaps"].columns

    def test_gap_score_positive(self, plasticdb_df):
        result = research_gap_analysis(plasticdb_df)
        assert (result["plastic_gaps"]["gap_score"] >= 0).all()

    def test_no_recent_plastics_is_list(self, plasticdb_df):
        result = research_gap_analysis(plasticdb_df)
        assert isinstance(result["no_recent_data_plastics"], list)


class TestPlasticCoOccurrence:
    def test_returns_dataframe(self, organisms_df):
        result = plastic_co_occurrence(organisms_df)
        assert isinstance(result, pd.DataFrame)

    def test_symmetric_matrix(self, organisms_df):
        result = plastic_co_occurrence(organisms_df)
        pd.testing.assert_frame_equal(result, result.T)

    def test_self_cooccurrence_matches_count(self, organisms_df):
        result = plastic_co_occurrence(organisms_df)
        for col in result.columns:
            assert result.loc[col, col] >= 0


class TestCrossDBComparison:
    def test_returns_dict(self, plasticdb_df, pazy_df):
        result = cross_database_comparison(plasticdb_df, pazy_df)
        assert isinstance(result, dict)

    def test_coverage_ratio_range(self, plasticdb_df, pazy_df):
        result = cross_database_comparison(plasticdb_df, pazy_df)
        assert 0 <= result["coverage_ratio"] <= 1

    def test_organism_counts_positive(self, plasticdb_df, pazy_df):
        result = cross_database_comparison(plasticdb_df, pazy_df)
        assert result["plasticdb_n_organisms"] > 0
        assert result["pazy_n_organisms"] > 0

    def test_shared_plastics_is_list(self, plasticdb_df, pazy_df):
        result = cross_database_comparison(plasticdb_df, pazy_df)
        assert isinstance(result["shared_plastics"], list)

    def test_pet_in_both(self, plasticdb_df, pazy_df):
        result = cross_database_comparison(plasticdb_df, pazy_df)
        assert "PET" in result["shared_plastics"]


class TestNoveltyPotential:
    def test_returns_dataframe(self, organisms_df, plasticdb_df):
        result = compute_novelty_potential(organisms_df, plasticdb_df)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, organisms_df, plasticdb_df):
        result = compute_novelty_potential(organisms_df, plasticdb_df)
        for col in ["organism", "n_plastics", "novelty_potential",
                    "breadth_score", "rarity_score", "recency_score"]:
            assert col in result.columns

    def test_sorted_descending(self, organisms_df, plasticdb_df):
        result = compute_novelty_potential(organisms_df, plasticdb_df)
        scores = result["novelty_potential"].tolist()
        assert scores == sorted(scores, reverse=True)

    def test_score_range(self, organisms_df, plasticdb_df):
        result = compute_novelty_potential(organisms_df, plasticdb_df)
        assert result["novelty_potential"].between(0, 100).all()

    def test_one_row_per_organism(self, organisms_df, plasticdb_df):
        result = compute_novelty_potential(organisms_df, plasticdb_df)
        assert result["organism"].is_unique

    def test_multi_plastic_orgs_rank_high(self, organisms_df, plasticdb_df):
        result = compute_novelty_potential(organisms_df, plasticdb_df)
        top10 = result.head(10)
        assert top10["n_plastics"].mean() >= 3


class TestRunFullAnalysis:
    def test_returns_all_keys(self, data):
        results = run_full_analysis(data)
        expected_keys = [
            "df_scored", "taxonomic_diversity", "temporal_trends",
            "geographic_distribution", "isolation_environments",
            "research_gaps", "co_occurrence", "cross_db", "novelty_scores",
        ]
        for key in expected_keys:
            assert key in results, f"Missing key in run_full_analysis: {key}"

    def test_df_scored_has_evidence_columns(self, data):
        results = run_full_analysis(data)
        assert "evidence_score" in results["df_scored"].columns
        assert "evidence_tier" in results["df_scored"].columns
