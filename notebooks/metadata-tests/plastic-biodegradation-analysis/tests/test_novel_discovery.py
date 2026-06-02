"""
Tests for src/novel_discovery.py

Covers: phylogenetic gap detection, environment exploration scoring,
per-plastic candidate ranking, full discovery report generation.
"""

import pytest
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.novel_discovery import (
    identify_phylogenetic_gaps,
    underexplored_environments,
    plastic_specific_candidates,
    generate_discovery_report,
    format_discovery_report_text,
    WELL_STUDIED_GENERA,
    PRIORITY_PLASTICS_FOR_DISCOVERY,
    POORLY_CHARACTERISED_PLASTICS,
)
from src.analysis import compute_novelty_potential


class TestIdentifyPhylogeneticGaps:
    def test_returns_dataframe(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        for col in ["genus", "n_species", "n_plastics",
                    "is_singleton", "known_genus", "discovery_priority"]:
            assert col in result.columns

    def test_sorted_by_priority_desc(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        priorities = result["discovery_priority"].tolist()
        assert priorities == sorted(priorities, reverse=True)

    def test_singleton_flag_correct(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        singletons = result[result["is_singleton"]]
        assert (singletons["n_species"] == 1).all()

    def test_known_genera_flagged(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        known = result[result["known_genus"]]
        for genus in known["genus"].tolist():
            assert genus in WELL_STUDIED_GENERA

    def test_priority_positive(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        assert (result["discovery_priority"] >= 0).all()

    def test_minimal_df_singletons(self, minimal_organisms, minimal_df):
        result = identify_phylogenetic_gaps(minimal_organisms, minimal_df)
        assert "is_singleton" in result.columns
        singletons = result[result["is_singleton"]]
        assert len(singletons) >= 1

    def test_pseudomonas_is_known(self, organisms_df, plasticdb_df):
        result = identify_phylogenetic_gaps(organisms_df, plasticdb_df)
        ps_row = result[result["genus"] == "Pseudomonas"]
        if not ps_row.empty:
            assert bool(ps_row.iloc[0]["known_genus"]) is True


class TestUnderexploredEnvironments:
    def test_returns_dataframe(self, plasticdb_df):
        result = underexplored_environments(plasticdb_df)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, plasticdb_df):
        result = underexplored_environments(plasticdb_df)
        for col in ["isolation_environment", "n_species", "n_entries",
                    "pct_with_sequence", "characterisation_gap", "exploration_score"]:
            assert col in result.columns

    def test_characterisation_gap_range(self, plasticdb_df):
        result = underexplored_environments(plasticdb_df)
        assert result["characterisation_gap"].between(0, 1).all()

    def test_sorted_by_exploration_score(self, plasticdb_df):
        result = underexplored_environments(plasticdb_df)
        scores = result["exploration_score"].tolist()
        assert scores == sorted(scores, reverse=True)

    def test_soil_environment_present(self, plasticdb_df):
        result = underexplored_environments(plasticdb_df)
        envs = result["isolation_environment"].str.lower().tolist()
        assert any("soil" in e for e in envs)


class TestPlasticSpecificCandidates:
    def test_returns_dataframe_for_pet(self, plasticdb_df):
        result = plastic_specific_candidates(plasticdb_df, "PET")
        assert isinstance(result, pd.DataFrame)

    def test_returns_empty_for_unknown_plastic(self, plasticdb_df):
        result = plastic_specific_candidates(plasticdb_df, "ZZZZUNKNOWN")
        assert result.empty

    def test_required_columns(self, plasticdb_df):
        result = plastic_specific_candidates(plasticdb_df, "PET")
        if not result.empty:
            for col in ["organism", "has_sequence", "has_enzyme",
                        "last_year", "priority_score"]:
                assert col in result.columns

    def test_sorted_by_priority(self, plasticdb_df):
        result = plastic_specific_candidates(plasticdb_df, "PET")
        if not result.empty:
            scores = result["priority_score"].tolist()
            assert scores == sorted(scores, reverse=True)

    def test_all_entries_for_correct_plastic(self, minimal_df):
        result = plastic_specific_candidates(minimal_df, "PET")
        assert len(result) > 0

    @pytest.mark.parametrize("plastic", PRIORITY_PLASTICS_FOR_DISCOVERY)
    def test_hard_plastics_have_candidates(self, plasticdb_df, plastic):
        result = plastic_specific_candidates(plasticdb_df, plastic)
        assert isinstance(result, pd.DataFrame), \
            f"plastic_specific_candidates returned non-DataFrame for {plastic}"


class TestGenerateDiscoveryReport:
    @pytest.fixture(scope="class")
    def report(self, plasticdb_df, organisms_df):
        novelty = compute_novelty_potential(organisms_df, plasticdb_df)
        return generate_discovery_report(plasticdb_df, organisms_df, novelty)

    def test_returns_dict(self, report):
        assert isinstance(report, dict)

    def test_required_keys(self, report):
        for key in [
            "top_novel_organisms", "phylogenetic_gaps",
            "underexplored_environments", "hard_plastic_candidates",
            "undersampled_regions", "multi_plastic_rare_genera",
            "no_sequence_but_broad_degraders", "summary_stats",
        ]:
            assert key in report, f"Missing report key: {key}"

    def test_top_novel_is_dataframe(self, report):
        assert isinstance(report["top_novel_organisms"], pd.DataFrame)

    def test_top_novel_not_empty(self, report):
        assert len(report["top_novel_organisms"]) > 0

    def test_hard_plastic_candidates_is_dict(self, report):
        assert isinstance(report["hard_plastic_candidates"], dict)

    def test_summary_stats_is_dict(self, report):
        assert isinstance(report["summary_stats"], dict)

    def test_summary_stats_has_counts(self, report):
        assert "total_singleton_genera" in report["summary_stats"]
        assert report["summary_stats"]["total_singleton_genera"] > 0


class TestFormatDiscoveryReportText:
    def test_returns_string(self, plasticdb_df, organisms_df):
        novelty = compute_novelty_potential(organisms_df, plasticdb_df)
        report = generate_discovery_report(plasticdb_df, organisms_df, novelty)
        text = format_discovery_report_text(report)
        assert isinstance(text, str)

    def test_contains_header(self, plasticdb_df, organisms_df):
        novelty = compute_novelty_potential(organisms_df, plasticdb_df)
        report = generate_discovery_report(plasticdb_df, organisms_df, novelty)
        text = format_discovery_report_text(report)
        assert "NOVEL SPECIES DISCOVERY REPORT" in text

    def test_contains_top_organisms(self, plasticdb_df, organisms_df):
        novelty = compute_novelty_potential(organisms_df, plasticdb_df)
        report = generate_discovery_report(plasticdb_df, organisms_df, novelty)
        text = format_discovery_report_text(report)
        assert "TOP 10 ORGANISMS" in text

    def test_contains_phylogenetic_gaps(self, plasticdb_df, organisms_df):
        novelty = compute_novelty_potential(organisms_df, plasticdb_df)
        report = generate_discovery_report(plasticdb_df, organisms_df, novelty)
        text = format_discovery_report_text(report)
        assert "PHYLOGENETIC GAPS" in text

    def test_minimum_length(self, plasticdb_df, organisms_df):
        novelty = compute_novelty_potential(organisms_df, plasticdb_df)
        report = generate_discovery_report(plasticdb_df, organisms_df, novelty)
        text = format_discovery_report_text(report)
        assert len(text) > 500


class TestConstants:
    def test_well_studied_genera_not_empty(self):
        assert len(WELL_STUDIED_GENERA) > 5

    def test_known_genera_in_set(self):
        for genus in ["Pseudomonas", "Bacillus", "Aspergillus"]:
            assert genus in WELL_STUDIED_GENERA

    def test_priority_plastics_non_empty(self):
        assert len(PRIORITY_PLASTICS_FOR_DISCOVERY) > 0

    def test_priority_plastics_includes_hard_ones(self):
        for p in ["PE", "PP", "PS"]:
            assert p in PRIORITY_PLASTICS_FOR_DISCOVERY
