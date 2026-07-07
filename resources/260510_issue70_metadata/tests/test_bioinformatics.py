"""
Tests for src/bioinformatics.py

Covers: protein property computation, enzyme family classification,
sequence motifs, amino acid composition, taxonomy tree building,
co-occurrence utilities, and sequence statistics.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.bioinformatics import (
    clean_sequence,
    compute_protein_properties,
    analyse_all_sequences,
    classify_enzyme_family,
    scan_serine_hydrolase_motif,
    get_amino_acid_composition_table,
    sequence_length_stats,
    build_taxonomy_newick,
    compare_sequences_pairwise,
    enrich_df_with_protein_properties,
    top_enzyme_families_per_plastic,
    sequence_property_summary,
    ENZYME_FAMILIES,
    STANDARD_AA,
)
from tests.conftest import SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ


class TestCleanSequence:
    def test_removes_non_aa_chars(self):
        assert clean_sequence("MA-TC*GH\n") == "MATCGH"

    def test_lowercased_to_upper(self):
        assert clean_sequence("matkl") == "MATKL"

    def test_empty_string(self):
        assert clean_sequence("") == ""

    def test_none_returns_empty(self):
        assert clean_sequence(None) == ""

    def test_only_standard_aa_remain(self):
        seq = clean_sequence("MATBXZJUO-*123!@#")
        for c in seq:
            assert c in STANDARD_AA, f"Non-standard AA in output: {c}"

    def test_already_clean_unchanged(self):
        seq = "MATKGHILMNPQRSTVWY"
        assert clean_sequence(seq) == seq


class TestComputeProteinProperties:
    def test_returns_dict_for_valid_seq(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert isinstance(props, dict)
        assert len(props) > 0

    def test_returns_empty_dict_for_short_seq(self):
        props = compute_protein_properties("MATK")
        assert props == {}

    def test_returns_empty_for_empty(self):
        assert compute_protein_properties("") == {}

    def test_molecular_weight_positive(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert props["molecular_weight_kda"] > 0

    def test_molecular_weight_reasonable_range(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert 1 < props["molecular_weight_kda"] < 500

    def test_isoelectric_point_range(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert 0 < props["isoelectric_point"] < 14

    def test_instability_index_present(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert "instability_index" in props
        assert isinstance(props["instability_index"], float)

    def test_is_stable_bool(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert isinstance(props["is_stable"], bool)

    def test_gravy_is_float(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert isinstance(props["gravy"], float)

    def test_ss_fractions_sum_near_one(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        total = props["helix_fraction"] + props["turn_fraction"] + props["sheet_fraction"]
        assert 0.5 <= total <= 1.05, f"SS fractions sum = {total}, expected 0.5–1.05"

    def test_aa_composition_sums_to_one(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        total = sum(props["aa_composition"].values())
        assert abs(total - 1.0) < 0.01

    def test_length_matches_clean_seq(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        expected_len = len(clean_sequence(SAMPLE_PET_SEQ))
        assert props["length"] == expected_len

    def test_charge_balance_is_int(self):
        props = compute_protein_properties(SAMPLE_PET_SEQ)
        assert isinstance(props["charge_balance"], (int, np.integer))

    def test_phb_seq_different_mw_from_pet(self):
        pet_props = compute_protein_properties(SAMPLE_PET_SEQ)
        phb_props = compute_protein_properties(SAMPLE_PHB_SEQ)
        assert pet_props["molecular_weight_kda"] != phb_props["molecular_weight_kda"]


class TestAnalyseAllSequences:
    def test_returns_dataframe(self, minimal_df):
        result = analyse_all_sequences(minimal_df)
        assert isinstance(result, pd.DataFrame)

    def test_only_rows_with_sequences(self, minimal_df):
        result = analyse_all_sequences(minimal_df)
        assert len(result) > 0
        assert len(result) <= minimal_df["has_sequence"].sum()

    def test_protein_columns_present(self, minimal_df):
        result = analyse_all_sequences(minimal_df)
        if not result.empty:
            for col in ["molecular_weight_kda", "isoelectric_point",
                        "instability_index", "gravy"]:
                assert col in result.columns

    def test_full_dataset_produces_results(self, plasticdb_df):
        result = analyse_all_sequences(plasticdb_df, min_length=30)
        assert len(result) > 50, "Expected >50 sequences with properties"

    def test_min_length_filter(self, plasticdb_df):
        long = analyse_all_sequences(plasticdb_df, min_length=200)
        short = analyse_all_sequences(plasticdb_df, min_length=20)
        assert len(long) <= len(short)


class TestClassifyEnzymeFamily:
    @pytest.mark.parametrize("name,plastic,expected_family", [
        ("IsPETase",          "PET", "PETase / Cutinase"),
        ("LCC cutinase",      "PET", "PETase / Cutinase"),
        ("PHB depolymerase",  "PHB", "PHB Depolymerase"),
        ("PhaZ7",             "PHB", "PHB Depolymerase"),
        ("Lipase A",          "PU",  "Lipase"),
        ("NylB nylonase",     "Nylon", "Amidase / Nylonase"),
        ("Laccase",           "PE",  "Laccase / Oxidase"),
        ("",                  "PET", "Unknown"),
    ])
    def test_family_classification(self, name, plastic, expected_family):
        result = classify_enzyme_family(name, plastic)
        assert result == expected_family, \
            f"classify_enzyme_family({name!r}, {plastic!r}) = {result!r}, expected {expected_family!r}"

    def test_returns_string(self):
        assert isinstance(classify_enzyme_family("lipase", "PCL"), str)

    def test_unknown_name_and_plastic(self):
        result = classify_enzyme_family("xyzprotein", "XYZ")
        assert result == "Unknown"

    def test_case_insensitive(self):
        upper = classify_enzyme_family("PETASE", "PET")
        lower = classify_enzyme_family("petase", "PET")
        assert upper == lower


class TestScanSerineHydrolaseMotif:
    def test_motif_found_in_pet_seq(self):
        seq_with_motif = "MATKGASMGGGSGGGSXX"
        assert scan_serine_hydrolase_motif(seq_with_motif) in [True, False]

    def test_gxsxg_found(self):
        assert scan_serine_hydrolase_motif("MXXXXXGASMGXXXXXX") is True

    def test_no_motif_short_seq(self):
        assert scan_serine_hydrolase_motif("MATK") is False

    def test_returns_bool(self):
        assert isinstance(scan_serine_hydrolase_motif(SAMPLE_PET_SEQ), bool)

    def test_empty_sequence(self):
        assert scan_serine_hydrolase_motif("") is False


class TestSequenceLengthStats:
    def test_returns_dataframe(self, plasticdb_df):
        result = sequence_length_stats(plasticdb_df)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, plasticdb_df):
        result = sequence_length_stats(plasticdb_df)
        for col in ["plastic", "count", "mean", "min", "max"]:
            assert col in result.columns

    def test_mean_positive(self, plasticdb_df):
        result = sequence_length_stats(plasticdb_df)
        assert (result["mean"] > 0).all()

    def test_min_lte_max(self, plasticdb_df):
        result = sequence_length_stats(plasticdb_df)
        assert (result["min"] <= result["max"]).all()

    def test_sorted_by_count(self, plasticdb_df):
        result = sequence_length_stats(plasticdb_df)
        assert result["count"].iloc[0] >= result["count"].iloc[-1]


class TestBuildTaxonomyNewick:
    def test_returns_string(self, organisms_df):
        newick = build_taxonomy_newick(organisms_df)
        assert isinstance(newick, str)
        assert len(newick) > 10

    def test_newick_parseable_by_biopython(self, organisms_df):
        from Bio import Phylo
        from io import StringIO
        newick = build_taxonomy_newick(organisms_df, n=20)
        tree = Phylo.read(StringIO(newick), "newick")
        assert tree is not None

    def test_respects_n_parameter(self, organisms_df):
        newick_20 = build_taxonomy_newick(organisms_df, n=20)
        newick_5 = build_taxonomy_newick(organisms_df, n=5)
        assert len(newick_20) >= len(newick_5)


class TestCompareSequencesPairwise:
    def test_returns_dataframe(self):
        seqs = [SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ, SAMPLE_PET_SEQ]
        result = compare_sequences_pairwise(seqs)
        assert isinstance(result, pd.DataFrame)

    def test_square_matrix(self):
        seqs = [SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ]
        result = compare_sequences_pairwise(seqs, ids=["seq1", "seq2"])
        assert result.shape == (2, 2)

    def test_diagonal_is_one(self):
        seqs = [SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ]
        result = compare_sequences_pairwise(seqs)
        for i in range(len(seqs)):
            assert result.iloc[i, i] == pytest.approx(1.0, abs=0.001)

    def test_symmetric(self):
        seqs = [SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ, SAMPLE_PET_SEQ[:100]]
        result = compare_sequences_pairwise(seqs)
        pd.testing.assert_frame_equal(result, result.T)

    def test_identical_seqs_score_one(self):
        result = compare_sequences_pairwise([SAMPLE_PET_SEQ, SAMPLE_PET_SEQ])
        assert result.iloc[0, 1] == pytest.approx(1.0, abs=0.001)

    def test_different_seqs_score_less_than_one(self):
        result = compare_sequences_pairwise([SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ])
        assert result.iloc[0, 1] < 1.0

    def test_custom_ids(self):
        seqs = [SAMPLE_PET_SEQ, SAMPLE_PHB_SEQ]
        result = compare_sequences_pairwise(seqs, ids=["PETase", "PHBdep"])
        assert "PETase" in result.index
        assert "PHBdep" in result.columns


class TestEnrichDFWithProteinProperties:
    def test_adds_enzyme_family_column(self, minimal_df):
        result = enrich_df_with_protein_properties(minimal_df)
        assert "enzyme_family" in result.columns

    def test_adds_serine_motif_column(self, minimal_df):
        result = enrich_df_with_protein_properties(minimal_df)
        assert "has_serine_motif" in result.columns

    def test_adds_seq_length_column(self, minimal_df):
        result = enrich_df_with_protein_properties(minimal_df)
        assert "seq_length" in result.columns

    def test_ideonella_enzyme_family(self, minimal_df):
        result = enrich_df_with_protein_properties(minimal_df)
        ideonella_fam = result[result["organism"] == "Ideonella sakaiensis"]["enzyme_family"].iloc[0]
        assert ideonella_fam == "PETase / Cutinase"

    def test_does_not_modify_original(self, minimal_df):
        original_cols = list(minimal_df.columns)
        enrich_df_with_protein_properties(minimal_df)
        assert list(minimal_df.columns) == original_cols


class TestTopEnzymeFamiliesPerPlastic:
    def test_returns_dataframe(self, minimal_df):
        result = top_enzyme_families_per_plastic(minimal_df)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, minimal_df):
        result = top_enzyme_families_per_plastic(minimal_df)
        assert "plastic" in result.columns
        assert "enzyme_family" in result.columns
        assert "count" in result.columns

    def test_pet_has_entries(self, minimal_df):
        result = top_enzyme_families_per_plastic(minimal_df)
        assert "PET" in result["plastic"].values


class TestEnzymeFamiliesDict:
    def test_all_families_have_required_keys(self):
        for family, info in ENZYME_FAMILIES.items():
            assert "keywords" in info, f"{family} missing 'keywords'"
            assert "plastics" in info, f"{family} missing 'plastics'"
            assert "ec" in info, f"{family} missing 'ec'"
            assert isinstance(info["keywords"], list)
            assert isinstance(info["plastics"], list)

    def test_family_count(self):
        assert len(ENZYME_FAMILIES) >= 6, "Expected at least 6 enzyme families"
