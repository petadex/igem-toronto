"""
Tests for src/data_loader.py

Covers: loading, cleaning, column presence, data types, aggregation utilities.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_loader import (
    load_plasticdb,
    fetch_pazy_proteins,
    get_unique_organisms,
    get_plastic_summary,
    load_all,
    PLASTIC_FULL_NAMES,
    PLASTIC_CATEGORIES,
)


class TestLoadPlasticDB:
    def test_returns_dataframe(self, plasticdb_df):
        assert isinstance(plasticdb_df, pd.DataFrame)

    def test_minimum_row_count(self, plasticdb_df):
        assert len(plasticdb_df) > 2000, "Expected >2000 entries from PlasticDB"

    def test_required_columns_present(self, plasticdb_df):
        required = [
            "organism", "tax_id", "plastic", "year",
            "has_sequence", "has_genbank", "has_enzyme",
            "genus", "species", "plastic_category", "plastic_full_name", "decade",
        ]
        for col in required:
            assert col in plasticdb_df.columns, f"Missing column: {col}"

    def test_year_is_numeric(self, plasticdb_df):
        assert pd.api.types.is_float_dtype(plasticdb_df["year"]) or \
               pd.api.types.is_integer_dtype(plasticdb_df["year"])

    def test_year_range_sensible(self, plasticdb_df):
        years = plasticdb_df["year"].dropna()
        assert years.min() >= 1960
        assert years.max() <= 2030

    def test_boolean_flag_columns(self, plasticdb_df):
        for col in ["has_sequence", "has_genbank", "has_enzyme"]:
            assert plasticdb_df[col].dtype in [bool, np.bool_, object], \
                f"{col} should be boolean"

    def test_plastic_column_not_all_null(self, plasticdb_df):
        assert plasticdb_df["plastic"].notna().sum() > 1000

    def test_organism_column_not_all_null(self, plasticdb_df):
        assert plasticdb_df["organism"].notna().sum() > 1000

    def test_genus_extracted(self, plasticdb_df):
        assert plasticdb_df["genus"].notna().mean() > 0.8

    def test_known_plastic_types_present(self, plasticdb_df):
        found = set(plasticdb_df["plastic"].dropna().unique())
        for expected in ["PET", "PHA", "PHB", "LDPE", "PU"]:
            assert expected in found, f"Expected plastic '{expected}' missing"

    def test_plastic_category_assigned(self, plasticdb_df):
        cats = plasticdb_df["plastic_category"].unique().tolist()
        assert "Biodegradable/Bio-based" in cats
        assert "Commodity Thermoplastics" in cats

    def test_decade_column(self, plasticdb_df):
        decades = plasticdb_df["decade"].dropna().unique()
        assert 1990 in decades or 2000 in decades or 2010 in decades

    def test_no_duplicate_required_columns(self, plasticdb_df):
        assert len(plasticdb_df.columns) == len(set(plasticdb_df.columns))

    def test_minimal_df_columns(self, minimal_df):
        assert "organism" in minimal_df.columns
        assert "plastic" in minimal_df.columns
        assert len(minimal_df) == 6


class TestFetchPazyProteins:
    def test_returns_dataframe(self, pazy_df):
        assert isinstance(pazy_df, pd.DataFrame)

    def test_has_entries(self, pazy_df):
        assert len(pazy_df) > 0

    def test_has_expected_columns(self, pazy_df):
        for col in ["enzyme_name", "organism", "plastic"]:
            assert col in pazy_df.columns, f"PAZy missing column: {col}"

    def test_known_enzymes_present(self, pazy_df):
        names = pazy_df["enzyme_name"].str.lower().tolist()
        assert any("petase" in n or "lcc" in n for n in names), \
            "Expected IsPETase or LCC in PAZy"

    def test_pet_plastic_covered(self, pazy_df):
        assert "PET" in pazy_df["plastic"].values


class TestGetUniqueOrganisms:
    def test_returns_dataframe(self, organisms_df):
        assert isinstance(organisms_df, pd.DataFrame)

    def test_one_row_per_organism(self, organisms_df):
        assert organisms_df["organism"].is_unique

    def test_n_plastics_column(self, organisms_df):
        assert "n_plastics" in organisms_df.columns
        assert organisms_df["n_plastics"].min() >= 1

    def test_plastics_degraded_is_list(self, organisms_df):
        sample = organisms_df["plastics_degraded"].dropna().iloc[0]
        assert isinstance(sample, list)

    def test_has_sequence_flag_aggregated(self, organisms_df):
        assert "has_sequence" in organisms_df.columns

    def test_minimal_organisms_unique(self, minimal_organisms):
        assert minimal_organisms["organism"].is_unique
        assert len(minimal_organisms) == 5

    def test_multi_plastic_organism_counted(self, minimal_organisms):
        ideonella = minimal_organisms[
            minimal_organisms["organism"] == "Ideonella sakaiensis"
        ]
        assert len(ideonella) == 1
        assert ideonella.iloc[0]["n_plastics"] == 2


class TestGetPlasticSummary:
    def test_returns_dataframe(self, plastics_df):
        assert isinstance(plastics_df, pd.DataFrame)

    def test_required_columns(self, plastics_df):
        for col in ["plastic", "n_entries", "n_unique_organisms", "pct_with_sequence"]:
            assert col in plastics_df.columns

    def test_sorted_by_entries_descending(self, plastics_df):
        entries = plastics_df["n_entries"].tolist()
        assert entries == sorted(entries, reverse=True)

    def test_pct_with_sequence_range(self, plastics_df):
        assert plastics_df["pct_with_sequence"].between(0, 100).all()


class TestLoadAll:
    def test_returns_dict_with_all_keys(self, data):
        for key in ["plasticdb", "pazy", "organisms", "plastics"]:
            assert key in data

    def test_all_values_are_dataframes(self, data):
        for key, val in data.items():
            assert isinstance(val, pd.DataFrame), f"data['{key}'] is not a DataFrame"

    def test_plastic_full_names_mapping(self):
        assert "PET" in PLASTIC_FULL_NAMES
        assert "PHB" in PLASTIC_FULL_NAMES
        assert PLASTIC_FULL_NAMES["PET"] == "Polyethylene terephthalate"

    def test_plastic_categories_mapping(self):
        assert "Biodegradable/Bio-based" in PLASTIC_CATEGORIES
        assert "PHB" in PLASTIC_CATEGORIES["Biodegradable/Bio-based"]
        assert "PET" in PLASTIC_CATEGORIES["Commodity Thermoplastics"]
