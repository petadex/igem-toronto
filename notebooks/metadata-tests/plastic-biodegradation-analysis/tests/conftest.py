"""
Shared pytest fixtures for the plastic biodegradation analysis test suite.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_loader import load_all, get_unique_organisms, get_plastic_summary


@pytest.fixture(scope="session")
def data():
    """Load all datasets once per test session (cached)."""
    return load_all()


@pytest.fixture(scope="session")
def plasticdb_df(data):
    return data["plasticdb"]


@pytest.fixture(scope="session")
def pazy_df(data):
    return data["pazy"]


@pytest.fixture(scope="session")
def organisms_df(data):
    return data["organisms"]


@pytest.fixture(scope="session")
def plastics_df(data):
    return data["plastics"]


@pytest.fixture
def minimal_df():
    """A tiny synthetic DataFrame with known values for unit testing."""
    return pd.DataFrame({
        "organism": ["Ideonella sakaiensis", "Ideonella sakaiensis",
                     "Pseudomonas putida", "Bacillus subtilis",
                     "Aspergillus niger", "Thermobifida fusca"],
        "tax_id": [1547922, 1547922, 303, 1423, 5061, 2061],
        "plastic": ["PET", "PHB", "LDPE", "PET", "LDPE", "PET"],
        "year": [2016, 2016, 2018, 2020, 2021, 2012],
        "genus": ["Ideonella", "Ideonella", "Pseudomonas",
                  "Bacillus", "Aspergillus", "Thermobifida"],
        "species": ["sakaiensis", "sakaiensis", "putida",
                    "subtilis", "niger", "fusca"],
        "has_sequence": [True, True, False, True, False, True],
        "has_genbank": [True, True, False, True, False, True],
        "has_enzyme": [True, True, False, True, False, True],
        "analytical_grade": [True, True, False, True, False, True],
        "extrapolated_from_enzyme": [False, False, True, False, True, False],
        "isolation_environment": ["Plastic waste", "Soil", "Soil", "Compost", "Soil", "Compost"],
        "isolation_location": ["Japan", "USA", "India", "Germany", "India", "Germany"],
        "enzyme_name": ["IsPETase", "PHB depolymerase", "", "PETase", "", "TfCut2"],
        "sequence": [
            "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDPNSYFGQNLHPYPAQDDLSGHLMGNTVEQIAQLRQEF"
            "QAAIAQRGTITIDQQPGHPHTYIQSYSDFQDAFQHYLPNVSDDQTLDDGYLFHVNAKYRDYETLMPSGKY"
            "RNVIADYQNIVKNNDLEISPDQFAGMIQDIMTADLQNFVSQYPENTLIYIIGHSMGGGLVSRTAFDQIGA",
            "MKHPYGYRWHWLYALVVTLMTALATFSAHAAVTAGPGAWSSQQTWAADTVNGGNLTGYFYWPASQPTTPNG"
            "KRALVLVLHGCLQTASGDVIDNANGAGFNWKTIAEQYGAVVLAPNATGNVYSNHCWDYANTSPSRTSG",
            "",
            "MKHPYGYRWHWLYALVVTLMTALATFSAHAAVTAGPGAWSSQQTWAADTVNG",
            "",
            "MSSGRPGAAAGLLALLAAVLAFSGAVAHADAADTPPAATAAPAASAAPAASAAPAASAAPASAAPSAS"
            "AASGSTGSTVGSTVSSTVGSTVSSTVGSTVSSTVGSTVSSTVGSAAPAAPAAQVTLNLGYPAASGKV",
        ],
        "plastic_category": ["Commodity Thermoplastics", "Biodegradable/Bio-based",
                              "Commodity Thermoplastics", "Commodity Thermoplastics",
                              "Commodity Thermoplastics", "Commodity Thermoplastics"],
        "plastic_full_name": ["Polyethylene terephthalate", "Polyhydroxybutyrate",
                               "Low-density polyethylene", "Polyethylene terephthalate",
                               "Low-density polyethylene", "Polyethylene terephthalate"],
        "decade": [2010, 2010, 2010, 2020, 2020, 2010],
    })


@pytest.fixture
def minimal_organisms(minimal_df):
    return get_unique_organisms(minimal_df)


SAMPLE_PET_SEQ = (
    "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDPNSYFGQNLHPYPAQDDLSGHLMGNTVEQIAQLRQEF"
    "QAAIAQRGTITIDQQPGHPHTYIQSYSDFQDAFQHYLPNVSDDQTLDDGYLFHVNAKYRDYETLMPSGKY"
    "RNVIADYQNIVKNNDLEISPDQFAGMIQDIMTADLQNFVSQYPENTLIYIIGHSMGGGLVSRTAFDQIGA"
)

SAMPLE_PHB_SEQ = (
    "MKHPYGYRWHWLYALVVTLMTALATFSAHAAVTAGPGAWSSQQTWAADTVNGGNLTGYFYWPASQPTTPNG"
    "KRALVLVLHGCLQTASGDVIDNANGAGFNWKTIAEQYGAVVLAPNATGNVYSNHCWDYANTSPSRTSG"
)
