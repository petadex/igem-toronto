"""
PlasticDB + PAZy Metadata Analysis — Issue #70
Run this script to reproduce all analysis outputs.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "scripts"))

import pandas as pd
from data_loader import load_plasticdb, load_pazy
from analysis import (
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

DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "outputs"
OUT_DIR.mkdir(exist_ok=True)


def main():
    df = load_plasticdb(DATA_DIR / "plasticdb_microorganisms.tsv")
    pazy = load_pazy(DATA_DIR / "pazy_proteins.csv")

    df_scored = evidence_quality_score(df)
    df_scored.to_csv(OUT_DIR / "plasticdb_scored.csv", index=False)
    print(f"PlasticDB entries: {len(df)}")

    td = taxonomic_diversity(df)
    print(f"Unique genera: {td['n_unique_genera']}  species: {td['n_unique_species']}")
    print(f"Shannon diversity (genus): {td['genus_shannon_diversity']}")

    gaps = research_gap_analysis(df)
    gaps["plastic_gaps"].to_csv(OUT_DIR / "research_gaps.csv", index=False)
    print(f"Research gaps written to {OUT_DIR / 'research_gaps.csv'}")

    cdb = cross_database_comparison(df, pazy)
    print(f"Shared plastics (PlasticDB ∩ PAZy): {sorted(cdb['shared_plastics'])}")
    print(f"PlasticDB-only plastics: {len(cdb['plasticdb_only_plastics'])}")

    print("Done.")


if __name__ == "__main__":
    main()
