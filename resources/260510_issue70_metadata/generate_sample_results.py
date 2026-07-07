"""
Generate sample analysis results, figures, and report files.

Run: python generate_sample_results.py

Outputs:
  outputs/figures/     — PNG charts (matplotlib, publication-ready)
  outputs/reports/     — CSV tables
  outputs/discovery_report.txt
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import json

from src.data_loader import load_all
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
    run_full_analysis,
)
from src.visualization import save_all_figures, FIGURE_DIR
from src.novel_discovery import generate_discovery_report, format_discovery_report_text

REPORT_DIR = Path(__file__).parent / "outputs" / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("=" * 60)
    print("  Plastic Biodegradation Meta-Analysis")
    print("  Generating sample results...")
    print("=" * 60)

    print("\n[1/6] Loading data from PlasticDB and PAZy...")
    data = load_all()
    df = data["plasticdb"]
    organisms = data["organisms"]
    pazy = data["pazy"]
    print(f"      PlasticDB: {len(df):,} entries, {df['organism'].nunique()} species, "
          f"{df['plastic'].nunique()} plastics")
    print(f"      PAZy:      {len(pazy):,} characterised enzyme entries")

    print("\n[2/6] Running full analysis...")
    results = run_full_analysis(data)

    print("\n[3/6] Saving CSV reports...")
    df_scored = results["df_scored"]
    df_scored.drop(columns=["sequence"], errors="ignore").to_csv(
        REPORT_DIR / "plasticdb_scored.csv", index=False)
    data["plastics"].to_csv(REPORT_DIR / "plastic_summary.csv", index=False)
    results["temporal_trends"].to_csv(REPORT_DIR / "temporal_trends.csv", index=False)
    results["geographic_distribution"].to_csv(REPORT_DIR / "geographic_distribution.csv", index=False)
    results["isolation_environments"].to_csv(REPORT_DIR / "isolation_environments.csv", index=False)
    results["research_gaps"]["plastic_gaps"].to_csv(REPORT_DIR / "research_gaps.csv", index=False)
    results["novelty_scores"].to_csv(REPORT_DIR / "novelty_scores.csv", index=False)
    pazy.to_csv(REPORT_DIR / "pazy_proteins.csv", index=False)

    td = results["taxonomic_diversity"]
    with open(REPORT_DIR / "taxonomic_diversity.json", "w") as f:
        json.dump({k: v for k, v in td.items() if not isinstance(v, dict)} |
                  {"top_10_genera": td["top_10_genera"],
                   "top_10_species": td["top_10_species"]}, f, indent=2)

    cross = results["cross_db"]
    with open(REPORT_DIR / "cross_database_comparison.json", "w") as f:
        json.dump({k: v for k, v in cross.items() if isinstance(v, (int, float, str, list))},
                  f, indent=2)

    print(f"      Saved to {REPORT_DIR}/")

    print("\n[4/6] Saving sample figures (matplotlib / publication quality)...")
    paths = save_all_figures(results, backend="matplotlib")
    for p in paths:
        if p:
            print(f"      {p}")

    print("\n[5/6] Generating novel discovery report...")
    discovery = generate_discovery_report(df, organisms, results["novelty_scores"])
    report_text = format_discovery_report_text(discovery)
    report_path = Path(__file__).parent / "outputs" / "discovery_report.txt"
    report_path.write_text(report_text)
    print(report_text)

    print("\n[6/6] Printing key findings summary...")
    _print_summary(results, td, cross)

    print(f"\nAll outputs saved to: {Path(__file__).parent / 'outputs'}/")
    print("Run the Streamlit app: streamlit run app.py --server.port 5000")


def _print_summary(results, td, cross):
    print("\n" + "=" * 60)
    print("  KEY FINDINGS SUMMARY")
    print("=" * 60)
    print(f"\nTaxonomy:")
    print(f"  • {td['n_unique_genera']} unique genera, {td['n_unique_species']} species")
    print(f"  • Shannon diversity (genus level): {td['genus_shannon_diversity']:.3f}")
    print(f"  • Singleton genera: {td['singleton_genera']} ({td['pct_singleton_genera']}%) — "
          f"phylogenetically undersampled")

    tt = results["temporal_trends"]
    peak_year = tt.loc[tt["n_entries"].idxmax(), "year"]
    print(f"\nTemporal:")
    print(f"  • Peak publication year: {peak_year:.0f} "
          f"({int(tt.loc[tt['n_entries'].idxmax(), 'n_entries'])} entries)")
    recent = tt[tt["year"] >= 2020]["n_entries"].sum()
    print(f"  • Entries since 2020: {int(recent)} "
          f"({100*recent/len(results['df_scored']):.0f}% of all entries)")

    top_gaps = results["research_gaps"]["plastic_gaps"].head(5)
    print(f"\nTop 5 Research Gaps (by gap score):")
    for _, row in top_gaps.iterrows():
        print(f"  • {row['plastic']:<8s}  organisms={int(row['n_organisms']):<4d}  "
              f"seq_coverage={row['pct_with_sequence']:.0%}  "
              f"gap_score={row['gap_score']:.1f}")

    print(f"\nCross-database:")
    print(f"  • PlasticDB: {cross['plasticdb_n_organisms']} species vs "
          f"PAZy: {cross['pazy_n_organisms']} characterised species")
    print(f"  • Coverage ratio: {cross['coverage_ratio']:.1%} of plastics "
          f"appear in both databases")
    print(f"  • PlasticDB-only plastics (need characterisation): "
          f"{', '.join(sorted(cross['plasticdb_only_plastics'])[:10])}")

    top5_novel = results["novelty_scores"].head(5)
    print(f"\nTop 5 Novel Discovery Candidates:")
    for i, row in top5_novel.iterrows():
        plastics = ", ".join(row["plastics_degraded"]) if isinstance(row["plastics_degraded"], list) else str(row["plastics_degraded"])
        print(f"  {i+1}. {row['organism']:<40s}  "
              f"score={row['novelty_potential']:.1f}  "
              f"plastics=[{plastics}]")
    print("=" * 60)


if __name__ == "__main__":
    main()
