# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 05 — Novel Species Discovery Pipeline
#
# Applies novelty scoring, phylogenetic gap detection, and environment-based
# prioritisation to suggest organisms most worth investigating for new plastic-
# degrading capabilities.

# +
import sys; sys.path.insert(0, '..')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from src.data_loader import load_all
from src.analysis import compute_novelty_potential, evidence_quality_score
from src.novel_discovery import (
    identify_phylogenetic_gaps,
    underexplored_environments,
    plastic_specific_candidates,
    generate_discovery_report,
    format_discovery_report_text,
    PRIORITY_PLASTICS_FOR_DISCOVERY,
)

data = load_all()
df = data['plasticdb']
organisms = data['organisms']
# -

# ## Novelty potential scores
novelty = compute_novelty_potential(organisms, df)
print("Top 20 organisms by novelty potential:")
top20 = novelty.head(20)
for i, row in top20.iterrows():
    plastics = ', '.join(row['plastics_degraded']) if isinstance(row['plastics_degraded'], list) else str(row['plastics_degraded'])
    print(f"  {i+1:2}. {row['organism']:<42s}  score={row['novelty_potential']:.1f}  [{plastics}]")

# ## Novelty scatter plot
fig, ax = plt.subplots(figsize=(11, 7))
top30 = novelty.head(30)
sc = ax.scatter(top30['breadth_score'], top30['rarity_score'],
                s=top30['evidence_gap_score'] * 6 + 30,
                c=top30['novelty_potential'], cmap='viridis', alpha=0.85, zorder=3)
plt.colorbar(sc, ax=ax, label='Novelty Potential')
for _, row in top30.iterrows():
    label = ' '.join(row['organism'].split()[:2])
    ax.annotate(label, (row['breadth_score'], row['rarity_score']),
                xytext=(4, 4), textcoords='offset points', fontsize=7)
ax.set_xlabel('Plastic Breadth Score', fontsize=12)
ax.set_ylabel('Taxonomic Rarity Score', fontsize=12)
ax.set_title('Novelty Potential — Top 30 Organisms\n(bubble size = evidence gap)', fontsize=13)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('../outputs/figures/05_novelty_scatter.png', dpi=150, bbox_inches='tight')
plt.show()

# ## Phylogenetic gaps
phylo = identify_phylogenetic_gaps(organisms, df)
print("\nTop 15 phylogenetic discovery targets (singleton genera):")
print(phylo.head(15).to_string(index=False))

# ## Underexplored environments
env_gaps = underexplored_environments(df)
print("\nUnderexplored isolation environments:")
print(env_gaps.head(12)[['isolation_environment','n_species','pct_with_sequence',
                           'characterisation_gap','exploration_score']].to_string(index=False))

# ## Priority candidates for hard plastics
for plastic in PRIORITY_PLASTICS_FOR_DISCOVERY[:4]:
    candidates = plastic_specific_candidates(df, plastic)
    if not candidates.empty:
        print(f"\n--- {plastic} — Top 5 candidates ---")
        print(candidates.head(5)[['organism','has_sequence','has_enzyme','last_year',
                                   'isolation_environments','priority_score']].to_string(index=False))

# ## Full discovery report
report = generate_discovery_report(df, organisms, novelty)
print(format_discovery_report_text(report))

