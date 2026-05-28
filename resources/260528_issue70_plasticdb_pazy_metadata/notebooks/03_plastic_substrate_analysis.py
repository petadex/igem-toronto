# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 03 — Plastic Substrate Analysis
#
# Investigates which plastic types are best covered, which share degraders,
# and which present major research gaps.

# +
import sys; sys.path.insert(0, '..')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from src.data_loader import load_all, PLASTIC_CATEGORIES
from src.analysis import plastic_co_occurrence, research_gap_analysis, evidence_quality_score
from src.visualization import plot_co_occurrence_heatmap, plot_research_gaps

data = load_all()
df = data['plasticdb']
organisms = data['organisms']
# -

# ## Plastic type summary
ps = data['plastics']
print("Plastic summary (top 20):")
print(ps.head(20).to_string(index=False))

# ## Category breakdown
cat_counts = df.groupby('plastic_category').size().sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(8, 4))
cat_counts.plot.bar(ax=ax, color=sns.color_palette('Set2', len(cat_counts)))
ax.set_title('Entries by Plastic Category')
ax.set_ylabel('Entries')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.savefig('../outputs/figures/03_plastic_categories.png', dpi=150, bbox_inches='tight')
plt.show()

# ## Co-occurrence matrix
co = plastic_co_occurrence(organisms)
top15 = co.sum().nlargest(15).index
sub = co.loc[top15, top15].copy()
np.fill_diagonal(sub.values, 0)
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(sub, annot=True, fmt='d', cmap='Blues', ax=ax,
            linewidths=0.4, cbar_kws={'label': 'Shared Organisms'})
ax.set_title('Plastic Co-occurrence (organisms degrading both types)', fontsize=13)
plt.tight_layout()
plt.savefig('../outputs/figures/03_cooccurrence.png', dpi=150, bbox_inches='tight')
plt.show()

# ## Research gaps
gaps = research_gap_analysis(df)
print("\nTop 10 plastics by research gap score:")
print(gaps['plastic_gaps'].head(10)[['plastic','n_organisms','pct_with_sequence',
                                      'pct_extrapolated','last_year','gap_score']].to_string(index=False))

# ## Evidence quality by plastic
df_sc = evidence_quality_score(df)
ev_by_plastic = df_sc.groupby('plastic')['evidence_score'].mean().sort_values().head(20)
fig, ax = plt.subplots(figsize=(10, 5))
ev_by_plastic.sort_values().plot.barh(ax=ax, color='coral')
ax.set_xlabel('Mean Evidence Score')
ax.set_title('Mean Evidence Quality Score by Plastic Type (lowest = most under-characterised)')
plt.tight_layout()
plt.savefig('../outputs/figures/03_evidence_by_plastic.png', dpi=150, bbox_inches='tight')
plt.show()

