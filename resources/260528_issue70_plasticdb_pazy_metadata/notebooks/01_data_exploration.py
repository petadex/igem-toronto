# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 01 — Data Exploration
# **Sources:** PlasticDB (plasticdb.org) · PAZy (pazy.eu)
#
# This notebook loads and explores both databases, showing raw distributions,
# completeness, and basic structure.

# +
import sys; sys.path.insert(0, '..')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from src.data_loader import load_all, PLASTIC_FULL_NAMES

data = load_all()
df = data['plasticdb']
pazy = data['pazy']
print(f"PlasticDB: {len(df):,} entries | {df['organism'].nunique()} species | {df['plastic'].nunique()} plastics")
print(f"PAZy:      {len(pazy):,} entries")
# -

# ## Column overview
df.dtypes

# ## Null / completeness analysis
nulls = df.isnull().mean().sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(10, 5))
nulls.plot.bar(ax=ax, color='steelblue')
ax.set_ylabel('Fraction missing')
ax.set_title('Data completeness — PlasticDB')
ax.axhline(0.5, color='red', ls='--', label='50% missing')
ax.legend()
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('../outputs/figures/01_completeness.png', dpi=150, bbox_inches='tight')
plt.show()
print(nulls.round(3).to_string())

# ## Plastic type distribution
plastic_counts = df['plastic'].value_counts().head(25)
fig, ax = plt.subplots(figsize=(12, 5))
plastic_counts.plot.bar(ax=ax, color=sns.color_palette('Blues_r', len(plastic_counts)))
ax.set_xlabel('Plastic Type')
ax.set_ylabel('Number of Entries')
ax.set_title('Plastic Type Distribution — PlasticDB (top 25)')
plt.tight_layout()
plt.savefig('../outputs/figures/01_plastic_dist.png', dpi=150, bbox_inches='tight')
plt.show()

# ## Evidence quality flags
for col in ['has_sequence', 'has_genbank', 'has_enzyme', 'analytical_grade', 'thermophilic']:
    val = df[col].sum() if df[col].dtype == bool else df[col].fillna(False).sum()
    print(f"  {col:<30s}: {val:4d} ({100*val/len(df):.1f}%)")

# ## PAZy overview
print("\nPAZy enzyme entries:")
print(pazy.to_string())

