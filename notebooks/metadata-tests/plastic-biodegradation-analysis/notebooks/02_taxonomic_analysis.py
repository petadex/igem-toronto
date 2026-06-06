# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 02 — Taxonomic Analysis
#
# Explores diversity at genus and species level, identifies under-sampled clades,
# and scores phylogenetic coverage.

# +
import sys; sys.path.insert(0, '..')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from src.data_loader import load_all
from src.analysis import taxonomic_diversity
from src.novel_discovery import identify_phylogenetic_gaps

data = load_all()
df = data['plasticdb']
organisms = data['organisms']
# -

# ## Diversity metrics
td = taxonomic_diversity(df)
for k, v in td.items():
    if not isinstance(v, dict):
        print(f"  {k}: {v}")

# ## Top genera
genus_counts = df['genus'].value_counts().head(25)
fig, ax = plt.subplots(figsize=(12, 5))
genus_counts.sort_values().plot.barh(ax=ax,
    color=sns.color_palette('teal', len(genus_counts)))
ax.set_xlabel('Number of Entries')
ax.set_title('Top 25 Genera — Plastic Biodegradation')
plt.tight_layout()
plt.savefig('../outputs/figures/02_top_genera.png', dpi=150, bbox_inches='tight')
plt.show()

# ## Species per genus distribution
spp_per_genus = df.groupby('genus')['organism'].nunique()
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(spp_per_genus, bins=range(1, 20), color='steelblue', edgecolor='white')
ax.set_xlabel('Number of Species per Genus')
ax.set_ylabel('Number of Genera')
ax.set_title('Species-per-Genus Distribution (most genera have only 1 species)')
plt.tight_layout()
plt.savefig('../outputs/figures/02_species_per_genus.png', dpi=150, bbox_inches='tight')
plt.show()
print(f"Genera with exactly 1 species: {(spp_per_genus == 1).sum()} / {len(spp_per_genus)}")

# ## Phylogenetic gaps
phylo = identify_phylogenetic_gaps(organisms, df)
print("\nTop 20 Singleton Genera (Discovery Priority):")
print(phylo.head(20).to_string(index=False))

# ## Multi-plastic generalists
multi = organisms.nlargest(20, 'n_plastics')[
    ['organism', 'genus', 'n_plastics', 'has_sequence', 'has_enzyme',
     'first_year', 'last_year']]
print("\nTop 20 multi-plastic degrading organisms:")
print(multi.to_string(index=False))

