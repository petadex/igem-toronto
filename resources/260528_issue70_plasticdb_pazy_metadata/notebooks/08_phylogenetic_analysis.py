# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 08 — Phylogenetic Analysis (Biopython Bio.Phylo)
#
# Constructs phylogenetic trees from PlasticDB taxonomy using Biopython's
# **Bio.Phylo** module. Analyses clade diversity, branch-length proxies
# for taxonomic novelty, and visualises the evolutionary landscape of
# plastic-degrading organisms.

# +
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).parent.parent))
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
from collections import defaultdict

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

from src.data_loader import load_all
from src.bioinformatics import build_taxonomy_newick
from src.novel_discovery import identify_phylogenetic_gaps

data = load_all()
df = data['plasticdb']
organisms = data['organisms']

print(f"Organisms: {len(organisms):,}")
print(f"Unique genera: {df['genus'].nunique()}")
# -

# ## 1. Build a Newick tree for the top 60 multi-plastic degraders

newick_str = build_taxonomy_newick(organisms, n=60)
print(f"Newick string length: {len(newick_str)} characters")
print("First 300 chars:", newick_str[:300])

tree = Phylo.read(StringIO(newick_str), "newick")
print(f"\nTree: {tree.count_terminals()} terminal nodes (leaves)")
print(f"       {len(list(tree.find_clades()))} total clades")

# ## 2. Visualise the tree

fig, ax = plt.subplots(figsize=(14, 20))
Phylo.draw(tree, axes=ax, do_show=False, branch_labels=None)
ax.set_title("Phylogenetic Tree of Top 60 Multi-Plastic Degrading Organisms\n"
             "(genus-level grouping, branch lengths = taxonomic proximity proxy)",
             fontsize=13)
plt.tight_layout()
plt.savefig('outputs/figures/08_phylo_tree.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 3. Genus-level clade statistics

print("\n--- Genus-level clade summary ---")
internal_clades = [c for c in tree.find_clades() if c.clades]
for clade in sorted(internal_clades, key=lambda c: -len(c.get_terminals()))[:15]:
    terminals = clade.get_terminals()
    print(f"  {clade.name or 'root':<30s}: {len(terminals):3d} leaves")

# ## 4. Phylogenetic gap analysis

phylo_gaps = identify_phylogenetic_gaps(organisms, df)
print(f"\nTotal genera in DB:        {len(phylo_gaps)}")
print(f"Singleton genera:          {phylo_gaps['is_singleton'].sum()} "
      f"({100*phylo_gaps['is_singleton'].mean():.1f}%)")
print(f"Unknown/rare genera:       {(~phylo_gaps['known_genus']).sum()}")

# ## 5. Phylogenetic gap priority scatter

fig, ax = plt.subplots(figsize=(10, 6))
colors = phylo_gaps['is_singleton'].map({True: '#E94F37', False: '#2E86AB'})
sc = ax.scatter(
    phylo_gaps['n_plastics'],
    phylo_gaps['discovery_priority'],
    c=colors, s=40 + phylo_gaps['n_plastics'] * 3,
    alpha=0.7
)
for _, row in phylo_gaps.head(12).iterrows():
    ax.annotate(row['genus'],
                (row['n_plastics'], row['discovery_priority']),
                xytext=(4, 4), textcoords='offset points', fontsize=7)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#E94F37', label='Singleton genus'),
                   Patch(color='#2E86AB', label='Multi-species genus')])
ax.set_xlabel('Number of Plastic Types Degraded')
ax.set_ylabel('Discovery Priority Score')
ax.set_title('Phylogenetic Gap Analysis — Discovery Priority by Plastic Breadth')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('outputs/figures/08_phylo_gaps.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 6. Species distribution within top genera

top_genera = df['genus'].value_counts().head(10).index.tolist()
print("\nSpecies count within top 10 genera:")
for genus in top_genera:
    species_list = df[df['genus'] == genus]['organism'].unique()
    print(f"  {genus:<25s}: {len(species_list):3d} species  "
          f"| e.g. {', '.join(list(species_list)[:3])}")

# ## 7. Clade depth analysis — average plasticity by phylogenetic position

genus_plastic_counts = df.groupby('genus')['plastic'].nunique().rename('n_plastics')
genus_species_counts = df.groupby('genus')['organism'].nunique().rename('n_species')
genus_df = pd.concat([genus_plastic_counts, genus_species_counts], axis=1).reset_index()
genus_df['plastics_per_species'] = genus_df['n_plastics'] / genus_df['n_species']
genus_df = genus_df.sort_values('plastics_per_species', ascending=False)

print("\nTop genera by plastic-degradation breadth per species:")
print(genus_df.head(20).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
top20_gen = genus_df.head(20)
axes[0].barh(top20_gen['genus'].iloc[::-1], top20_gen['n_plastics'].iloc[::-1],
             color=sns.color_palette('Blues', 20))
axes[0].set_xlabel('Number of Plastic Types')
axes[0].set_title('Top 20 Genera by Plastic Breadth')

axes[1].scatter(genus_df['n_species'], genus_df['n_plastics'],
                alpha=0.5, c='steelblue', s=30)
axes[1].set_xlabel('Number of Species in DB')
axes[1].set_ylabel('Number of Plastic Types')
axes[1].set_title('Species Count vs Plastic Breadth (per genus)')
for _, row in genus_df.head(8).iterrows():
    axes[1].annotate(row['genus'],
                     (row['n_species'], row['n_plastics']),
                     xytext=(4, 2), textcoords='offset points', fontsize=7)
plt.tight_layout()
plt.savefig('outputs/figures/08_genus_breadth.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 8. Phylum-level groupings (proxy from genus names)

KNOWN_PHYLA = {
    'Pseudomonas': 'Proteobacteria', 'Ralstonia': 'Proteobacteria',
    'Ideonella': 'Proteobacteria', 'Burkholderia': 'Proteobacteria',
    'Bacillus': 'Firmicutes', 'Paenibacillus': 'Firmicutes',
    'Streptomyces': 'Actinobacteria', 'Rhodococcus': 'Actinobacteria',
    'Nocardia': 'Actinobacteria', 'Thermobifida': 'Actinobacteria',
    'Aspergillus': 'Ascomycota (Fungi)', 'Penicillium': 'Ascomycota (Fungi)',
    'Fusarium': 'Ascomycota (Fungi)', 'Trichoderma': 'Ascomycota (Fungi)',
    'Phanerochaete': 'Basidiomycota (Fungi)',
    'Clostridium': 'Firmicutes', 'Lactobacillus': 'Firmicutes',
}

df['phylum_proxy'] = df['genus'].map(KNOWN_PHYLA).fillna('Other / Unclassified')
phylum_counts = df.groupby('phylum_proxy').agg(
    n_entries=('organism', 'count'),
    n_species=('organism', 'nunique'),
    n_plastics=('plastic', 'nunique'),
).reset_index().sort_values('n_entries', ascending=False)

print("\nPhylum-level breakdown:")
print(phylum_counts.to_string(index=False))

fig, ax = plt.subplots(figsize=(9, 5))
ax.barh(phylum_counts['phylum_proxy'].iloc[::-1],
        phylum_counts['n_entries'].iloc[::-1],
        color=sns.color_palette('Set2', len(phylum_counts)))
ax.set_xlabel('Number of Entries')
ax.set_title('Research Coverage by Phylum (proxy mapping)')
plt.tight_layout()
plt.savefig('outputs/figures/08_phylum_coverage.png', dpi=150, bbox_inches='tight')
plt.show()
