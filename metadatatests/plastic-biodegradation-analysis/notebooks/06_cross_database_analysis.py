# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 06 — Cross-Database Analysis: PlasticDB vs PAZy
#
# Compares the broad PlasticDB coverage with the rigorously characterised
# PAZy enzyme set to identify where evidence is strong and where gaps exist.

# +
import sys; sys.path.insert(0, '..')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from src.data_loader import load_all
from src.analysis import cross_database_comparison, evidence_quality_score

data = load_all()
df = data['plasticdb']
pazy = data['pazy']
# -

# ## Cross-database comparison summary
cdb = cross_database_comparison(df, pazy)
for k, v in cdb.items():
    if not isinstance(v, list):
        print(f"  {k}: {v}")
print(f"\n  Shared plastics: {', '.join(sorted(cdb['shared_plastics']))}")
print(f"  PlasticDB-only:  {', '.join(sorted(cdb['plasticdb_only_plastics']))}")
print(f"  PAZy-only:       {', '.join(sorted(cdb['pazy_only_plastics']))}")

# ## Coverage comparison chart
all_plastics = sorted(
    set(cdb['shared_plastics']) |
    set(cdb['plasticdb_only_plastics']) |
    set(cdb['pazy_only_plastics'])
)
in_pdb  = [1 if p in cdb['shared_plastics'] + cdb['plasticdb_only_plastics'] else 0 for p in all_plastics]
in_pazy = [1 if p in cdb['shared_plastics'] + cdb['pazy_only_plastics'] else 0 for p in all_plastics]

x = np.arange(len(all_plastics))
fig, ax = plt.subplots(figsize=(14, 4))
ax.bar(x - 0.2, in_pdb,  0.4, label='PlasticDB', color='#2E86AB', alpha=0.85)
ax.bar(x + 0.2, in_pazy, 0.4, label='PAZy',      color='#A23B72', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(all_plastics, rotation=45, ha='right')
ax.set_title('Plastic Type Coverage: PlasticDB vs PAZy')
ax.legend()
ax.set_ylabel('Covered (1=Yes)')
plt.tight_layout()
plt.savefig('../outputs/figures/06_db_coverage.png', dpi=150, bbox_inches='tight')
plt.show()

# ## PAZy characterised enzyme table
print("\nPAZy — Thoroughly Characterised Plastic-Active Enzymes:")
print(pazy.to_string(index=False))

# ## Evidence quality distribution in PlasticDB
df_sc = evidence_quality_score(df)
tier_counts = df_sc['evidence_tier'].value_counts()
fig, ax = plt.subplots(figsize=(7, 5))
colors = ['#2E86AB','#44BBA4','#F18F01','#E94F37']
wedges, texts, autotexts = ax.pie(
    tier_counts.reindex(['Excellent','High','Medium','Low']).fillna(0),
    labels=['Excellent','High','Medium','Low'],
    autopct='%1.1f%%', colors=colors, startangle=90,
    wedgeprops={'width': 0.55}
)
ax.set_title('Evidence Quality Tier — PlasticDB', fontsize=13)
plt.tight_layout()
plt.savefig('../outputs/figures/06_evidence_tiers.png', dpi=150, bbox_inches='tight')
plt.show()

# ## PlasticDB entries per year with enzyme sequence data vs without
df_sc['year_int'] = df_sc['year'].fillna(0).astype(int)
yearly_seq = df_sc[df_sc['year_int'] > 1990].groupby(['year_int','has_sequence']).size().unstack(fill_value=0)
yearly_seq.columns = ['No Sequence', 'Has Sequence']
fig, ax = plt.subplots(figsize=(12, 4))
yearly_seq.plot.bar(ax=ax, stacked=True, color=['#E94F37','#2E86AB'], alpha=0.85)
ax.set_xlabel('Year'); ax.set_ylabel('Entries')
ax.set_title('Entries with vs without Protein Sequence Data per Year')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('../outputs/figures/06_sequence_coverage_over_time.png', dpi=150, bbox_inches='tight')
plt.show()

