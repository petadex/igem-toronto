# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 04 — Geographic & Temporal Trends
#
# Tracks the growth of plastic biodegradation research over time and
# maps where discoveries are being made.

# +
import sys; sys.path.insert(0, '..')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from src.data_loader import load_all
from src.analysis import temporal_trend_analysis, geographic_distribution, isolation_environment_profile

data = load_all()
df = data['plasticdb']
# -

# ## Temporal trends
tt = temporal_trend_analysis(df)
fig, ax1 = plt.subplots(figsize=(13, 5))
ax2 = ax1.twinx()
ax1.bar(tt['year'], tt['n_entries'], color='#2E86AB', alpha=0.65, label='Entries/yr')
ax1.plot(tt['year'], tt['rolling_3yr'], color='#F18F01', lw=2.5, label='3-yr avg')
ax2.plot(tt['year'], tt['cumulative_species'], color='#A23B72', lw=2, ls='--', label='Cumulative species')
ax1.set_xlabel('Year'); ax1.set_ylabel('Entries per Year')
ax2.set_ylabel('Cumulative Unique Species')
ax1.set_title('Plastic Biodegradation Research Growth (1990–2025)', fontsize=14)
lines1, labs1 = ax1.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
ax1.legend(lines1+lines2, labs1+labs2, loc='upper left')
plt.tight_layout()
plt.savefig('../outputs/figures/04_temporal_trends.png', dpi=150, bbox_inches='tight')
plt.show()
print(tt[['year','n_entries','n_unique_species','rolling_3yr','cumulative_species']].tail(15).to_string(index=False))

# ## Growth in hard-to-degrade plastics (PE, PP, PS, PVC)
hard = ['PE','LDPE','HDPE','PP','PS','PVC']
hard_df = df[df['plastic'].isin(hard)]
hard_tt = hard_df.groupby(['year','plastic']).size().reset_index(name='n')
pivot = hard_tt.pivot(index='year', columns='plastic', values='n').fillna(0)
fig, ax = plt.subplots(figsize=(12, 5))
pivot.plot(ax=ax, marker='o', linewidth=1.5)
ax.set_title('Research Growth — Hard-to-Degrade Plastics')
ax.set_xlabel('Year'); ax.set_ylabel('Entries per Year')
plt.tight_layout()
plt.savefig('../outputs/figures/04_hard_plastic_trends.png', dpi=150, bbox_inches='tight')
plt.show()

# ## Geographic distribution
geo = geographic_distribution(df)
top30 = geo.head(30)
fig, ax = plt.subplots(figsize=(10, 11))
ax.barh(top30['isolation_location'].iloc[::-1], top30['n_entries'].iloc[::-1],
        color=sns.color_palette('teal', len(top30)))
ax.set_xlabel('Number of Entries')
ax.set_title('Research Activity by Country/Region (Top 30)')
plt.tight_layout()
plt.savefig('../outputs/figures/04_geographic.png', dpi=150, bbox_inches='tight')
plt.show()
print(geo.head(20)[['isolation_location','n_entries','n_species','n_plastics']].to_string(index=False))

# ## Isolation environment
env = isolation_environment_profile(df)
fig, ax = plt.subplots(figsize=(10, 6))
top_env = env.head(20)
ax.barh(top_env['isolation_environment'].iloc[::-1], top_env['n_entries'].iloc[::-1],
        color=sns.color_palette('muted', len(top_env)))
ax.set_xlabel('Number of Entries')
ax.set_title('Entries by Isolation Environment')
plt.tight_layout()
plt.savefig('../outputs/figures/04_environments.png', dpi=150, bbox_inches='tight')
plt.show()

