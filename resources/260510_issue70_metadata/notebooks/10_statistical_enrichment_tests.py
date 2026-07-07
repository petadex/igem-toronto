# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 10 — Statistical Enrichment & Hypothesis Testing
#
# Applies rigorous statistical tests to detect non-random patterns in the
# plastic biodegradation literature:
# - **Chi-square** tests for geographic and taxonomic enrichment
# - **Fisher's exact** test for 2×2 plastic × evidence contingency tables
# - **Mann-Whitney U** test for comparing evidence scores across plastic categories
# - **Spearman correlation** between temporal activity and geographic diversity
# - **Kruskal-Wallis** ANOVA for physicochemical property differences across families
# - **Bonferroni / FDR** multiple-testing correction

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
from scipy import stats
from scipy.stats import (
    chi2_contingency, fisher_exact, mannwhitneyu,
    spearmanr, kruskal, shapiro
)
from statsmodels.stats.multitest import multipletests

from src.data_loader import load_all
from src.analysis import evidence_quality_score
from src.bioinformatics import analyse_all_sequences

data = load_all()
df = data['plasticdb']
df_sc = evidence_quality_score(df)
# -

# ## 1. Normality check — evidence scores

scores = df_sc['evidence_score'].dropna()
stat, p_shapiro = shapiro(scores.sample(min(5000, len(scores)), random_state=42))
print(f"Shapiro-Wilk normality test on evidence scores:")
print(f"  W = {stat:.4f},  p = {p_shapiro:.4e}")
print(f"  Distribution is {'NOT ' if p_shapiro < 0.05 else ''}normal at α=0.05")
print(f"  → Using non-parametric tests throughout\n")

# ## 2. Chi-square: are certain plastic types over-represented in top-productivity countries?

top_countries = ['India', 'Japan', 'South Korea', 'France', 'Thailand', 'China']
top_plastics_chi = df_sc['plastic'].value_counts().head(6).index.tolist()
sub_chi = df_sc[
    df_sc['isolation_location'].isin(top_countries) &
    df_sc['plastic'].isin(top_plastics_chi)
]
contingency = pd.crosstab(sub_chi['isolation_location'], sub_chi['plastic'])
chi2, p_chi2, dof, expected = chi2_contingency(contingency)
print(f"Chi-square test: country × plastic type")
print(f"  χ² = {chi2:.3f},  df = {dof},  p = {p_chi2:.4e}")
print(f"  {'Significant' if p_chi2 < 0.05 else 'Not significant'} at α=0.05")
print(f"\nObserved counts:")
print(contingency.to_string())

fig, ax = plt.subplots(figsize=(11, 6))
sns.heatmap(contingency, annot=True, fmt='d', cmap='Blues', ax=ax,
            linewidths=0.3)
ax.set_title(f'Country × Plastic Type Contingency\n(χ²={chi2:.1f}, p={p_chi2:.2e})')
plt.tight_layout()
plt.savefig('outputs/figures/10_country_plastic_contingency.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 3. Fisher's exact: is PHB better characterised than LDPE?

def fisher_evidence(df_in, plastic_a, plastic_b):
    """2×2 Fisher exact: has_sequence vs plastic_a vs plastic_b."""
    a_seq = ((df_in['plastic'] == plastic_a) & df_in['has_sequence']).sum()
    a_no  = ((df_in['plastic'] == plastic_a) & ~df_in['has_sequence'].fillna(False)).sum()
    b_seq = ((df_in['plastic'] == plastic_b) & df_in['has_sequence']).sum()
    b_no  = ((df_in['plastic'] == plastic_b) & ~df_in['has_sequence'].fillna(False)).sum()
    table = [[a_seq, a_no], [b_seq, b_no]]
    odds_ratio, p_val = fisher_exact(table, alternative='two-sided')
    return {'plastic_a': plastic_a, 'plastic_b': plastic_b,
            'odds_ratio': round(odds_ratio, 3), 'p_value': round(p_val, 6),
            f'{plastic_a}_seq%': round(100*a_seq/(a_seq+a_no), 1) if (a_seq+a_no) > 0 else 0,
            f'{plastic_b}_seq%': round(100*b_seq/(b_seq+b_no), 1) if (b_seq+b_no) > 0 else 0}

pairs = [('PHB', 'LDPE'), ('PET', 'PS'), ('PCL', 'PVC'), ('PHB', 'PP'), ('PLA', 'PE')]
fisher_results = [fisher_evidence(df_sc, a, b) for a, b in pairs]
fisher_df = pd.DataFrame(fisher_results)
print("\nFisher's Exact Tests — Sequence Evidence Rate Comparison:")
print(fisher_df.to_string(index=False))

# ## 4. Mann-Whitney U — evidence scores: biodegradable vs commodity plastics

biodeg = df_sc[df_sc['plastic_category'] == 'Biodegradable/Bio-based']['evidence_score'].dropna()
commodity = df_sc[df_sc['plastic_category'] == 'Commodity Thermoplastics']['evidence_score'].dropna()
u_stat, p_mw = mannwhitneyu(biodeg, commodity, alternative='two-sided')
r_effect = u_stat / (len(biodeg) * len(commodity))
print(f"\nMann-Whitney U — Biodegradable vs Commodity plastic evidence scores:")
print(f"  Biodegradable  mean={biodeg.mean():.2f}  n={len(biodeg)}")
print(f"  Commodity      mean={commodity.mean():.2f}  n={len(commodity)}")
print(f"  U = {u_stat:.0f},  p = {p_mw:.4e}")
print(f"  Effect size r = {r_effect:.4f}  ({'large' if r_effect > 0.3 else 'medium' if r_effect > 0.1 else 'small'})")

fig, ax = plt.subplots(figsize=(8, 5))
ax.boxplot([biodeg, commodity], labels=['Biodegradable/Bio-based', 'Commodity Thermoplastics'],
           patch_artist=True,
           boxprops=dict(facecolor='lightblue'),
           medianprops=dict(color='red', lw=2))
ax.set_ylabel('Evidence Score')
ax.set_title(f'Evidence Quality: Biodegradable vs Commodity Plastics\n'
             f'Mann-Whitney U p={p_mw:.2e}')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('outputs/figures/10_evidence_by_category.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 5. Kruskal-Wallis ANOVA across top isolation environments

top_envs = df_sc['isolation_environment'].value_counts().head(6).index.tolist()
groups = [
    df_sc[df_sc['isolation_environment'] == env]['evidence_score'].dropna().values
    for env in top_envs
]
groups = [g for g in groups if len(g) > 5]
if len(groups) >= 3:
    h_stat, p_kw = kruskal(*groups)
    print(f"\nKruskal-Wallis ANOVA — evidence score across isolation environments:")
    print(f"  H = {h_stat:.3f},  p = {p_kw:.4e}")
    print(f"  {'Significant' if p_kw < 0.05 else 'Not significant'} differences")

# ## 6. Multiple pairwise comparisons (all plastic pairs), FDR correction

plastics_for_pairs = df_sc['plastic'].value_counts().head(10).index.tolist()
pvalues = []
test_labels = []
for i in range(len(plastics_for_pairs)):
    for j in range(i+1, len(plastics_for_pairs)):
        pa, pb = plastics_for_pairs[i], plastics_for_pairs[j]
        ga = df_sc[df_sc['plastic'] == pa]['evidence_score'].dropna().values
        gb = df_sc[df_sc['plastic'] == pb]['evidence_score'].dropna().values
        if len(ga) > 5 and len(gb) > 5:
            _, pv = mannwhitneyu(ga, gb, alternative='two-sided')
            pvalues.append(pv)
            test_labels.append(f"{pa} vs {pb}")

if pvalues:
    reject, pvals_corrected, _, _ = multipletests(pvalues, method='fdr_bh')
    pairwise_df = pd.DataFrame({
        'comparison': test_labels,
        'p_raw': [round(p, 5) for p in pvalues],
        'p_fdr': [round(p, 5) for p in pvals_corrected],
        'significant': reject,
    }).sort_values('p_fdr')
    print(f"\nPairwise evidence score comparisons ({len(pvalues)} tests), FDR corrected:")
    print(pairwise_df.to_string(index=False))

    sig_count = reject.sum()
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.scatter(range(len(pvalues)), sorted(pvals_corrected), s=20, c=reject,
               cmap='RdYlGn_r', alpha=0.8)
    ax.axhline(0.05, color='red', ls='--', lw=1.5, label='FDR threshold 0.05')
    ax.set_xlabel('Test rank')
    ax.set_ylabel('FDR-corrected p-value')
    ax.set_title(f'Multiple Testing Correction (Benjamini-Hochberg)\n'
                 f'{sig_count}/{len(pvalues)} significant comparisons')
    ax.legend()
    plt.tight_layout()
    plt.savefig('outputs/figures/10_pairwise_fdr.png', dpi=150, bbox_inches='tight')
    plt.show()

# ## 7. Spearman correlation: year vs geographic diversity

yearly_geo = df_sc.groupby('year').agg(
    n_countries=('isolation_location', 'nunique'),
    n_species=('organism', 'nunique'),
    n_plastics=('plastic', 'nunique'),
).reset_index().dropna()

rho, p_sp = spearmanr(yearly_geo['year'], yearly_geo['n_countries'])
print(f"\nSpearman ρ (year vs geographic diversity): {rho:.4f},  p={p_sp:.4e}")
rho2, p_sp2 = spearmanr(yearly_geo['year'], yearly_geo['n_species'])
print(f"Spearman ρ (year vs unique species/yr):    {rho2:.4f},  p={p_sp2:.4e}")

fig, ax = plt.subplots(figsize=(9, 4))
ax.scatter(yearly_geo['year'], yearly_geo['n_countries'],
           s=yearly_geo['n_species'] / 2, alpha=0.6, c='steelblue')
ax.set_xlabel('Year')
ax.set_ylabel('Number of Countries Reporting')
ax.set_title(f'Geographic Reach vs Year  (Spearman ρ={rho:.3f}, p={p_sp:.2e})')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('outputs/figures/10_temporal_geo_correlation.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 8. Effect size summary table

print("\n=== Statistical Results Summary ===")
print(f"{'Test':<50s}  {'Statistic':>12s}  {'p-value':>12s}  {'Significant':>12s}")
print("-" * 95)
print(f"  Chi-square (country × plastic){'':20s}  χ²={chi2:>7.1f}    {p_chi2:>10.2e}   {'Yes' if p_chi2<0.05 else 'No':>12s}")
print(f"  Mann-Whitney (biodeg vs commodity evid.){'':9s}  U={u_stat:>8.0f}    {p_mw:>10.2e}   {'Yes' if p_mw<0.05 else 'No':>12s}")
print(f"  Spearman (year vs geo diversity){'':17s}  ρ={rho:>8.4f}    {p_sp:>10.2e}   {'Yes' if p_sp<0.05 else 'No':>12s}")
if len(groups) >= 3:
    print(f"  Kruskal-Wallis (evid. by environment){'':12s}  H={h_stat:>8.3f}    {p_kw:>10.2e}   {'Yes' if p_kw<0.05 else 'No':>12s}")
