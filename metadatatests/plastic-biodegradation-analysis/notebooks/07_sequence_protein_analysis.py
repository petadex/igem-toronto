# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 07 — Sequence & Protein Property Analysis (Biopython)
#
# Applies **Bio.SeqUtils.ProtParam** to every protein sequence stored in PlasticDB
# to compute physicochemical properties: molecular weight, isoelectric point,
# instability index, GRAVY score, aromaticity, and secondary-structure propensity.
# Results are compared across plastic substrate types to identify systematic
# biochemical patterns in plastic-degrading enzymes.

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
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from src.data_loader import load_all
from src.bioinformatics import (
    analyse_all_sequences,
    compute_protein_properties,
    sequence_length_stats,
    sequence_property_summary,
)

data = load_all()
df = data['plasticdb']

print(f"PlasticDB: {len(df):,} entries")
print(f"Entries with sequence: {df['has_sequence'].sum():,}")
print(f"Entries with GenBank ID: {df['has_genbank'].sum():,}")
# -

# ## 1. Walk-through — IsPETase (Ideonella sakaiensis)
#
# IsPETase is the landmark PET-degrading enzyme. We compute its full property profile.

ISPETASE_SEQ = (
    "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDPNSYFGQNLHPYPAQDDLSGHLMGNTVEQIAQLRQEF"
    "QAAIAQRGTITIDQQPGHPHTYIQSYSDFQDAFQHYLPNVSDDQTLDDGYLFHVNAKYRDYETLMPSGKY"
    "RNVIADYQNIVKNNDLEISPDQFAGMIQDIMTADLQNFVSQYPENTLIYIIGHSMGGGLVSRTAFDQIGA"
    "AVDLEHPFVSKLADSIGDPIGKPSEGVSQHPQYVKTIFQNPANPLDTTAPIVTLNNTDYFLGSQMILHRY"
    "APAQGGLIFLGSRSGSAFSSEGGDRLIDVAESILSGSGDPTATVAMKNGYQIPSLAAAQMLQELYQAACRE"
)

pa = ProteinAnalysis(ISPETASE_SEQ)
props = compute_protein_properties(ISPETASE_SEQ)

print("=== IsPETase Physicochemical Profile ===")
for k, v in props.items():
    if k != 'aa_composition':
        print(f"  {k:<30s}: {v}")

print("\n  Amino acid composition:")
for aa, frac in sorted(props['aa_composition'].items(), key=lambda x: -x[1])[:10]:
    print(f"    {aa}: {frac:.3f} ({frac*100:.1f}%)")

# ## 2. Analyse all sequences in PlasticDB

print("\nRunning ProtParam on all sequences (this takes ~30s)...")
seq_df = analyse_all_sequences(df, min_length=30)
print(f"Analysed {len(seq_df):,} sequences with valid properties")
seq_df.head(5)

# ## 3. Sequence length distribution by plastic type

len_stats = sequence_length_stats(df)
print("\nSequence length statistics by plastic (top 15 by count):")
print(len_stats.head(15).to_string(index=False))

fig, ax = plt.subplots(figsize=(10, 5))
top_plastics = len_stats.head(12)['plastic'].tolist()
sub = seq_df[seq_df['plastic'].isin(top_plastics)]
sub_lens = sub.copy()
sub_lens['seq_length'] = sub_lens['length']
order = top_plastics
sns.boxplot(data=sub_lens, x='plastic', y='seq_length', order=order,
            palette='Blues', ax=ax, linewidth=0.8)
ax.set_xlabel('Plastic Type')
ax.set_ylabel('Protein Sequence Length (aa)')
ax.set_title('Protein Sequence Length Distribution by Plastic Substrate')
ax.tick_params(axis='x', rotation=45)
plt.tight_layout()
plt.savefig('outputs/figures/07_seq_length_by_plastic.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 4. Molecular weight distribution

if not seq_df.empty:
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    props_to_plot = [
        ('molecular_weight_kda', 'Molecular Weight (kDa)'),
        ('isoelectric_point', 'Isoelectric Point (pI)'),
        ('instability_index', 'Instability Index'),
        ('gravy', 'GRAVY Score'),
        ('aromaticity', 'Aromaticity'),
        ('helix_fraction', 'Helix Fraction'),
    ]
    for ax, (col, title) in zip(axes.flat, props_to_plot):
        if col in seq_df.columns:
            ax.hist(seq_df[col].dropna(), bins=40, color='steelblue',
                    edgecolor='white', alpha=0.85)
            ax.axvline(seq_df[col].mean(), color='red', ls='--', lw=1.5,
                       label=f'mean={seq_df[col].mean():.2f}')
            ax.set_xlabel(title)
            ax.set_ylabel('Count')
            ax.set_title(title)
            ax.legend(fontsize=8)
    plt.suptitle('Physicochemical Property Distributions — PlasticDB Sequences', y=1.02)
    plt.tight_layout()
    plt.savefig('outputs/figures/07_protparam_distributions.png', dpi=150, bbox_inches='tight')
    plt.show()

    print("\nOverall statistics:")
    for col, title in props_to_plot:
        if col in seq_df.columns:
            print(f"  {title:<35s}: mean={seq_df[col].mean():.3f}  "
                  f"std={seq_df[col].std():.3f}  "
                  f"min={seq_df[col].min():.3f}  max={seq_df[col].max():.3f}")

# ## 5. Stability analysis — instability index threshold < 40 = stable

if not seq_df.empty:
    stable = seq_df[seq_df['instability_index'] < 40]
    unstable = seq_df[seq_df['instability_index'] >= 40]
    print(f"\nStable enzymes (II < 40):   {len(stable):4d} ({100*len(stable)/len(seq_df):.1f}%)")
    print(f"Unstable enzymes (II >= 40): {len(unstable):4d} ({100*len(unstable)/len(seq_df):.1f}%)")

    stable_by_plastic = seq_df.groupby('plastic')['is_stable'].mean().sort_values(ascending=False)
    print("\nFraction stable by plastic (top 15):")
    print(stable_by_plastic.head(15).round(3).to_string())

# ## 6. GRAVY score: hydrophobic vs hydrophilic bias

if not seq_df.empty:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter: pI vs GRAVY coloured by stability
    colors = seq_df['is_stable'].map({True: '#2E86AB', False: '#E94F37'})
    axes[0].scatter(seq_df['gravy'], seq_df['isoelectric_point'],
                    c=colors, alpha=0.5, s=15)
    axes[0].axvline(0, color='gray', ls='--', lw=0.8)
    axes[0].set_xlabel('GRAVY Score (+ = hydrophobic)')
    axes[0].set_ylabel('Isoelectric Point (pI)')
    axes[0].set_title('GRAVY vs pI — Colour = Stability')
    from matplotlib.patches import Patch
    axes[0].legend(handles=[Patch(color='#2E86AB', label='Stable (II<40)'),
                              Patch(color='#E94F37', label='Unstable (II≥40)')])

    # MW vs instability index
    axes[1].scatter(seq_df['molecular_weight_kda'], seq_df['instability_index'],
                    alpha=0.4, s=15, c='steelblue')
    axes[1].axhline(40, color='red', ls='--', lw=1.5, label='Stability threshold (40)')
    axes[1].set_xlabel('Molecular Weight (kDa)')
    axes[1].set_ylabel('Instability Index')
    axes[1].set_title('MW vs Instability Index')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig('outputs/figures/07_stability_scatter.png', dpi=150, bbox_inches='tight')
    plt.show()

# ## 7. Property summary by plastic type

if not seq_df.empty:
    summary = sequence_property_summary(seq_df)
    print("\nProperty summary by plastic (top 10 by sequence count):")
    cols_to_show = [c for c in summary.columns if 'count' in c or 'mean' in c]
    print(summary.nlargest(10, 'molecular_weight_kda_count')[['plastic'] + cols_to_show].to_string(index=False))

# ## 8. Secondary structure propensity heatmap

if not seq_df.empty:
    ss_cols = ['helix_fraction', 'turn_fraction', 'sheet_fraction']
    top_plas = seq_df['plastic'].value_counts().head(12).index.tolist()
    ss_mean = seq_df[seq_df['plastic'].isin(top_plas)].groupby('plastic')[ss_cols].mean()
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(ss_mean.round(3), annot=True, fmt='.3f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Mean Fraction'})
    ax.set_title('Mean Secondary Structure Propensity by Plastic Type')
    ax.set_xlabel('Structure Element')
    plt.tight_layout()
    plt.savefig('outputs/figures/07_secondary_structure.png', dpi=150, bbox_inches='tight')
    plt.show()
