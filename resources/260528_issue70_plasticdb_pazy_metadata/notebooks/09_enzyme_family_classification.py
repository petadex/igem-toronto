# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 09 — Enzyme Family Classification & Active-Site Analysis
#
# Classifies all enzymes in PlasticDB into mechanistic families (PETase/Cutinase,
# Lipase, PHB Depolymerase, Laccase/Oxidase, Amidase/Nylonase, etc.) using
# keyword matching on enzyme names. Scans sequences for the canonical
# **GxSxG serine hydrolase** active-site motif using Biopython regex utilities.
# Compares enzyme family distributions across plastic types and databases.

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
import re
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from src.data_loader import load_all
from src.bioinformatics import (
    classify_enzyme_family,
    enrich_df_with_protein_properties,
    top_enzyme_families_per_plastic,
    scan_serine_hydrolase_motif,
    compute_protein_properties,
    ENZYME_FAMILIES,
    SERINE_HYDROLASE_MOTIF,
)

data = load_all()
df = data['plasticdb']
pazy = data['pazy']
# -

# ## 1. Enzyme family overview

print("Defined enzyme families:")
for family, info in ENZYME_FAMILIES.items():
    print(f"  {family:<30s}  keywords={info['keywords'][:3]}  EC={info['ec']}")

# ## 2. Classify all PlasticDB entries

enriched = enrich_df_with_protein_properties(df)
family_counts = enriched['enzyme_family'].value_counts()
print("\nEnzyme family distribution:")
print(family_counts.to_string())

fig, ax = plt.subplots(figsize=(10, 5))
family_counts.sort_values().plot.barh(
    ax=ax,
    color=sns.color_palette('Set2', len(family_counts))
)
ax.set_xlabel('Number of Entries')
ax.set_title('Enzyme Family Distribution — PlasticDB')
plt.tight_layout()
plt.savefig('outputs/figures/09_enzyme_families.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 3. Enzyme family heatmap: plastic type vs enzyme family

top_plastics = df['plastic'].value_counts().head(15).index.tolist()
sub = enriched[enriched['plastic'].isin(top_plastics)]
pivot = sub.pivot_table(index='plastic', columns='enzyme_family',
                         aggfunc='size', fill_value=0)

fig, ax = plt.subplots(figsize=(14, 8))
sns.heatmap(pivot, annot=True, fmt='d', cmap='Blues', ax=ax,
            linewidths=0.3, cbar_kws={'label': 'Entry Count'})
ax.set_title('Enzyme Family Distribution by Plastic Type (Top 15 Plastics)', fontsize=13)
ax.set_xlabel('Enzyme Family')
ax.set_ylabel('Plastic Type')
plt.xticks(rotation=40, ha='right')
plt.tight_layout()
plt.savefig('outputs/figures/09_family_by_plastic_heatmap.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 4. Serine hydrolase motif scan (GxSxG)

enriched_with_seq = enriched[enriched['has_sequence'] & enriched['sequence'].notna()].copy()
n_with_motif = enriched_with_seq['has_serine_motif'].sum()
n_total_seq = len(enriched_with_seq)
print(f"\nSerine hydrolase motif (GxSxG) scan:")
print(f"  Sequences scanned:       {n_total_seq:,}")
print(f"  Sequences with motif:    {n_with_motif:,}  ({100*n_with_motif/n_total_seq:.1f}%)")
print(f"  Sequences without motif: {n_total_seq - n_with_motif:,}")

motif_by_family = enriched_with_seq.groupby('enzyme_family')['has_serine_motif'].agg(
    total='count', with_motif='sum'
)
motif_by_family['pct_with_motif'] = (motif_by_family['with_motif'] /
                                      motif_by_family['total'] * 100).round(1)
print("\nMotif frequency by enzyme family:")
print(motif_by_family.to_string())

# ## 5. Known serine hydrolase active-site residue context

def extract_serine_context(sequence: str, window: int = 5) -> list[str]:
    """Return all GxSxG motif occurrences with ±window flanking residues."""
    seq = sequence.upper()
    contexts = []
    for m in SERINE_HYDROLASE_MOTIF.finditer(seq):
        start = max(0, m.start() - window)
        end = min(len(seq), m.end() + window)
        contexts.append(f"...{seq[start:m.start()]}[{m.group()}]{seq[m.end():end]}...")
    return contexts

print("\nActive-site contexts in IsPETase:")
for ctx in extract_serine_context(
    "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDPNSYFGQNLHPYPAQDDLSGHLMGNTVEQIAQLRQEF"
    "QAAIAQRGTITIDQQPGHPHTYIQSYSDFQDAFQHYLPNVSDDQTLDDGYLFHVNAKYRDYETLMPSGKY"
    "RNVIADYQNIVKNNDLEISPDQFAGMIQDIMTADLQNFVSQYPENTLIYIIGHSMGGGLVSRTAFDQIGA"
):
    print(f"  {ctx}")

# ## 6. EC number analysis in PAZy

if 'ec_number' in pazy.columns:
    ec_counts = pazy['ec_number'].value_counts()
    print("\nPAZy EC number distribution:")
    print(ec_counts.to_string())
    fig, ax = plt.subplots(figsize=(9, 4))
    ec_counts.plot.bar(ax=ax, color=sns.color_palette('Pastel1', len(ec_counts)))
    ax.set_title('EC Number Distribution — PAZy Characterised Enzymes')
    ax.set_ylabel('Count')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig('outputs/figures/09_pazy_ec_numbers.png', dpi=150, bbox_inches='tight')
    plt.show()

# ## 7. Physicochemical comparison: PETase vs PHB-depolymerase family

def get_family_sequences(enriched_df, family):
    sub = enriched_df[
        (enriched_df['enzyme_family'] == family) &
        enriched_df['has_sequence'] &
        enriched_df['sequence'].notna()
    ]
    return sub['sequence'].tolist()

petase_seqs = get_family_sequences(enriched, 'PETase / Cutinase')
phb_seqs    = get_family_sequences(enriched, 'PHB Depolymerase')

def batch_properties(seqs):
    rows = []
    for s in seqs:
        p = compute_protein_properties(str(s))
        if p:
            rows.append(p)
    return pd.DataFrame(rows)

pet_props = batch_properties(petase_seqs[:50])
phb_props = batch_properties(phb_seqs[:50])

if not pet_props.empty and not phb_props.empty:
    compare_cols = ['molecular_weight_kda', 'isoelectric_point',
                    'instability_index', 'gravy']
    print("\n=== Physicochemical Comparison: PETase vs PHB Depolymerase ===")
    print(f"{'Property':<30s} {'PETase mean':>15s} {'PHB Depol. mean':>18s}")
    print("-" * 65)
    for col in compare_cols:
        if col in pet_props.columns and col in phb_props.columns:
            print(f"  {col:<28s} {pet_props[col].mean():>12.3f}   {phb_props[col].mean():>15.3f}")

    fig, axes = plt.subplots(1, len(compare_cols), figsize=(16, 4))
    for ax, col in zip(axes, compare_cols):
        ax.hist(pet_props[col].dropna(), bins=20, alpha=0.6,
                color='#2E86AB', label='PETase/Cutinase')
        ax.hist(phb_props[col].dropna(), bins=20, alpha=0.6,
                color='#A23B72', label='PHB Depolymerase')
        ax.set_title(col.replace('_', ' ').title())
        ax.legend(fontsize=7)
    plt.suptitle('Physicochemical Profiles: PETase vs PHB Depolymerase', y=1.02)
    plt.tight_layout()
    plt.savefig('outputs/figures/09_family_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()

# ## 8. Enzyme name frequency in DB

top_enzymes = enriched['enzyme_name'].dropna().value_counts().head(20)
print("\nTop 20 most frequent enzyme names in PlasticDB:")
print(top_enzymes.to_string())
