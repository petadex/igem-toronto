# ---
# jupyter:
#   jupytext:
#     formats: py:light
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Notebook 12 — Protein Biochemistry Deep Dive (Biopython ProtParam + Sequence Analysis)
#
# Deep biochemical characterisation of all plastic-degrading enzymes in PlasticDB:
# - Amino acid composition comparison across plastic substrate types
# - k-mer sequence similarity between enzyme families
# - Hydrophobicity profiles and predicted signal peptide analysis
# - Correlation between protein properties and evidence quality
# - Sequence identity clustering using Biopython tools
# - Comparative analysis of well-studied (PAZy) vs broader (PlasticDB) enzyme sets

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
from Bio.Seq import Seq
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy import stats

from src.data_loader import load_all
from src.bioinformatics import (
    analyse_all_sequences,
    get_amino_acid_composition_table,
    compare_sequences_pairwise,
    compute_protein_properties,
    clean_sequence,
    STANDARD_AA,
)
from src.analysis import evidence_quality_score

data = load_all()
df = data['plasticdb']
pazy = data['pazy']
df_sc = evidence_quality_score(df)

print("Running protein analysis (may take ~45s)...")
seq_df = analyse_all_sequences(df, min_length=30)
print(f"Analysed {len(seq_df):,} enzyme sequences from PlasticDB")
# -

# ## 1. Amino acid composition heatmap across plastic types

aa_comp_table = get_amino_acid_composition_table(df)
if not aa_comp_table.empty:
    top_plas = df['plastic'].value_counts().head(12).index.tolist()
    aa_cols = sorted(STANDARD_AA)
    aa_mean = (
        aa_comp_table[aa_comp_table['plastic'].isin(top_plas)]
        .groupby('plastic')[aa_cols]
        .mean()
    )
    fig, ax = plt.subplots(figsize=(18, 7))
    sns.heatmap(
        aa_mean.round(3),
        cmap='RdYlBu_r', ax=ax, annot=True, fmt='.3f',
        linewidths=0.2, annot_kws={'size': 7},
        cbar_kws={'label': 'Mean Fraction'},
    )
    ax.set_title('Amino Acid Composition by Plastic Substrate Type\n'
                 '(fraction of total residues, top 12 plastics)', fontsize=13)
    ax.set_xlabel('Amino Acid (1-letter code)')
    ax.set_ylabel('Plastic Type')
    plt.tight_layout()
    plt.savefig('outputs/figures/12_aa_composition_heatmap.png', dpi=150, bbox_inches='tight')
    plt.show()

    # Which AAs differ most across plastic types?
    aa_variance = aa_mean.var().sort_values(ascending=False)
    print("\nAmino acids with highest variance across plastic types:")
    print(aa_variance.head(10).round(5).to_string())

# ## 2. Sliding window hydrophobicity profile
#
# Kyte-Doolittle hydrophobicity scale — positive = hydrophobic

KD_SCALE = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

def kd_window(sequence: str, window: int = 9) -> np.ndarray:
    """Kyte-Doolittle hydrophobicity plot with sliding window."""
    seq = clean_sequence(sequence)
    scores = [KD_SCALE.get(aa, 0) for aa in seq]
    averaged = []
    half = window // 2
    for i in range(len(seq)):
        start, end = max(0, i - half), min(len(seq), i + half + 1)
        averaged.append(np.mean(scores[start:end]))
    return np.array(averaged)

PET_SEQ = (
    "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDPNSYFGQNLHPYPAQDDLSGHLMGNTVEQIAQLRQEF"
    "QAAIAQRGTITIDQQPGHPHTYIQSYSDFQDAFQHYLPNVSDDQTLDDGYLFHVNAKYRDYETLMPSGKY"
    "RNVIADYQNIVKNNDLEISPDQFAGMIQDIMTADLQNFVSQYPENTLIYIIGHSMGGGLVSRTAFDQIGA"
)
LCC_SEQ = (
    "VDAFRNANGAAAGSATSNPSPYKVNLVFNGSIHCSSAGQYGGQNYTLNVTPQNRGIFDNYYLDGSKQIRY"
    "AWSKNADSIKDTLRVNSTADLAQSGTYLNAQTLPWQNLNAWTASQNVSSPTQNQMIPVSGFQTLSDNHDT"
    "AMRGGGSTSDGTASGQGSTLNIQSAGKAVFEIPENVKGYKPTTDVFIGYHSGQGNAGLVNSSFYVDSASNV"
    "GLMGGASSGPAYATSITPPQTLFNEFLNLYQQSGGLGGATAGIPAQFLPNTPITLSSGASPTNTRGSVKTY"
    "ASTLNNTTSYEIQGLMQSGQNVSGSGPLTLNLPNGGVQNGTFTTGKGGSATVSQQVQNNNISTISSGTGF"
    "YRSSATPGGASSATATLSTSATNNVTSATSGASGSTVSATSGASGFATSGASGAATSGASGFATSGASGSAT"
)

fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=False)
for ax, seq, label in zip(axes, [PET_SEQ, LCC_SEQ], ['IsPETase', 'LCC (leaf-branch compost)']):
    profile = kd_window(seq)
    x = np.arange(len(profile))
    ax.fill_between(x, profile, 0, where=profile > 0, color='#E94F37', alpha=0.6, label='Hydrophobic')
    ax.fill_between(x, profile, 0, where=profile < 0, color='#2E86AB', alpha=0.6, label='Hydrophilic')
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xlim(0, len(profile))
    ax.set_ylabel('KD Score')
    ax.set_title(f'Kyte-Doolittle Hydrophobicity Profile — {label}')
    ax.legend(fontsize=8)
plt.suptitle('Sliding-Window Hydrophobicity (window=9)', fontsize=13)
plt.tight_layout()
plt.savefig('outputs/figures/12_hydrophobicity_profiles.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 3. k-mer similarity matrix for representative PET-degrading enzymes

KNOWN_PETASES = {
    'IsPETase':   PET_SEQ,
    'LCC':        LCC_SEQ[:200],
    'TfCut2':     ("MSGTATSRSGLPPPAPARRLAAALALAASAVLLAAPAFAADTLNPIGSSLTYSKLLAKVPTPSYHNFSPNTLG"
                   "NIFNIGNDAYVEHDAATKKGVMFDFKDTFTSLLQQAGGFDGAFKDMIQNYHLDPSQDAGRIEIEYLGASGG"),
    'PHB_dep':    ("MKHPYGYRWHWLYALVVTLMTALATFSAHAAVTAGPGAWSSQQTWAADTVNGGNLTGYFYWPASQPTTPNG"
                   "KRALVLVLHGCLQTASGDVIDNANGAGFNWKTIAEQYGAVVLAPNATGNVYSNHCWDYANTSPSRTSGHV"),
    'PCL_dep':    ("MKYSLLALVITFAAASAQAADTASAVNAATSPAQTAATVSQAPTPDGTPQTTQEGIDFHGLTYSPAQAKLS"
                   "AIAATDQLASAAKQQMTADIQQAYAAAQDAKAKAAADAKAKAAADLKAKAAADAKAKAAADAK"),
}
labels = list(KNOWN_PETASES.keys())
seqs   = list(KNOWN_PETASES.values())

sim_matrix = compare_sequences_pairwise(seqs, ids=labels)
print("\nk-mer similarity matrix (k=5) for representative enzymes:")
print(sim_matrix.round(4).to_string())

fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(
    sim_matrix.round(3), annot=True, fmt='.3f', cmap='Blues',
    ax=ax, vmin=0, vmax=1, linewidths=0.5,
    cbar_kws={'label': 'k-mer Jaccard Similarity (k=5)'},
)
ax.set_title('Sequence Similarity (k-mer Jaccard, k=5)\nRepresentative Plastic-Active Enzymes')
plt.tight_layout()
plt.savefig('outputs/figures/12_kmer_similarity.png', dpi=150, bbox_inches='tight')
plt.show()

# ## 4. Hierarchical clustering of enzyme sequences (from PlasticDB)

if not seq_df.empty:
    petase_rows = seq_df[seq_df['enzyme_name'].str.lower().str.contains(
        'petase|cutinase|lcc|phl|tfcut', na=False
    )].head(30)
    if len(petase_rows) >= 5:
        seq_list = petase_rows['sequence'].tolist()
        id_list  = [f"{row['enzyme_name'][:12]}|{row['organism'][:15]}"
                    for _, row in petase_rows.iterrows()]
        sim = compare_sequences_pairwise(seq_list, ids=id_list)
        dist_matrix = 1 - sim.values
        np.fill_diagonal(dist_matrix, 0)
        condensed = squareform(dist_matrix, checks=False)
        Z = linkage(condensed, method='ward')
        fig, ax = plt.subplots(figsize=(12, 7))
        dendrogram(Z, labels=id_list, ax=ax, leaf_rotation=45, leaf_font_size=8)
        ax.set_title('Hierarchical Clustering of PETase/Cutinase Sequences (k-mer distance)')
        ax.set_ylabel('Distance')
        plt.tight_layout()
        plt.savefig('outputs/figures/12_petase_dendrogram.png', dpi=150, bbox_inches='tight')
        plt.show()

# ## 5. Correlation between protein properties and evidence score

if not seq_df.empty:
    merged = seq_df.merge(
        df_sc[['organism', 'plastic', 'evidence_score']].drop_duplicates(
            subset=['organism', 'plastic']),
        on=['organism', 'plastic'], how='left'
    )
    prop_cols = ['molecular_weight_kda', 'isoelectric_point',
                 'instability_index', 'gravy', 'aromaticity', 'length']
    print("\nSpearman correlations with evidence score:")
    corr_rows = []
    for col in prop_cols:
        if col in merged.columns:
            mask = merged[[col, 'evidence_score']].dropna()
            if len(mask) > 10:
                rho, pval = stats.spearmanr(mask[col], mask['evidence_score'])
                corr_rows.append({'property': col, 'rho': round(rho, 4),
                                  'p_value': round(pval, 5),
                                  'significant': pval < 0.05})
    corr_df = pd.DataFrame(corr_rows).sort_values('rho', key=abs, ascending=False)
    print(corr_df.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['#2E86AB' if sig else '#AAAAAA' for sig in corr_df['significant']]
    ax.barh(corr_df['property'], corr_df['rho'], color=colors)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlabel('Spearman ρ')
    ax.set_title('Protein Property Correlation with Evidence Score\n(blue = significant at p<0.05)')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/figures/12_property_evidence_correlation.png', dpi=150, bbox_inches='tight')
    plt.show()

# ## 6. PAZy vs PlasticDB protein property comparison

def pazy_properties(pazy_df):
    rows = []
    if 'sequence' not in pazy_df.columns:
        return pd.DataFrame()
    for _, row in pazy_df[pazy_df['sequence'].notna()].iterrows():
        props = compute_protein_properties(str(row['sequence']))
        if props:
            props['source'] = 'PAZy'
            props['plastic'] = row.get('plastic', '')
            rows.append(props)
    return pd.DataFrame(rows)

pazy_props = pazy_properties(pazy)
if not seq_df.empty:
    pdb_sample = seq_df.sample(min(100, len(seq_df)), random_state=42).copy()
    pdb_sample['source'] = 'PlasticDB'
    compare_df = pd.concat([pdb_sample, pazy_props], ignore_index=True) if not pazy_props.empty \
        else pdb_sample

    if len(compare_df['source'].unique()) > 1:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for ax, col in zip(axes, ['molecular_weight_kda', 'isoelectric_point', 'instability_index']):
            for src, grp in compare_df.groupby('source'):
                ax.hist(grp[col].dropna(), bins=25, alpha=0.6, label=src)
            ax.set_title(col.replace('_', ' ').title())
            ax.legend()
        plt.suptitle('PAZy vs PlasticDB Protein Properties', fontsize=13)
        plt.tight_layout()
        plt.savefig('outputs/figures/12_pazy_vs_plasticdb_props.png', dpi=150, bbox_inches='tight')
        plt.show()

# ## 7. Signal peptide prediction (N-terminal hydrophobicity proxy)

def predict_signal_peptide(sequence: str, h_region_window: int = 8,
                            threshold: float = 1.5) -> dict:
    """
    Heuristic signal peptide detection based on N-terminal hydrophobicity.
    A true signal peptide has: n-region (charged), h-region (hydrophobic), c-region.
    This is a proxy — for definitive prediction use SignalP.
    """
    seq = clean_sequence(sequence)
    if len(seq) < 25:
        return {'has_signal_peptide': False, 'confidence': 'N/A'}
    n_term = seq[:30]
    scores = [KD_SCALE.get(aa, 0) for aa in n_term]
    max_hydro = max(
        np.mean(scores[i:i+h_region_window])
        for i in range(len(scores) - h_region_window + 1)
    )
    return {
        'has_signal_peptide': max_hydro > threshold,
        'n_terminal_max_hydro': round(max_hydro, 3),
        'confidence': 'high' if max_hydro > 2.5 else 'low',
    }

if not seq_df.empty:
    sp_results = seq_df.head(200).apply(
        lambda row: pd.Series(predict_signal_peptide(str(row.get('sequence', '')))),
        axis=1
    )
    sp_rate = sp_results['has_signal_peptide'].mean()
    print(f"\nSignal peptide prediction (N-terminal hydrophobicity heuristic):")
    print(f"  {100*sp_rate:.1f}% of sequences predicted to have signal peptide")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(sp_results['n_terminal_max_hydro'].dropna(), bins=30,
            color='steelblue', edgecolor='white', alpha=0.8)
    ax.axvline(1.5, color='red', ls='--', lw=1.5, label='Signal peptide threshold')
    ax.set_xlabel('Max sliding-window hydrophobicity (N-terminal 30 aa)')
    ax.set_ylabel('Count')
    ax.set_title('N-terminal Hydrophobicity Distribution (signal peptide proxy)')
    ax.legend()
    plt.tight_layout()
    plt.savefig('outputs/figures/12_signal_peptide_hydro.png', dpi=150, bbox_inches='tight')
    plt.show()
