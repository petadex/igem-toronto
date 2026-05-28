"""
Advanced bioinformatics analysis using Biopython.

Covers:
- Protein physicochemical properties (ProtParam)
- Amino acid composition analysis
- Sequence-based enzyme family classification
- FASTA parsing and sequence statistics
- Phylogenetic tree utilities from taxonomy strings
- GenBank record utilities
- Sequence motif scanning
- Comparative sequence statistics across plastic types
"""

from __future__ import annotations
import re
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from io import StringIO

warnings.filterwarnings("ignore")

DATA_DIR = Path(__file__).parent.parent / "data"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

ENZYME_FAMILIES = {
    "PETase / Cutinase": {
        "keywords": ["petase", "cutinase", "cut", "phl"],
        "plastics": ["PET", "PES", "PBSA", "PBAT", "PEF"],
        "ec": ["3.1.1.74", "3.1.1.101"],
        "motif": "GxSxG",
    },
    "Lipase": {
        "keywords": ["lipase", "lip", "pue", "pfl", "pla"],
        "plastics": ["PU", "PCL", "PLA", "PBSA", "PES"],
        "ec": ["3.1.1.3"],
        "motif": "GxSxG",
    },
    "PHB Depolymerase": {
        "keywords": ["phb", "pha", "depolymerase", "phaz", "poly"],
        "plastics": ["PHB", "PHA", "PHO", "PHBV", "PCL"],
        "ec": ["3.1.1.75", "3.1.1.76"],
        "motif": "GxSxG",
    },
    "Laccase / Oxidase": {
        "keywords": ["laccase", "oxidase", "lac", "mnp", "lmco"],
        "plastics": ["PE", "PS", "PVC", "LDPE", "HDPE"],
        "ec": ["1.10.3.2"],
        "motif": None,
    },
    "Protease": {
        "keywords": ["protease", "serine protease", "subtilisin", "alkaline protease"],
        "plastics": ["Nylon", "PU"],
        "ec": ["3.4.-.-"],
        "motif": None,
    },
    "Esterase": {
        "keywords": ["esterase", "est", "arylesterase", "carboxylesterase"],
        "plastics": ["PET", "PCL", "PBS", "PBSA", "PLA"],
        "ec": ["3.1.1.1", "3.1.1.6"],
        "motif": "GxSxG",
    },
    "Amidase / Nylonase": {
        "keywords": ["nylonase", "nyl", "amidase", "aminopeptidase"],
        "plastics": ["Nylon"],
        "ec": ["3.5.1.84", "3.5.1.-"],
        "motif": None,
    },
    "Peroxidase / MnP": {
        "keywords": ["peroxidase", "mnp", "manganese", "ligninase", "ldp"],
        "plastics": ["PE", "PS", "PVC", "Nylon", "LDPE"],
        "ec": ["1.11.1.13", "1.11.1.16"],
        "motif": None,
    },
}

SERINE_HYDROLASE_MOTIF = re.compile(r"G[ACDEFGHIKLMNPQRSTVWY]S[ACDEFGHIKLMNPQRSTVWY]G")


def clean_sequence(seq: str) -> str:
    """Remove non-amino-acid characters and return uppercase."""
    if not isinstance(seq, str):
        return ""
    return "".join(c for c in seq.upper() if c in STANDARD_AA)


def compute_protein_properties(sequence: str) -> dict:
    """
    Compute physicochemical properties via Bio.SeqUtils.ProtParam.

    Returns dict with:
      - length, molecular_weight_kda, isoelectric_point
      - instability_index, gravy (grand average of hydropathicity)
      - aromaticity, helix_fraction, turn_fraction, sheet_fraction
      - charge_at_ph7, n_positive, n_negative
      - aa_composition (dict of fractions)
    """
    seq = clean_sequence(sequence)
    if len(seq) < 10:
        return {}
    try:
        pa = ProteinAnalysis(seq)
        mw = pa.molecular_weight()
        pi = pa.isoelectric_point()
        ii = pa.instability_index()
        gravy = pa.gravy()
        arom = pa.aromaticity()
        ss = pa.secondary_structure_fraction()
        aa_counts = pa.count_amino_acids()
        total_aa = sum(aa_counts.values()) or 1
        aa_comp = {k: v / total_aa for k, v in aa_counts.items()}

        pos_aa = sum(seq.count(aa) for aa in "KRH")
        neg_aa = sum(seq.count(aa) for aa in "DE")

        return {
            "length": len(seq),
            "molecular_weight_kda": round(mw / 1000, 2),
            "isoelectric_point": round(pi, 3),
            "instability_index": round(ii, 2),
            "is_stable": ii < 40,
            "gravy": round(gravy, 4),
            "is_hydrophobic": gravy > 0,
            "aromaticity": round(arom, 4),
            "helix_fraction": round(ss[0], 4),
            "turn_fraction": round(ss[1], 4),
            "sheet_fraction": round(ss[2], 4),
            "n_positive_residues": pos_aa,
            "n_negative_residues": neg_aa,
            "charge_balance": pos_aa - neg_aa,
            "aa_composition": {k: round(v, 4) for k, v in aa_comp.items()},
        }
    except Exception:
        return {}


def analyse_all_sequences(df: pd.DataFrame, min_length: int = 20) -> pd.DataFrame:
    """
    Run ProtParam on every sequence in the PlasticDB DataFrame.
    Returns a per-entry DataFrame with physicochemical properties.
    """
    rows = []
    has_seq = df[df["has_sequence"] & df["sequence"].notna()].copy()
    for _, row in has_seq.iterrows():
        props = compute_protein_properties(str(row["sequence"]))
        if props:
            record = {
                "organism": row["organism"],
                "plastic": row["plastic"],
                "enzyme_name": row.get("enzyme_name", ""),
                "genbank_id": row.get("genbank_id", ""),
                "year": row.get("year"),
                "sequence": str(row["sequence"]),
            }
            record.update(props)
            rows.append(record)
    if not rows:
        return pd.DataFrame()
    result = pd.DataFrame(rows)
    result = result[result["length"] >= min_length].reset_index(drop=True)
    return result


def classify_enzyme_family(enzyme_name: str, plastic: str) -> str:
    """
    Assign an enzyme family based on enzyme name keywords and plastic substrate.
    Returns the best matching family name or 'Unknown'.
    """
    if not isinstance(enzyme_name, str):
        enzyme_name = ""
    name_lower = enzyme_name.lower()
    plastic_str = str(plastic).upper() if plastic else ""

    scores = {}
    for family, info in ENZYME_FAMILIES.items():
        score = 0
        for kw in info["keywords"]:
            if kw in name_lower:
                score += 3
        if any(p in plastic_str for p in info["plastics"]):
            score += 1
        scores[family] = score

    best = max(scores, key=scores.get)
    keyword_score = sum(
        3 for kw in ENZYME_FAMILIES[best]["keywords"] if kw in name_lower
    )
    return best if keyword_score > 0 else "Unknown"


def scan_serine_hydrolase_motif(sequence: str) -> bool:
    """Return True if sequence contains the GxSxG serine hydrolase active-site motif."""
    seq = clean_sequence(sequence)
    return bool(SERINE_HYDROLASE_MOTIF.search(seq))


def get_amino_acid_composition_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute mean amino acid composition per plastic type across all sequences.
    """
    rows = []
    for _, row in df[df["has_sequence"] & df["sequence"].notna()].iterrows():
        seq = clean_sequence(str(row["sequence"]))
        if len(seq) < 20:
            continue
        total = len(seq)
        comp = {aa: seq.count(aa) / total for aa in sorted(STANDARD_AA)}
        comp["plastic"] = row["plastic"]
        comp["organism"] = row["organism"]
        rows.append(comp)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def sequence_length_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Return per-plastic length statistics for all sequences."""
    has_seq = df[df["has_sequence"] & df["sequence"].notna()].copy()
    has_seq["seq_length"] = has_seq["sequence"].apply(
        lambda s: len(clean_sequence(str(s)))
    )
    stats = (
        has_seq[has_seq["seq_length"] > 10]
        .groupby("plastic")["seq_length"]
        .agg(
            count="count",
            mean=lambda x: round(x.mean(), 1),
            median="median",
            std=lambda x: round(x.std(), 1),
            min="min",
            max="max",
        )
        .reset_index()
        .sort_values("count", ascending=False)
    )
    return stats


def build_taxonomy_newick(organisms_df: pd.DataFrame, n: int = 50) -> str:
    """
    Build a simple Newick-format tree from organism names using the
    genus-level grouping as proxy clades (no NCBI lookup — offline-safe).

    Groups species under their genus, genera under their first-letter group.
    Returns a Newick string suitable for Bio.Phylo.
    """
    top = organisms_df.nlargest(n, "n_plastics").copy()
    top["genus"] = top["organism"].str.extract(r"^(\w+)")

    genus_groups: dict[str, list[str]] = defaultdict(list)
    for _, row in top.iterrows():
        g = row["genus"] if pd.notna(row["genus"]) else "Unknown"
        species = "_".join(str(row["organism"]).split()[:2])
        genus_groups[g].append(species)

    clades = []
    for genus, species_list in genus_groups.items():
        if len(species_list) == 1:
            leaf = Clade(name=species_list[0], branch_length=1.0)
            clades.append(leaf)
        else:
            inner = Clade(
                name=genus,
                branch_length=1.0,
                clades=[Clade(name=s, branch_length=0.5) for s in species_list],
            )
            clades.append(inner)

    root = Clade(name="root", branch_length=0.0, clades=clades)
    tree = Tree(root=root, rooted=True)

    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue().strip()


def parse_plasticdb_fasta(fasta_path: str | Path) -> pd.DataFrame:
    """
    Parse the PlasticDB FASTA file (PlasticDB.fasta from plasticdb.org/downloaddata).
    Returns a DataFrame with id, description, and sequence per record.
    """
    records = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        records.append({
            "id": rec.id,
            "description": rec.description,
            "sequence": str(rec.seq),
            "length": len(rec.seq),
        })
    return pd.DataFrame(records)


def compare_sequences_pairwise(seqs: list[str], ids: list[str] | None = None) -> pd.DataFrame:
    """
    Compute pairwise sequence identity between a list of sequences using
    a simple k-mer overlap method (k=5) — no alignment needed.
    Returns a square similarity DataFrame.
    """
    n = len(seqs)
    clean = [clean_sequence(s) for s in seqs]
    labels = ids if ids else [f"seq_{i}" for i in range(n)]

    def kmer_similarity(a: str, b: str, k: int = 5) -> float:
        if len(a) < k or len(b) < k:
            return 0.0
        ka = set(a[i: i + k] for i in range(len(a) - k + 1))
        kb = set(b[i: i + k] for i in range(len(b) - k + 1))
        return len(ka & kb) / len(ka | kb) if ka | kb else 0.0

    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            sim = kmer_similarity(clean[i], clean[j])
            matrix[i, j] = sim
            matrix[j, i] = sim
    return pd.DataFrame(matrix, index=labels, columns=labels)


def enrich_df_with_protein_properties(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add enzyme_family, has_serine_motif, and (where possible) basic ProtParam
    columns to the main PlasticDB DataFrame.
    """
    df = df.copy()
    df["enzyme_family"] = df.apply(
        lambda r: classify_enzyme_family(r.get("enzyme_name", ""), r.get("plastic", "")),
        axis=1,
    )
    df["has_serine_motif"] = df["sequence"].apply(
        lambda s: scan_serine_hydrolase_motif(str(s)) if pd.notna(s) else False
    )
    df["seq_length"] = df["sequence"].apply(
        lambda s: len(clean_sequence(str(s))) if pd.notna(s) else 0
    )
    return df


def top_enzyme_families_per_plastic(df: pd.DataFrame) -> pd.DataFrame:
    """Return enzyme family distribution for each plastic type."""
    enriched = enrich_df_with_protein_properties(df)
    result = (
        enriched.groupby(["plastic", "enzyme_family"])
        .size()
        .reset_index(name="count")
        .sort_values(["plastic", "count"], ascending=[True, False])
    )
    return result


def sequence_property_summary(seq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate physicochemical properties by plastic type.
    seq_df is output of analyse_all_sequences().
    """
    if seq_df.empty:
        return pd.DataFrame()
    numeric_cols = [
        "molecular_weight_kda", "isoelectric_point", "instability_index",
        "gravy", "aromaticity", "helix_fraction", "sheet_fraction",
        "length",
    ]
    available = [c for c in numeric_cols if c in seq_df.columns]
    result = (
        seq_df.groupby("plastic")[available]
        .agg(["mean", "std", "count"])
        .round(3)
    )
    result.columns = ["_".join(c) for c in result.columns]
    return result.reset_index()
