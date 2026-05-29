# ESM Landscape Controls

Generate negative/positive control sequence sets for the ESM-2 embedding landscape
(issue #34) to test whether the clusters seen in the issue #13 / #32 UMAP carry
biological signal or are artefacts of amino-acid composition.

Six control sets are derived from the 64,730 family-representative sequences and
written as CSVs. **Embedding and UMAP projection of these controls is tracked
separately in issue #75.**

---

## Files

| File | Purpose |
|---|---|
| `main.py` | Runs every script in `scripts/`, writing one control CSV per set to `controls/` |
| `scripts/generate_rand_matched.py` | Random AA string per centroid, **length matched 1:1** to the real sequence |
| `scripts/generate_rand_empirical.py` | Random AA strings whose lengths are **sampled from the observed length distribution** |
| `scripts/generate_rand_95th.py` | Random AA strings with lengths drawn uniformly **between the 5th–95th length percentiles** |
| `scripts/generate_shuffle.py` | Each centroid **shuffled** — preserves length and AA composition, destroys order |
| `scripts/generate_rand_fragments.py` | Random **30% / 60% / 90% fragments** of each centroid (3 per centroid) |
| `scripts/generate_rand_uniprot.py` | **n random real UniProt** proteins, fetched by accession via the UniProt REST API |

All scripts use a fixed seed (`471829`) for reproducibility.

---

## Inputs

Place under `data/` (pulled from the PETadex S3 archive):

- `data/family-representatives.csv` — columns `family_id`, `sequence` (64,730 rows)
- `data/uniprot_ids.tsv` — column `Entry`, the accession pool sampled by `generate_rand_uniprot.py`

---

## Usage

```bash
# Generate all six control sets -> controls/*.csv
python3 main.py
```

The resulting CSVs are the input to issue #75 (embed + UMAP).

---

## Control sets

| Set | Variant label | Count | Holds constant | Varies |
|---|---|---|---|---|
| Random matched | `rand_matched` | 64,730 | length (per centroid) | composition + order (uniform random AA) |
| Random empirical | `rand_empirical_family` | 64,730 | length *distribution* | composition + order |
| Random 95th | `rand_95th_family` | 64,730 | length in [p5, p95] | composition + order |
| Shuffled | `shuffled` | 64,730 | length **and composition** | order only |
| Fragments | `rand_fragment` | 194,190 | real subsequence | truncated to 30/60/90% |
| Random UniProt | `rand_uniprot` | 64,730 | (real proteins) | identity (non-PETase) |
| **Total controls** | | **517,840** | | |

---

## Prerequisites

- `pandas`, `numpy`, `requests`
- Network access (UniProt REST API) for `generate_rand_uniprot.py`
