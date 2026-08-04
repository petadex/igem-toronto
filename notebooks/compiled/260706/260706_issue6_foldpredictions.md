## Petadex Centroid Structural Predictions

Lead     : `AlexLeonardos`

Issue    : [Github Issue #6](https://github.com/petadex/igem-toronto/issues/) — _Generate a structural prediction for each family centroid_

Start    : `2026-07-06`

Complete : `N/A`

Files    : `~/igem-toronto/resources/260706_issue6_foldpredictions/` — folding script

S3 files : `s3://petadex-protein-structures/esmfold2-centroids/` - eventual path for these structures.

## Step 1: Access the Centroids

### Data Accessed: YYYY-MM-DD
```
# commands used to pull input data - UPDATE THIS WHEN SQL COMMANDS ARE DONE
```

## Step 2: Folding Script (Run over all Centroids)

Folds every centroid with **ESMFold2** (`esmfold2-fast-2026-05`, https://biohub.ai/models/esmfold2).

Script: `resources/260706_issue6_foldpredictions/esmfold2_predictor.py`

**Input** — the CSV from Step 1 (rows = centroid AA sequences). Columns are auto-detected
(`sequence`/`centroid_id`); pin them with `--seq-col` / `--id-col` if needed.

**Output**
- `<outdir>/structures/<id>.cif` — all-atom mmCIF straight from the SDK (`result.complex.to_mmcif()`)
- `<outdir>/arrays/<id>.npz` — `per_residue_plddt[L]`, `pae[L,L]`, `residue_index[L]`
- `<outdir>/metrics.csv` — one row per centroid: `wall_s`, `mean_plddt`, `ptm`, `has_pae`, `seq_len`, … (also the **resume ledger**)
- `structures/`, `arrays/` **and** `metrics.csv` are uploaded to **`s3://petadex-protein-structures/esmfold2-centroids/`** by default (pass `--s3 ''` to keep output local-only). `metrics.csv` is mirrored periodically + at exit and re-seeded from S3 on resume, so the ledger survives an ephemeral EC2 box.

**Inference is remote** — the biohub SDK (`esm.sdk.esmfold2_client` → `client.fold_all_atom(...)`) runs on biohub's GPUs, so the VM needs **no GPU**; it is an API client + S3 writer (needs `BIOHUB_TOKEN`). The tutorial exposes **no local-weights loader**, so `--backend local` is disabled.

**FoldingConfig knobs** (accuracy vs runtime/cost): `--num-loops` (default 20), `--num-sampling-steps` (100), `--lm-mask-pct`, `--msa-max-depth`, `--no-pae`.

**Benchmark arms** (`--arm`): `single` = single-sequence fold *(implemented)*; `msa` = **petadex MSA+ESM** *(TODO — mechanism is now known (`MSA.from_a3m(...)` / `MSA.from_sequences([...])` → `ProteinInput(msa=...)`), but the petadex MSA source is still TBD, so `--arm msa` stays disabled and `build_msa()` is a documented stub)*.

```bash
export BIOHUB_TOKEN=<biohub.ai token>

# smoke-test 5 centroids against the remote API, kept local (no upload)
python esmfold2_predictor.py centroids.csv --limit 5 --s3 ''

# full run; uploads to the default S3 sink
python esmfold2_predictor.py centroids.csv --outdir out

# faster/cheaper pass with fewer refinement loops
python esmfold2_predictor.py centroids.csv --num-loops 10
```

> ✅ SDK calls are aligned to the official biohub ESMFold2 tutorial (`esmfold2.ipynb`): `esmfold2_client` + `FoldingConfig`, `input_builder.ProteinInput`/`StructurePredictionInput`, `client.fold_all_atom(...)`, and the result's `.plddt` / `.ptm` / `.pae` / `.complex.to_mmcif()`. The earlier `# VERIFY` placeholders are resolved.

## Structure-Prediction Quality Metrics

Every folded structure — from **ESMFold2** (`esmfold2_predictor.py`) and from the **ESM Atlas** download (`esm_atlas_download.py`) — carries the same quality data, written as a `{structures/<id>.cif, arrays/<id>.npz, metrics.csv}` triple.

| Metric | Shape / range | Where it's stored | What it tells you |
|---|---|---|---|
| **pLDDT** (predicted lDDT-Cα) | per-residue, 0–100 | CIF **B-factor** + `npz:per_residue_plddt` | Per-residue *local* confidence — how reliable each residue's local geometry is. **>90** very high (side chains trustworthy); **70–90** confident backbone; **50–70** low, treat cautiously; **<50** likely disordered/unreliable. |
| **mean pLDDT** | scalar, 0–100 | `metrics.csv:mean_plddt` | Chain-average of pLDDT — a single global confidence number for ranking/filtering whole models. |
| **pTM** (predicted TM-score) | scalar, 0–1 | `metrics.csv:ptm` | Global *fold/topology* confidence — estimated TM-score vs. the true structure. **>0.5** means the overall fold is probably correct; higher = more confident topology. Computed from the PAE. |
| **PAE** (Predicted Aligned Error) | L×L, Ångströms | `npz:pae` | Pairwise confidence. Entry (i, j) = expected position error at residue *j* when aligned on residue *i*. Low off-diagonal blocks ⇒ two regions/domains are confidently placed *relative to each other*; high ⇒ each domain may be fine internally but their arrangement is uncertain. Best signal for inter-domain confidence and for carving out rigid domains. |
| **has_pae** | 0/1 flag | `metrics.csv:has_pae` | Whether a PAE matrix was present for that structure (guards downstream code against missing PAE). |

**Notes**
- **pLDDT vs pTM** — pLDDT is *local* (is this residue placed well?), pTM is *global* (is the whole fold right?). A model can have high pLDDT everywhere yet mediocre pTM if two well-formed domains are mis-oriented — which is exactly what the PAE off-diagonal reveals.
- **ipTM** (interface pTM) only applies to multi-chain complexes; the centroids are single chains, so it is not relevant here.
- **Distogram / PDE** — ESMFold2 can also emit a *distogram* (predicted inter-residue distance distribution); richer than PAE but bulky, and **not currently saved** by either script.
- Earlier versions of `esm_atlas_download.py` kept only pLDDT (in the CIF B-factor) and discarded pTM/PAE; both are now retained so Atlas and freshly-folded structures are directly comparable.


## Local Inference (planned — next implementation)

`esmfold2_predictor.py` uses the **remote** biohub API. But local GPU inference is viable and will likely be the **primary** path (a GPU VM with lots of disk), so next session we build a **separate local script**.

**Finding (HuggingFace `biohub/ESMFold2-Fast`, ~0.2B params, F32; license mit/other):** weights are downloadable and run locally via `transformers`:

```python
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

model = ESMFold2Model.from_pretrained("biohub/ESMFold2-Fast").cuda().eval()
out = model.infer_protein(sequence, num_loops=3, num_sampling_steps=50)
```

Install: `pip install esm@git+https://github.com/Biohub/esm.git@main` (+ a CUDA `torch` and `transformers`).

**Plan for the local script (next session)**
- **Reuse from `esmfold2_predictor.py`** (identical output contract → local, API, and Atlas all emit the same `{structures/<id>.cif, arrays/<id>.npz, metrics.csv}` triple): the CSV reader (`read_centroids`), the `Sink` (triple + S3 mirror/seed + resume ledger), the `metrics.csv` schema. Concurrency differs — local is GPU-serial (`workers=1`, consider batching sequences).
- **Swap the backend**: load `ESMFold2Model` once, loop centroids, call `infer_protein(seq, num_loops, num_sampling_steps)`.
- **Open questions to resolve on the VM** (the local analog of the old `# VERIFY`):
  1. What does `infer_protein(...)` **return**? A structure object with `.to_mmcif()` / `.plddt` / `.pae` like the remote result, or raw tensors? (Inspect `type(out)` and `dir(out)` on the box; that determines the local `_extract_result`.)
  2. Does local `infer_protein` expose **PAE**? The signature shows `num_loops`/`num_sampling_steps` but no `include_pae` — check whether PAE comes back and how to request it.
  3. **VRAM / instance sizing** — not stated on the model card (0.2B params is small, but MSA/complex modes cost more). Confirm on the target instance.
- Once the return shape is confirmed, the local script drops straight into the same triple + S3 layout, so nothing downstream changes.
