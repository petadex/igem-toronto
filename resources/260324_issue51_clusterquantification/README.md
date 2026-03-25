# Cluster Quantification (Issue #51)

This folder documents exploratory cluster quality analysis performed in `cluster_quantification.ipynb`.

## Objective

Assess how well biological groupings are separated using silhouette scores under two label schemes:

- **Domain labels**
- **Component labels**

The analysis is run in two feature spaces:

- **640D embedding space** (`embeddings.npz`)
- **2D UMAP space** (`umap_x`, `umap_y` from `family_atlas.csv`)

## Data Sources

- `data/embeddings.npz`
	- `embeddings`: numeric vector representation per family
	- `family_ids`: IDs used for joins
- `data/family_atlas.csv`
	- family metadata including `family_id`, `component`, and UMAP coordinates
- `data/component.summary.tsv`
	- mapping from `component` to `domain`

Optional database pull is included in the notebook (`SAVED=False`) to regenerate `family_atlas.csv` directly from the Petadex RDS source.

## Notebook Workflow

### 1) Load + join data

- Load embeddings from `embeddings.npz`
- Load `family_atlas.csv` (or fetch from DB)
- Join embeddings to family/component metadata by `family_id`
- Join component-to-domain annotations
- Fill missing domain values with `Unknown` where needed

### 2) Silhouette on 640D embeddings (Domain labels)

- Features: embedding dimensions
- Labels: `domain`
- Metric used in notebook: `cosine`
- Outputs:
	- overall silhouette score
	- per-cluster summary (`mean`, `median`, `count`)
	- silhouette distribution plot by cluster

### 3) Silhouette on 640D embeddings (Component labels)

- Features: embedding dimensions
- Labels: `component`
- Metric used in notebook: `cosine`
- Outputs:
	- overall silhouette score
	- per-cluster summary (`mean`, `median`, `count`)
	- silhouette distribution plot by component
	- exported table: `component_silhouette.csv`

### 4) Silhouette on 2D UMAP (Domain labels)

- Features: `umap_x`, `umap_y`
- Labels: `domain`
- Metric: default (Euclidean)
- Outputs:
	- overall silhouette score
	- per-cluster summary table
	- silhouette distribution plot

### 5) Silhouette on 2D UMAP (Component labels)

- Features: `umap_x`, `umap_y`
- Labels: `component`
- Metric: default (Euclidean)
- Outputs:
	- overall silhouette score
	- per-cluster summary table
	- silhouette distribution plot

## Generated Artifacts

- `component_silhouette.csv`: component-level silhouette summary from the embedding/component run

## Interpretation Notes

- Silhouette values range from `-1` to `1`:
	- close to `1`: well-separated clusters
	- around `0`: overlapping boundaries
	- below `0`: likely misassignment or weak separation
- In high-dimensional embedding space, cosine-distance silhouette can appear strongly negative even when structure exists; compare against reduced-space (UMAP) and/or PCA+Euclidean validation.

## Reproducibility

Run notebook top-to-bottom from this directory so relative paths to `data/` resolve correctly:

```bash
cd resources/260324_issue51_clusterquantification
jupyter notebook cluster_quantification.ipynb
```
