# ESM Landscape

A protein embedding visualization and search tool leveraging ESM (Evolutionary Scale Modeling) embeddings, UMAP dimensionality reduction, and Qdrant vector search.

## Overview

This project processes protein embeddings to create an interactive landscape of protein families. It:

1. **Loads pre-computed ESM embeddings** from NumPy files
2. **Reduces dimensionality** using UMAP for 2D visualization
3. **Enriches data** with representative sequences and family metadata
4. **Uploads to Qdrant** for semantic similarity search
5. **Generates visualizations** of the protein embedding space

## Project Structure

```
esm-landscape/
├── data/
│   ├── embeddings.npz          # Pre-computed protein embeddings and family IDs
│   ├── representatives.csv     # Representative sequence per protein family
│   └── families.csv            # Enriched data with family sizes and coordinates
├── notebooks/
│   └── umap_vis.ipynb          # Main analysis notebook: UMAP visualization + enrichment
├── scripts/
│   └── upload_to_qdrant.py     # Upload embeddings to Qdrant vector database
└── plots/
    └── first_plot.png          # UMAP scatter plot visualization
```

## Data Files

### `embeddings.npz`
NumPy archive containing:
- `embeddings`: Dense matrix of ESM embeddings (N × D dimensions)
- `family_ids`: Array of protein family identifiers

### `representatives.csv`
CSV with columns:
- `family_id`: Unique protein family identifier
- `sequence`: Representative amino acid sequence for that family

### `families.csv` (Generated)
Enriched dataset with columns:
- `x`, `y`: UMAP 2D coordinates
- `family_id`: Protein family ID
- `sequence`: Representative sequence
- `size`: Number of enzymes in the family

## Workflow

### 1. UMAP Visualization (Jupyter Notebook)

The main analysis happens in `notebooks/umap_vis.ipynb`:

```python
# Load embeddings
embeds = np.load('data/embeddings.npz')

# Reduce to 2D using UMAP
reducer = umap.UMAP()
embedding_2d = reducer.fit_transform(embeds['embeddings'])

# Create visualization DataFrame
df = pd.DataFrame(embedding_2d, columns=['x', 'y'])

# Join with representative sequences
representatives_df = pd.read_csv('data/representatives.csv')
df["family_id"] = embeds["family_ids"]
df = df.merge(representatives_df[["family_id", "sequence"]], on="family_id", how="left")

# Add family size from database
# ... query enzyme_taxonomy table ...
df["size"] = df["family_id"].map(family_size_map)

# Save for downstream use
df.to_csv('data/families.csv', index=False)
```

### 2. Qdrant Vector Upload

The `scripts/upload_to_qdrant.py` script uploads embeddings to Qdrant for semantic search:

- Loads embeddings and metadata
- Creates a Qdrant collection with cosine distance metric
- Uploads points in batches to handle large payloads
- Each point includes: embedding vector, family_id, and representative sequence

## Requirements

- Python 3.12+

## Setup

1. **Clone and create virtual environment:**
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Ensure Qdrant is running:**
   ```bash
   # Qdrant should be running on localhost:6333, see documentation at https://qdrant.tech/documentation/quickstart/ for setting it up
   ```

## Usage

### Generate Visualizations & Analysis

Run the Jupyter notebook:
```bash
cd esm-landscape
jupyter notebook notebooks/umap_vis.ipynb
```

This will:
- Generate UMAP coordinates
- Create visualization plot in `plots/first_plot.png`
- Save family data to `data/families.csv`

### Upload to Qdrant

```bash
cd esm-landscape/scripts
python upload_to_qdrant.py
```

This will:
- Create a `esm_embeddings` collection in Qdrant
- Upload all embeddings with metadata in batches
- Print progress for each batch

## Notes

- ESM embeddings are 640-dimensional
- Qdrant uses cosine distance for similarity search
- Family sizes are cached from database at notebook runtime