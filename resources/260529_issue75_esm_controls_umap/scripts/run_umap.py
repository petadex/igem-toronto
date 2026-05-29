"""
Generate a 2D UMAP from one or more .npz embedding files and save the result.

Supported .npz formats:
  - family_embeddings.npz : keys family_ids + embeddings  → labeled by --labels arg
  - controls.npz          : keys accession_ids + embeddings + variants → labeled by variants column

Usage:
    python umap.py family_embeddings.npz controls.npz --output umap_output.csv --plot umap.png
"""

import argparse
import numpy as np
import pandas as pd
from umap import UMAP
import matplotlib.pyplot as plt

UMAP_PARAMS = dict(
    n_neighbors=15,
    n_components=2,
    metric="euclidean",
    min_dist=0.1,
    random_state=42,
)


def load_embeddings(path, label):
    data = np.load(path, allow_pickle=True)
    embeds = data["embeddings"].astype(np.float32)

    if "variants" in data.files:
        labels = data["variants"]
        ids = data["accession_ids"]
    else:
        labels = np.full(len(embeds), label)
        ids = data["family_ids"]

    df = pd.DataFrame({"id": ids, "label": labels})
    return df, embeds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("embeddings", nargs="+", help="Path(s) to .npz embedding files")
    parser.add_argument("--labels", nargs="*", help="Label per file (only used for files without a variants column)")
    parser.add_argument("--output", default="umap_output.csv", help="Output CSV path")
    parser.add_argument("--plot", default="umap.png", help="Output plot path")
    args = parser.parse_args()

    stems = [p.split("/")[-1].replace(".npz", "") for p in args.embeddings]
    provided = args.labels or []
    fallback_labels = [provided[i] if i < len(provided) else stems[i] for i in range(len(args.embeddings))]

    frames, embed_blocks = [], []
    for path, label in zip(args.embeddings, fallback_labels):
        df, embeds = load_embeddings(path, label)
        frames.append(df)
        embed_blocks.append(embeds)

    meta = pd.concat(frames, ignore_index=True)
    all_embeddings = np.vstack(embed_blocks)

    print(f"Fitting UMAP on {len(all_embeddings):,} sequences ({all_embeddings.shape[1]}D) ...")
    print(f"Labels: {sorted(meta['label'].unique())}")
    reducer = UMAP(**UMAP_PARAMS)
    coords = reducer.fit_transform(all_embeddings)

    meta["x"] = coords[:, 0]
    meta["y"] = coords[:, 1]
    meta.to_csv(args.output, index=False)
    print(f"Saved coordinates → {args.output}")

    fig, ax = plt.subplots(figsize=(12, 10))
    for label, group in meta.groupby("label"):
        ax.scatter(group["x"], group["y"], s=3, alpha=0.4, label=label)
    ax.set_title("UMAP Projection of Protein Embeddings")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=4, loc="best")
    fig.savefig(args.plot, dpi=300, bbox_inches="tight")
    print(f"Saved plot → {args.plot}")


if __name__ == "__main__":
    main()
