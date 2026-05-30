import sys
import numpy as np


def summarize_npz(path):
    data = np.load(path, allow_pickle=True)
    keys = list(data.files)

    print(f"File: {path}")
    print(f"Arrays stored: {keys}\n")

    for key in keys:
        arr = data[key]
        print(f"  [{key}]")
        print(f"    shape : {arr.shape}")
        print(f"    dtype : {arr.dtype}")

        # For object arrays (variable-length embeddings), report per-item shapes
        if arr.dtype == object:
            shapes = [a.shape for a in arr]
            unique_shapes = set(shapes)
            if len(unique_shapes) == 1:
                print(f"    item shape : {shapes[0]} (uniform)")
            else:
                print(f"    item shapes: variable — e.g. {shapes[:3]}")
        elif np.issubdtype(arr.dtype, np.floating) or np.issubdtype(arr.dtype, np.integer):
            print(f"    min/max : {arr.min():.4g} / {arr.max():.4g}")
        print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python npy_shape.py <file.npz>")
        sys.exit(1)
    summarize_npz(sys.argv[1])
