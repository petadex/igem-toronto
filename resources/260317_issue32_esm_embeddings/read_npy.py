import sys
import numpy as np


def describe_array(key, arr, n_samples=3):
    print(f"  [{key}]")
    print(f"    shape : {arr.shape}  |  dtype: {arr.dtype}")

    # String / name arrays
    if arr.dtype.kind in ("U", "S", "O") and arr.ndim == 1:
        try:
            samples = [str(x) for x in arr[:n_samples]]
            print(f"    count  : {len(arr)} entries")
            print(f"    sample : {samples}" + (" ..." if len(arr) > n_samples else ""))
            return
        except Exception:
            pass

    # Object arrays (variable-length per-residue embeddings)
    if arr.dtype == object:
        lengths = [a.shape[0] for a in arr]
        emb_dim = arr[0].shape[-1] if arr[0].ndim > 1 else None
        print(f"    count       : {len(arr)} sequences")
        print(f"    seq lengths : min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.1f} residues")
        if emb_dim:
            print(f"    emb dim     : {emb_dim}")
        # Stats across all embeddings flattened
        flat = np.concatenate([a.flatten().astype(np.float32) for a in arr])
        print(f"    value range : min={flat.min():.4f}, max={flat.max():.4f}")
        print(f"    mean ± std  : {flat.mean():.4f} ± {flat.std():.4f}")
        print(f"    sample[0]   : shape={arr[0].shape}, first values={arr[0].flatten()[:5]}")
        return

    # Dense numeric arrays (mean-pooled embeddings: shape [N, D])
    if np.issubdtype(arr.dtype, np.floating) or np.issubdtype(arr.dtype, np.integer):
        if arr.ndim == 2:
            print(f"    count   : {arr.shape[0]} sequences")
            print(f"    emb dim : {arr.shape[1]}")
        print(f"    value range : min={arr.min():.4f}, max={arr.max():.4f}")
        print(f"    mean ± std  : {arr.mean():.4f} ± {arr.std():.4f}")
        if arr.ndim >= 2:
            print(f"    sample[0]   : {arr[0, :5]}")
        return

    # Fallback
    print(f"    (no detailed summary for dtype={arr.dtype})")


def inspect_npz(path, n_samples=3):
    data = np.load(path, allow_pickle=True)
    keys = list(data.files)

    print(f"=== {path} ===")
    print(f"Arrays: {keys}")
    print()

    for key in keys:
        describe_array(key, data[key], n_samples=n_samples)
        print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python read_npy.py <file.npz> [n_samples]")
        sys.exit(1)
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    inspect_npz(sys.argv[1], n_samples=n)
