import argparse
import torch
import pandas as pd
import numpy as np
from transformers import AutoModelForMaskedLM, AutoTokenizer

BATCH_SIZE = 16
MAX_LENGTH = 1024

parser = argparse.ArgumentParser(description="Generate ESMC embeddings from a CSV of sequences.")
parser.add_argument("input", help="Path to input CSV file (must have a 'sequence' column)")
parser.add_argument("output", help="Path for output embeddings .npy file")
parser.add_argument("--index", help="Path for output index CSV (default: <output_stem>_index.csv)")
args = parser.parse_args()

CSV_PATH = args.input
OUTPUT_PATH = args.output
INDEX_PATH = args.index if args.index else OUTPUT_PATH.replace(".npy", "_index.csv")

model = AutoModelForMaskedLM.from_pretrained(
    "biohub/ESMC-600M",
    torch_dtype=torch.bfloat16,
).eval().to("cuda")
tokenizer = AutoTokenizer.from_pretrained("biohub/ESMC-600M")

df = pd.read_csv(CSV_PATH)
sequences = df["sequence"].tolist()
n = len(sequences)

# Infer hidden size from a single probe pass (config attribute name varies by model)
with torch.inference_mode():
    probe = tokenizer("A", return_tensors="pt")
    probe = {k: v.to(model.device) for k, v in probe.items()}
    probe_out = model(**probe, output_hidden_states=True)
    hidden_size = probe_out.hidden_states[-1].shape[-1]

print(f"Loaded {n} sequences, hidden_size={hidden_size}")

all_embeddings = np.zeros((n, hidden_size), dtype=np.float32)

for i in range(0, n, BATCH_SIZE):
    batch = sequences[i : i + BATCH_SIZE]
    inputs = tokenizer(
        batch,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=MAX_LENGTH,
    )
    inputs = {k: v.to(model.device) for k, v in inputs.items()}

    with torch.inference_mode():
        output = model(**inputs, output_hidden_states=True)

    # Mean-pool last hidden state over non-padding tokens
    last_hidden = output.hidden_states[-1].float()   # (B, L, H)
    mask = inputs["attention_mask"].unsqueeze(-1).float()  # (B, L, 1)
    embeddings = (last_hidden * mask).sum(1) / mask.sum(1).clamp(min=1e-9)

    all_embeddings[i : i + len(batch)] = embeddings.cpu().numpy()

    if i % (BATCH_SIZE * 50) == 0:
        print(f"  {i + len(batch)}/{n}")

np.save(OUTPUT_PATH, all_embeddings)
df[["30pid_superfamily_id", "centroid_orf_id"]].to_csv(INDEX_PATH, index=False)
print(f"Done. embeddings shape: {all_embeddings.shape}")
print(f"Saved to {OUTPUT_PATH} and {INDEX_PATH}")