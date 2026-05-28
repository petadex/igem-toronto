import pandas as pd
import numpy as np

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
SEED = 471829

rng = np.random.default_rng(SEED)

df = pd.read_csv("data/family-representatives.csv")
lengths = df["sequence"].str.len().values

p5  = int(np.percentile(lengths, 5))
p95 = int(np.percentile(lengths, 95))
print(f"Length bounds (5th–95th percentile): {p5}–{p95} aa")

n = len(df)
sampled_lengths = rng.integers(p5, p95 + 1, size=n)

records = []
for i, length in enumerate(sampled_lengths):
    seq = "".join(rng.choice(AMINO_ACIDS, size=length))
    records.append({"family_id": f"rand_95th_family_{i}", "sequence": seq})

out = pd.DataFrame(records)
out.to_csv("controls/rand_95th.csv", index=False)
print(f"Written {len(out):,} sequences → controls/rand_95th.csv  (seed={SEED})")