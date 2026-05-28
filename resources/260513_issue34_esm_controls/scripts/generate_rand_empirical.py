import pandas as pd
import numpy as np

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
SEED = 471829

rng = np.random.default_rng(SEED)

df = pd.read_csv("data/family-representatives.csv")
lengths = df["sequence"].str.len().values

n = len(df)
sampled_lengths = rng.choice(lengths, size=n, replace=True)

records = []
for i, length in enumerate(sampled_lengths):
    seq = "".join(rng.choice(AMINO_ACIDS, size=length))
    records.append({"family_id": f"rand_empirical_family_{i}", "sequence": seq})

out = pd.DataFrame(records)
out.to_csv("controls/rand_empirical.csv", index=False)
print(f"Written {len(out):,} sequences → controls/rand_empirical.csv  (seed={SEED})")