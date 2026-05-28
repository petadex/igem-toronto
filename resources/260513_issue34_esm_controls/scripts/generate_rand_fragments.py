import pandas as pd
import numpy as np

df = pd.read_csv("data/family-representatives.csv")
SEED = 471829
rng = np.random.default_rng(SEED)

records = []
for _, row in df.iterrows():
    seq = row["sequence"]
    # generates 3 fragments of length 30%, 60%, and 90% of original
    for frac in [0.3, 0.6, 0.9]:
        frag_len = int(len(seq) * frac)
        rand_start = rng.choice(len(seq) - frag_len + 1)
        rand_end = rand_start + frag_len
        records.append({"family_id": f"rand_fragment_{row['family_id']}", "sequence": seq[rand_start:rand_end]})

out = pd.DataFrame(records)
out.to_csv("controls/rand_fragments.csv", index=False)
print(f"Written {len(out):,} sequences → controls/rand_fragments.csv  (seed={SEED})")
