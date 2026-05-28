import pandas as pd
import numpy as np

SEED = 471829

rng = np.random.default_rng(SEED)

df = pd.read_csv("data/family-representatives.csv")

records = []
for _, row in df.iterrows():
    fam_id = row["family_id"]
    seq = list(row["sequence"])
    rng.shuffle(seq)
    records.append({"family_id": f"shuffled_{fam_id}", "sequence": "".join(seq)})

out = pd.DataFrame(records)
out.to_csv("controls/shuffled.csv", index=False)
print(f"Written {len(out):,} sequences → controls/shuffled.csv  (seed={SEED})")