import pandas as pd
import random

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

df = pd.read_csv("data/family-representatives.csv")

records = []
for _, row in df.iterrows():
    fam_id = row["family_id"]
    seq = row["sequence"]

    rand_seq = "".join(random.choices(AMINO_ACIDS, k=len(seq)))
    records.append({"family_id": f"rand_matched_{fam_id}", "sequence": rand_seq})

out = pd.DataFrame(records)
out.to_csv("controls/rand_matched.csv", index=False)
print(f"Written {len(out):,} sequences → controls/rand_matched.csv")