import time
import sys
import requests
import pandas as pd

UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/accessions"
BATCH_SIZE = 200  # UniProt recommends ≤200 IDs per POST request
SEED = 471829


def fetch_batch(accessions: list[str]) -> list[dict]:
    resp = requests.get(
        UNIPROT_REST,
        params={
            "accessions": ",".join(accessions),
            "format": "json",
            "fields": "accession,sequence",
        },
        timeout=60,
    )
    resp.raise_for_status()
    return resp.json().get("results", [])


def fetch_all(ids: list[str]) -> list[dict]:
    results = []
    total = len(ids)
    for i in range(0, total, BATCH_SIZE):
        batch = ids[i:i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        total_batches = (total + BATCH_SIZE - 1) // BATCH_SIZE
        print(f"  Batch {batch_num}/{total_batches} ({len(results) + len(batch)}/{total})...", flush=True)
        try:
            results.extend(fetch_batch(batch))
        except requests.RequestException as e:
            print(f"  Warning: batch {batch_num} failed ({e}), skipping.", flush=True)
        time.sleep(0.2)
    return results


df = pd.read_csv("data/family-representatives.csv")
n = len(df)

ids_df = pd.read_csv("data/uniprot_ids.tsv", sep="\t")
sampled_ids = ids_df["Entry"].sample(n=n, random_state=SEED).tolist()

print(f"Fetching {n} random UniProt sequences...")
entries = fetch_all(sampled_ids)

if not entries:
    print("No proteins retrieved. Exiting.", file=sys.stderr)
    sys.exit(1)

records = [
    {"family_id": f"rand_uniprot_{e['primaryAccession']}", "sequence": e["sequence"]["value"]}
    for e in entries
    if e.get("sequence", {}).get("value")
]

out = pd.DataFrame(records)
out.to_csv("controls/rand_uniprot.csv", index=False)
print(f"Written {len(out):,} sequences → controls/rand_uniprot.csv  (seed={SEED})")
