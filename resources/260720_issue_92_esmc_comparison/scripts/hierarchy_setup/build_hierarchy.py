import pandas as pd

# grab data pulled using format_data.py
centroids_df = pd.read_csv("./data/60pid_centroids.csv")
sequences_df = pd.read_csv("./data/60pid_family_clusters_with_sequences.csv")
df = centroids_df.merge(sequences_df, how="left", on="60pid_family_id")

# Build a dictionary mapping 30%pid centroids to all of the 60% centroids that make up the 30% superfamily
superfamilies_dict = {}
for index, row in df.iterrows():
    superfamily_id = str(row["30pid_superfamily_id"])
    family_stats = {
        "60pid_centroid_orf_id": str(row["orf_id"]),
        "60pid_family_id": str(row["60pid_family_id"]),
        "centroid_sequence": str(row["sequence"])
    }
    if superfamily_id not in superfamilies_dict:
        superfamilies_dict[superfamily_id] = []
    superfamilies_dict[superfamily_id].append(family_stats)

# Make a parquet for each superfamily, containing all of the 60% family centroids that make it up
keys = list(superfamilies_dict.keys())
for key in keys:
    df = pd.DataFrame(superfamilies_dict[key])
    df.to_parquet(f"./data/superfamily_hierarchy/{key}.parquet", index=False, engine="pyarrow")