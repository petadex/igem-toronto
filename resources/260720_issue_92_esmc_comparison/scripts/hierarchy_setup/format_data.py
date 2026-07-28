import psycopg2
import pandas as pd

# DB Creds

DB = {
    "host":     "petadex.ccz9y6yshbls.us-east-1.rds.amazonaws.com",
    "port":     5432,
    "database": "petadex",
    "user":     "readonly_user",
    "password": "petadex",
}

conn = psycopg2.connect(**DB)
cur  = conn.cursor()

# Helper function to fetch sequences from a FASTA file given a list of target IDs
def fetch_sequences(fasta_path, target_ids):
    target_ids = set(str(x) for x in target_ids)
    remaining = set(target_ids)
    sequences = {}
    with open(fasta_path) as f:
        orf_id = None
        seq_parts = []
        for line in f:
            if not remaining:
                break
            line = line.strip()
            if line.startswith('>'):
                if orf_id and orf_id in target_ids:
                    sequences[orf_id] = "".join(seq_parts)
                    remaining.discard(orf_id)
                orf_id = line[1:].split('|')[0]
                seq_parts = []
            elif orf_id in target_ids:
                seq_parts.append(line)
        if orf_id and orf_id in target_ids:
            sequences[orf_id] = "".join(seq_parts)
    return sequences


########## Get all 60% identity family clusters #################

cur.execute('SELECT * FROM "60pid_family_clusters";')
rows = cur.fetchall()
df = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])
num_clusters = max(df["60pid_family_id"])
sequences = fetch_sequences("./data/petadex.catalytic_orfs.v1.1.fa", df["centroid_orf_id"].tolist())
df["sequence"] = df["centroid_orf_id"].map(lambda x: sequences.get(str(x)))
print(df[["centroid_orf_id", "sequence"]].head())

# Save clusters with sequences to a parquet -> will be used for future embedding, and also avoids refetching
df.to_parquet("./data/60pid_family_clusters_with_sequences.parquet", index=False)
ids = df["centroid_orf_id"].tolist()

# Get clustering info for all centroids of the 60% clusters
cur.execute('SELECT * FROM "petadex_clustering" WHERE "orf_id" = ANY(%s);', (ids,))
rows = cur.fetchall()
df_clustering = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])
df_clustering.to_parquet("./data/60pid_centroids.parquet", index=False)

#################################################################



########## Get all 30% identity family clusters #################

cur.execute('SELECT * FROM "30pid_superfamily_clusters";')
rows = cur.fetchall()
df = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])
num_clusters = max(df["30pid_superfamily_id"])
sequences = fetch_sequences("./data/petadex.catalytic_orfs.v1.1.fa", df["centroid_orf_id"].tolist())
df["sequence"] = df["centroid_orf_id"].map(lambda x: sequences.get(str(x)))
print(df[["centroid_orf_id", "sequence"]].head())

# Save clusters with sequences to a parquet -> will be used for future embedding, and also avoids refetching
df.to_parquet("./data/30pid_superfamily_clusters_with_sequences.parquet", index=False)
ids = df["centroid_orf_id"].tolist()

# Get clustering info for all centroids of the 30% clusters
cur.execute('SELECT * FROM "petadex_clustering" WHERE "orf_id" = ANY(%s);', (ids,))
rows = cur.fetchall()
df_clustering = pd.DataFrame(rows, columns=[desc[0] for desc in cur.description])
df_clustering.to_parquet("./data/30pid_centroids.parquet", index=False)

#################################################################