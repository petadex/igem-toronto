"""Extract FASTA for family 182 only."""

import os
import psycopg2

DB = {
    "host":     "petadex.ccz9y6yshbls.us-east-1.rds.amazonaws.com",
    "port":     5432,
    "database": "petadex",
    "user":     "readonly_user",
    "password": "petadex",
}

FAMILY_ID  = 182
OUTPUT_DIR = "output"

conn = psycopg2.connect(**DB)
cur  = conn.cursor()

family_dir = os.path.join(OUTPUT_DIR, str(FAMILY_ID))
fasta_path = os.path.join(family_dir, "sequences.fasta")
os.makedirs(family_dir, exist_ok=True)

cur.execute("""
    SELECT
        et.enzyme_id,
        ef.translated_sequence
    FROM enzyme_taxonomy et
    JOIN enzyme_fastaa ef ON et.enzyme_id = ef.enzyme_id
    WHERE et.family = %s
      AND ef.translated_sequence IS NOT NULL
      AND ef.translated_sequence != '';
""", (FAMILY_ID,))

rows = cur.fetchall()

with open(fasta_path, "w") as fh:
    for enzyme_id, seq in rows:
        fh.write(f">{enzyme_id}\n{seq}\n")

conn.close()
print(f"Written {len(rows)} sequences to {fasta_path}")
