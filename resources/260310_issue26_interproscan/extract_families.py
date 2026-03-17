"""
extract_families.py
-------------------
Extracts per-family FASTA files from the PETadex database.
Tables used:
  - enzyme_taxonomy   : enzyme_id → family (int), family_pid, component
  - enzyme_fastaa     : enzyme_id → translated_sequence

Family size range: SIZE_MIN to SIZE_MAX (inclusive).
Output: output/<family_id>/sequences.fasta
"""

import os
import psycopg2

# ── Config ─────────────────────────────────────────────────────────────────
DB = {
    "host":     "petadex.ccz9y6yshbls.us-east-1.rds.amazonaws.com",
    "port":     5432,
    "database": "petadex",
    "user":     "readonly_user",
    "password": "petadex",
}
SIZE_MIN   = 4
SIZE_MAX   = 500
OUTPUT_DIR = "output"
# ───────────────────────────────────────────────────────────────────────────

conn = psycopg2.connect(**DB)
cur  = conn.cursor()

# ── Step 1: Get all families within the size range ─────────────────────────
print("Fetching family list...")
cur.execute("""
    SELECT family, COUNT(DISTINCT enzyme_id) AS enzyme_count
    FROM enzyme_taxonomy
    WHERE family IS NOT NULL
    GROUP BY family
    HAVING COUNT(DISTINCT enzyme_id) BETWEEN %s AND %s
    ORDER BY enzyme_count DESC;
""", (SIZE_MIN, SIZE_MAX))

families = cur.fetchall()
print(f"Families to process: {len(families):,}")

# ── Step 2: For each family, fetch sequences and write FASTA ───────────────
for i, (family_id, count) in enumerate(families, 1):
    family_dir = os.path.join(OUTPUT_DIR, str(family_id))
    fasta_path = os.path.join(family_dir, "sequences.fasta")

    # Skip if already extracted
    if os.path.exists(fasta_path):
        continue

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
    """, (family_id,))

    rows = cur.fetchall()

    with open(fasta_path, "w") as fh:
        for enzyme_id, seq in rows:
            fh.write(f">{enzyme_id}\n{seq}\n")

    if i % 500 == 0 or i == len(families):
        print(f"  [{i:,}/{len(families):,}] family {family_id}: {len(rows)} sequences")

conn.close()
print(f"\nDone. FASTAs written to: {OUTPUT_DIR}/")