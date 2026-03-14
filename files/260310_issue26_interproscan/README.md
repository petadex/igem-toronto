# InterProScan Domain Annotation Pipeline

End-to-end pipeline for generating InterProScan domain annotations for PETadex enzyme families.

---

## Files in This Folder

| File | Purpose |
|---|---|
| `extract_families.py` | Pull representative sequences for all families (size 4–500) from the PETadex database and write per-family FASTA files |
| `extract_family_182.py` | Same as above, but for a single family (182) — useful for quick tests |
| `split_fasta.sh` | Split a combined FASTA into fixed-size chunks (default: 5,000 sequences per chunk) |
| `run_iprscan.sh` | Run InterProScan over all chunks, skip already-completed chunks, and merge results into a single TSV |
| `interproscan_amazon_linux2.md` | Full installation and configuration guide for InterProScan on Amazon Linux 2 |

---

## Prerequisites

- InterProScan installed at `~/interproscan/interproscan-5.72-103.0/` (see `interproscan_amazon_linux2.md`)
- Python 3 with `psycopg2` (`pip install psycopg2-binary`)
- Access to the PETadex RDS instance

---

## Step-by-Step: Generating InterProScan Domains

### Step 1 — Extract sequences from PETadex

Run from the repo root (or any working directory — output goes to `./output/`):

```bash
python3 interproscan_pipeline/extract_families.py
```

This queries the PETadex database and writes one FASTA per family:

```
output/
  182/sequences.fasta
  204/sequences.fasta
  ...
```

To extract a single family for testing:

```bash
python3 interproscan_pipeline/extract_family_182.py
```

---

### Step 2 — Build a combined representative FASTA

Concatenate the per-family sequences into a single file (or use your existing representative FASTA directly):

```bash
cat output/*/sequences.fasta > ~/petadex-family-representatives.fasta
```

If you already have `~/petadex-family-representatives.fasta` from a prior run, skip this step.

---

### Step 3 — Split into chunks

InterProScan performs better on smaller input files. Split the combined FASTA into chunks of 5,000 sequences:

```bash
bash interproscan_pipeline/split_fasta.sh \
    ~/petadex-family-representatives.fasta \
    ~/iprscan_chunks \
    5000
```

Chunks are written to `~/iprscan_chunks/` as `chunk_0001.fasta`, `chunk_0002.fasta`, etc.

---

### Step 4 — Run InterProScan

```bash
bash interproscan_pipeline/run_iprscan.sh
```

Key settings at the top of `run_iprscan.sh` — edit before running if needed:

| Variable | Default | Meaning |
|---|---|---|
| `IPRSCAN` | `~/interproscan/interproscan-5.72-103.0/interproscan.sh` | Path to InterProScan executable |
| `CHUNKS_DIR` | `~/iprscan_chunks` | Directory of input chunk FASTAs |
| `RESULTS_DIR` | `~/iprscan_results` | Where TSV results are written |
| `APPL` | `Pfam,TIGRFAM,CDD,SUPERFAMILY,Gene3D,PANTHER` | Databases to search |
| `CPU` | `8` | Threads per run |
| `START_CHUNK` | `1` | First chunk to process (1-indexed) |
| `END_CHUNK` | `13` | Last chunk to process (inclusive) |

The script skips chunks whose output already exists, so it is safe to re-run after interruption.

**Run in a tmux session** to survive SSH disconnection:

```bash
tmux new -s iprscan
bash interproscan_pipeline/run_iprscan.sh 2>&1 | tee ~/iprscan_run.log
```

---

### Step 5 — Inspect the merged output

When all chunks complete, the script automatically merges results:

```
~/iprscan_results/all_results.tsv
```

TSV columns (one domain hit per row):

```
protein_id  md5  length  analysis  accession  description  start  stop  evalue  status  date  ipr_accession  ipr_description  go_terms
```

Quick summary:

```bash
wc -l ~/iprscan_results/all_results.tsv
cut -f4 ~/iprscan_results/all_results.tsv | sort | uniq -c | sort -rn   # hits per database
```

---

## Notes

- Full runs over tens of thousands of sequences take several hours. Monitor progress in `~/iprscan_run.log`.
- If you need to reprocess specific chunks, adjust `START_CHUNK` / `END_CHUNK` in `run_iprscan.sh`.
- See `interproscan_amazon_linux2.md` for troubleshooting installation issues (Java, libz, disk space, etc.).
