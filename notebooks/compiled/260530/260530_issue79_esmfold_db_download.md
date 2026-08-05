## ESMFold2 Structures Download

Lead     : `AlexLeonardos`

Issue    : [Github Issue #79](https://github.com/petadex/igem-toronto/issues/) — _ESMFold2 Structure Download_

Start    : `2026-05-28`

Complete : `2026-06-XX`

Files    : `~/igem-toronto/files/YYMMDD_issue##_<slug>/` — local working directory (active analysis, intermediate outputs)

---

Many folded structures for protein hashes are found on s3://esm-protein-atlas/v1/folds/

The purpose of these downloads is to avoid rerunning folding experiments on sequences that have already been folded.

### Data Accessed: 2026-05-30
```bash
# commands used to pull input data
aws s3 cp s3://petadex/<path>/<file> ./
```

## Step 1: Figure out the File Structure on S3

```
esmfold@esmfold-structure-downloader:~$ aws s3 ls s3://esm-protein-atlas/v1/folds/ --no-sign-request
                           PRE folds_1B.lance/
                           PRE folds_atlas.lance/
esmfold@esmfold-structure-downloader:~$ aws s3 ls s3://esm-protein-atlas/v1/folds/folds_1B.lance/ --no-sign-request
                           PRE _indices/
                           PRE _transactions/
                           PRE _versions/
                           PRE data/
esmfold@esmfold-structure-downloader:~$ aws s3 ls s3://esm-protein-atlas/v1/folds/folds_atlas.lance/ --no-sign-request
                           PRE _indices/
                           PRE _transactions/
                           PRE _versions/
                           PRE data/
esmfold@esmfold-structure-downloader:~$ aws s3 ls s3://esm-protein-atlas/v1/folds/folds_1B.lance/data/ --no-sign-request | wc -l
109554
esmfold@esmfold-structure-downloader:~$ aws s3 ls s3://esm-protein-atlas/v1/folds/folds_atlas.lance/data/ --no-sign-request | wc -l
667
esmfold@esmfold-structure-downloader:~$ aws s3 ls s3://esm-protein-atlas/v1/folds/folds_atlas.lance/data/ --no-sign-request | awk '{sum += $3} END {print sum/1024/1024/1024 " GB"}'
441.409 GB

aws s3 ls s3://esm-protein-atlas/v1/folds/folds_1B.lance/data/ --human-readable --summarize
Total Objects: 109554
   Total Size: 68.3 TiB
```

The indices folder seems to contain a lookup table, so this should be inspected first:

aws s3 cp s3://esm-protein-atlas/v1/folds/folds_1B.lance/_indices/29097da6-6d4e-4e15-8818-28e4926544f5/page_lookup.lance ./page_lookup.lance --no-sign-request

Schema for the versions and the indices:
```
python3 ~/check_schema.py
header: string
  -- field metadata --
  compression: 'zstd'
protein_hash: string
  -- field metadata --
  compression: 'zstd'
ptm: double
mean_plddt: double
per_residue_plddt: binary
  -- field metadata --
  compression: 'zstd'
pae: binary
  -- field metadata --
  compression: 'zstd'
structure_blob: binary
  -- field metadata --
  compression: 'zstd'
sequence: large_string
  -- field metadata --
  compression: 'zstd'
1095530880
```

Now that we know the schema, we can be sure of a next step: hash all of the proteins, use that as a query and pull them from the lance files.



```python
# B Trees for hash and accession:
esmfold@esmfold-structure-downloader:~$ python3 -c "
import lance
ds = lance.dataset('/home/esmfold/folds_1B.lance')
for idx in ds.list_indices():
    print({k: v for k, v in idx.items() if k != 'fragment_ids'})
"
{'name': 'protein_hash_idx', 'type': 'BTree', 'uuid': 'fdb71a9d-37fb-4edb-b1d1-5611914f9bb8', 'fields': ['protein_hash'], 'version': 1, 'base_id': None}
{'name': 'header_idx', 'type': 'BTree', 'uuid': '29097da6-6d4e-4e15-8818-28e4926544f5', 'fields': ['header'], 'version': 2, 'base_id': None}
```

# Plan:

1. Compute MD5 hashes for all 4.7M sequences → list of protein_hash values
2. Query folds_1B.lance directly filtering by those hashes using the  → get structure_blob
3. Decode and save as PDB files

This work is done in a script found in the resources folder.

## Steps Example:

This example is for one hash that we know: 

```
esmfold@esmfold-structure-downloader:~$ AWS_EC2_METADATA_DISABLED=true python3 -c "
import lance

ds = lance.dataset(
    's3://esm-protein-atlas/v1/folds/folds_1B.lance',
    storage_options={
        'skip_signature': 'true',
        'aws_region': 'us-west-2'
    }
)

target = '0000000c7d30e82e65f958222be96bf4'
result = ds.to_table(
    filter=f\"protein_hash = '{target}'\",
    columns=['protein_hash', 'mean_plddt', 'ptm', 'structure_blob']
).to_pandas()

print(result[['protein_hash', 'mean_plddt', 'ptm']])
print('structure_blob size:', len(result['structure_blob'][0]), 'bytes')

with open('./test_structure.bin', 'wb') as f:
    f.write(bytes(result['structure_blob'][0]))
print('Saved.')
"
                       protein_hash  mean_plddt       ptm
0  0000000c7d30e82e65f958222be96bf4    0.708267  0.546509
structure_blob size: 7033 bytes
Saved.
```

#### Note: Work on VM

1. Count how many of our sequences exist on the ESMAtlas.


parse: 100%|████████████████████████████████████████████████████████████████████████| 109G/109G [37:18<00:00, 52.1MB/s]
123966734 unique sequences parsed from s3://petadex/logan/petadex.catalytic_orfs.v1.1.fa
full orf_id -> protein_hash index written to orf_map.tsv.zst
count folds_atlas.lance:  42%|██████▋         | 13010/30992 [54:18<1:47:43,  2.78batch/s, hits=1765, unmatched=1.24e+8]

This run got stopped so it's not very useful.

Hash for Random Million ORFs:
HIT RATE: 5068/1000000 = 0.51% present in Atlas
       5068 from folds_1B.lance
     994932 not found in any dataset
So extrapolating from just this we'll have ~630k from running the whole thing


2. Do a full run and upload to s3:

So before the download, you'll need one of:

someone to attach a role with s3:PutObject on a bucket you own, or
temporary credentials (aws_access_key_id/secret) for an account that can write, or
fall back to writing CIFs to a local/EBS path (--outdir, no --s3) — but that needs a big volume for millions of files, which we sized away.

Note:

123966734 unique hashes loaded from orf_map.tsv.zst
querying datasets: 1b

This is many more hashes than catalytic cores. Are we trying to get folds for everything? What was the point of the AF DB? I'll ask this stuff and figure it out on Saturday?

Notes for Tuesday:

Results so Far for ESM Fold Atlas:
1. Hit Rate Average for 2 Runs of 1 Million Random Hashes:

1 HIT RATE: 5068/1000000 = 0.51% present in Atlas
       5068 from folds_1B.lance
     994932 not found in any dataset

2 HIT RATE: 5108/1000000 = 0.51% present in Atlas
       5108 from folds_1B.lance
     994892 not found in any dataset

**632,230.3434 Expected Hash Structures**

2. Time Benchmark

Both runs took 20+-1m. Extrapolating to full dataset that's **~41.333 Hours**.

3. Resources:
- Very cheap VM, doesn't need GPU or data since all of the accession and the searching is streamed.

## Need to Do:
1. Set up the upload so that I can upload to the bucket directly from the process.
  - Done, ask about long-running permissions
2. Run the full process? I think I need a bigger VM though cause I'm getting OOM errors on the one that's currently up.
3. Figure out how to compare the AF and ESMFold Structures
    - Redo the AF??

4.7 Million Catalytic Cores - accession number?

vs.

123966734 unique hashes loaded from orf_map.tsv.zst
 - Can be mapped from hash to ORF IDs

 Question about the ORFs, do we still need to do the preprocessing to just get the specific core that corresponds to the PETase?
 - This matters a lot for the hashing I would think.

Still need to redo the AF DB download:

## Output: Added Metrics + S3 File Structure (updated 2026-07-07)

The downloader now writes, per unique sequence (keyed by `protein_hash = md5(sequence)`), a **three-file set** that matches the ESMFold2 predictor's output (issue #6), so Atlas and freshly-folded structures are directly comparable.

**S3 layout** — default sink `s3://petadex-protein-structures/esmatlas/` (`--s3 ''` to write locally instead):

```
s3://petadex-protein-structures/esmatlas/
├── structures/<hash>.cif     # all-atom coords (atom37 → mmCIF); per-residue pLDDT in the B-factor
├── arrays/<hash>.npz         # per_residue_plddt[L], pae[L,L], residue_index[L]
└── metrics.csv               # one row per structure — resume ledger; mirrored to S3 + seeded back on resume
```

**Added quality metrics** — earlier the download kept only pLDDT (in the CIF B-factor) and **discarded pTM and PAE**. Both are now retained:

| datum | where it lands | notes |
|---|---|---|
| per-residue pLDDT | CIF B-factor + `npz:per_residue_plddt` | 0–100; from the `structure_blob`'s `confidence` field |
| mean pLDDT | `metrics.csv:mean_plddt` | chain average |
| pTM | `metrics.csv:ptm` | from the Atlas `ptm` column |
| PAE | `npz:pae` (`[L, L]`) | from the Atlas `pae` column |
| has_pae | `metrics.csv:has_pae` | `1` if a PAE matrix was present |

`metrics.csv` columns: `protein_hash, source_dataset, seq_len, mean_plddt, ptm, has_pae`.

**Gotcha (confirmed on a live row):** per the schema above, `pae` and `per_residue_plddt` are `binary` columns — each is a serialized **`np.savez` blob** (a zip whose member is `arr.npy`), *not* a list. `_coerce_pae()` `np.load()`s those bytes → `[L, L]`. pLDDT is taken from the `structure_blob`'s `confidence` field, so only the `pae` column needs this special decode.

**Verified end-to-end** on `test.fa` (`MNAKWITPVVTAFDAEGHVD`, 20 aa): reconstructed sequence matches, CIF B-factor == `npz` pLDDT (mean 71.33 == `metrics.csv`), PAE decodes to `(20, 20)` with the expected near-zero diagonal rising off-diagonal, `has_pae=1`.
