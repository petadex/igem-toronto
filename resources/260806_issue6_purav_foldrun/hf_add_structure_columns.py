#!/usr/bin/env python
"""Part III: add is_c60_centroid and structure_path to the orfs/ spine of petadex/petadex-corpus-full,
and open a Hugging Face pull request with the rewritten parts and dataset card.

    is_c60_centroid  true for the 1,814,296 orf_ids in dims/clusters60.centroid_orf_id, else false
    structure_path   full S3 URI of the centroid's structure (data/fold60pid/structure_index.tsv),
                     null for non-centroids and for the 1,054 centroids over the 2,500 aa cap

Both go directly after is_c90_centroid; sequence stays last. Rows keep their exact original order (the parts
are sorted by orf_id only within each row group), and every other column is unchanged: checked per part by
position (row i has the same orf_id as before) and by a row hash over all original columns.

    python hf_add_structure_columns.py lookup                 # build the 1.81 M-row lookup table
    python hf_add_structure_columns.py rewrite [--parts 0,1]  # download + rewrite + check parts
    python hf_add_structure_columns.py pr                     # upload rewritten parts + card as a PR
"""
import argparse
import json
import os
import sys

import duckdb

REPO = "petadex/petadex-corpus-full"
BASE = "d4cd275b0ca59ef409b7886aeeb706431c822668"   # main on 2026-10-19; the PR is checked against it
DST = "s3://petadex-protein-structures/petadex-v1.1-60pid-centroids-structures"
INDEX = "data/fold60pid/structure_index.tsv"
WORK = f"/scratch/{os.environ['USER']}/hf_corpus_full"
SRC, OUT = f"{WORK}/src", f"{WORK}/out"
LOOKUP = f"{WORK}/c60_structure_lookup.parquet"
NPARTS = 131
OLD = ["orf_id", "origin", "date_retrieval", "accession", "library_id", "contig", "orf_start", "orf_end",
       "orf_type", "length", "has_stop", "md5", "c90_id", "c60_id", "c30_id", "is_c90_centroid", "sequence"]


def con():
    c = duckdb.connect()
    c.sql(f"SET threads={os.environ.get('SLURM_CPUS_PER_TASK', 8)}; SET memory_limit='{int(int(os.environ.get('SLURM_MEM_PER_NODE', 16000)) * 0.75)}MB'; "
          f"SET preserve_insertion_order=true; SET temp_directory='{WORK}/tmp'")
    return c


def lookup():
    from huggingface_hub import hf_hub_download
    os.makedirs(WORK, exist_ok=True)
    c60 = hf_hub_download(REPO, "dims/clusters60.parquet", repo_type="dataset", revision=BASE, local_dir=SRC)
    c = con()
    c.sql(f"""create table L as
        select c.centroid_orf_id::bigint orf_id,
               case when i.structure_key <> '' then '{DST}/' || i.structure_key end structure_path
        from '{c60}' c left join read_csv('{INDEX}', delim='\t', header=true, all_varchar=true) i
             on i.orf_id::bigint = c.centroid_orf_id""")
    n, npath, ndup = c.sql("select count(*), count(structure_path), count(*) - count(distinct orf_id) from L").fetchone()
    nidx = c.sql(f"select count(*) from read_csv('{INDEX}', delim='\t', header=true)").fetchone()[0]
    print(f"lookup: {n:,} c60 centroids, {npath:,} with a structure_path, {ndup} duplicate ids; index rows {nidx:,}")
    assert n == nidx == 1_814_296 and ndup == 0 and npath == 1_813_242, "lookup does not match the structure set"
    c.sql(f"copy (select * from L order by orf_id) to '{LOOKUP}' (format parquet, compression zstd)")


def rewrite(parts):
    from huggingface_hub import hf_hub_download
    os.makedirs(f"{OUT}/orfs", exist_ok=True)
    os.makedirs(f"{WORK}/tmp", exist_ok=True)
    c = con()
    cols = ", ".join(f"o.{x}" for x in OLD[:-1])
    h_old = "hash(" + ", ".join(OLD) + ")"
    for k in parts:
        name = f"orfs/part-{k}.parquet"
        out = f"{OUT}/{name}"
        stat_f = f"{out}.check.json"
        if os.path.exists(stat_f):
            continue
        src = hf_hub_download(REPO, name, repo_type="dataset", revision=BASE, local_dir=SRC)
        c.sql(f"""copy (select {cols}, l.orf_id is not null as is_c60_centroid, l.structure_path, o.sequence
                        from read_parquet('{src}', file_row_number=true) o
                             left join read_parquet('{LOOKUP}') l using (orf_id)
                        order by o.file_row_number)
                  to '{out}.tmp' (format parquet, compression zstd, row_group_size 501760)""")
        # checks: same rows in the same positions, original columns byte-identical, flags consistent
        a = c.sql(f"select count(*), sum({h_old})::hugeint from read_parquet('{src}')").fetchone()
        b = c.sql(f"""select count(*), sum({h_old})::hugeint, null,
                             count(*) filter (where is_c60_centroid), count(structure_path),
                             count(*) filter (where is_c60_centroid and not is_c90_centroid),
                             count(*) filter (where structure_path is not null and not is_c60_centroid)
                      from read_parquet('{out}.tmp')""").fetchone()
        moved = c.sql(f"""select count(*) from read_parquet('{src}', file_row_number=true) s
                          join read_parquet('{out}.tmp', file_row_number=true) t using (file_row_number)
                          where s.orf_id <> t.orf_id""").fetchone()[0]
        schema = [r[0] for r in c.sql(f"describe select * from read_parquet('{out}.tmp')").fetchall()]
        want = OLD[:-1] + ["is_c60_centroid", "structure_path", "sequence"]
        ok = a[0] == b[0] and a[1] == b[1] and moved == 0 and b[5] == 0 and b[6] == 0 and schema == want
        s = dict(part=k, rows=a[0], rows_moved=moved, hash_equal=a[1] == b[1], is_c60=b[3], with_path=b[4],
                 c60_not_c90=b[5], path_not_c60=b[6], schema_ok=schema == want, ok=ok,
                 src_bytes=os.path.getsize(src), out_bytes=os.path.getsize(f"{out}.tmp"))
        print(json.dumps(s), flush=True)
        if not ok:
            sys.exit(f"part {k} failed its checks; {out}.tmp kept for inspection")
        os.replace(f"{out}.tmp", out)
        json.dump(s, open(stat_f, "w"))


def totals():
    rows = [json.load(open(f"{OUT}/orfs/part-{k}.parquet.check.json")) for k in range(NPARTS)]
    t = {x: sum(r[x] for r in rows) for x in ("rows", "is_c60", "with_path", "c60_not_c90", "path_not_c60",
                                               "src_bytes", "out_bytes")}
    t["all_ok"] = all(r["ok"] for r in rows)
    print(json.dumps(t))
    assert t["all_ok"] and t["rows"] == 307_155_746 and t["is_c60"] == 1_814_296 and t["with_path"] == 1_813_242
    return t


CARD_OLD = "`c90_id`, `c60_id`, `c30_id`, `is_c90_centroid`, `sequence`."
CARD_NEW = """`c90_id`, `c60_id`, `c30_id`, `is_c90_centroid`, `is_c60_centroid`, `structure_path`, `sequence`.

`is_c60_centroid` is true for the 1,814,296 centroids of the 60% clusters (`dims/clusters60.centroid_orf_id`);
every 60% centroid is also a 90% centroid. `structure_path` is the S3 URI of the predicted structure for
1,813,242 of them (gzip mmCIF, pLDDT in the B-factor column), from
`s3://petadex-protein-structures/petadex-v1.1-60pid-centroids-structures/`: ESMFold2 (Full under 150 aa,
Fast 150-2,500 aa) or, where the centroid's exact sequence is in AlphaFold DB, the AFDB v4 model. It is null for
every other ORF and for the 1,054 centroids longer than 2,500 aa. That folder's `README.md` and
`structure_index.tsv` give the model and metrics per structure (iGEM Toronto 2026, issue #6)."""


def pr():
    from huggingface_hub import CommitOperationAdd, HfApi, hf_hub_download
    t = totals()
    api = HfApi()
    card = open(hf_hub_download(REPO, "README.md", repo_type="dataset", revision=BASE)).read()
    assert card.count(CARD_OLD) == 1
    os.makedirs(OUT, exist_ok=True)
    open(f"{OUT}/README.md", "w").write(card.replace(CARD_OLD, CARD_NEW))
    state_f = f"{WORK}/pr.json"
    st = json.load(open(state_f)) if os.path.exists(state_f) else {"done": []}
    desc = f"""Adds two columns to the `orfs/` spine, directly after `is_c90_centroid`:

- `is_c60_centroid` (bool): true for the {t['is_c60']:,} orf_ids in `dims/clusters60.centroid_orf_id`.
- `structure_path` (string): S3 URI of the centroid's predicted structure, set for {t['with_path']:,} rows;
  null elsewhere (including the 1,054 centroids over 2,500 aa, which were not folded).

Structures: `{DST}/` (ESMFold2 Full/Fast, 32 steps, 3 loops, plus AFDB v4 models where the centroid's exact
sequence is in AFDB). Every other column is unchanged: each part was checked for identical row count, `orf_id`
order and an order-independent hash over all original columns. All {NPARTS} parts are rewritten
({t['rows']:,} rows). The dataset card documents the new columns.

Made by `resources/260806_issue6_purav_foldrun/hf_add_structure_columns.py` in petadex/igem-toronto
(notebook `260806_issue6_purav_foldrun.ipynb`, Part III)."""
    files = ["README.md"] + [f"orfs/part-{k}.parquet" for k in range(NPARTS)]
    todo = [f for f in files if f not in st["done"]]
    batch = 10
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        ops = [CommitOperationAdd(path_in_repo=f, path_or_fileobj=f"{OUT}/{f}") for f in chunk]
        msg = f"Add is_c60_centroid and structure_path to orfs/ ({chunk[0]} .. {chunk[-1]})"
        if "pr_ref" not in st:  # the first commit opens the PR; its message is the PR title
            msg = "Add is_c60_centroid and structure_path to orfs/"
            info = api.create_commit(REPO, ops, commit_message=msg, repo_type="dataset", create_pr=True,
                                     parent_commit=BASE, commit_description=desc)
            st["pr_ref"], st["pr_url"] = info.pr_revision, info.pr_url
            print("opened", info.pr_url, flush=True)
        else:
            api.create_commit(REPO, ops, commit_message=msg, repo_type="dataset", revision=st["pr_ref"])
        st["done"] += chunk
        json.dump(st, open(state_f, "w"))
        print(f"uploaded {len(st['done'])}/{len(files)}", flush=True)
    print("PR complete:", st["pr_url"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["lookup", "rewrite", "totals", "pr"])
    ap.add_argument("--parts", default=f"0-{NPARTS - 1}")
    a = ap.parse_args()
    if a.step == "lookup":
        lookup()
    elif a.step == "rewrite":
        lo, _, hi = a.parts.partition("-")
        rewrite([int(x) for x in a.parts.split(",")] if "," in a.parts else range(int(lo), int(hi or lo) + 1))
    elif a.step == "totals":
        totals()
    else:
        pr()
