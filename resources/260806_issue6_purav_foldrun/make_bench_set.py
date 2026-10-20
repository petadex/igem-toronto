#!/usr/bin/env python
"""Build benchmark set B: a fixed, stratified sample of the 90pid centroids for the
model / hyperparameter sweep. Reads petadex-corpus-centroids90 (HF, Dennis 2026-10-15).

Strata (seeded, reproducible):
  length bins x {NR, Logan complete (has stop), Logan truncated (no stop)}
  + AFDB-backed centroids (own AFDB structure) for an ESMFold2 vs AlphaFold2 arm
  + plasticome clusters (known plastic-degrading enzymes)

    python make_bench_set.py --parquet /project/6122149/purav/petadex-corpus/centroids90/centroids90.parquet
"""
import argparse
import os

LEN_BINS = [(29, 100), (100, 150), (150, 250), (250, 350), (350, 500), (500, 800), (800, 1281), (1281, 2049)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="data/bench")
    ap.add_argument("--per-cell", type=int, default=20, help="per length-bin x origin-class cell")
    ap.add_argument("--n-afdb", type=int, default=100)
    ap.add_argument("--n-plasticome", type=int, default=40)
    ap.add_argument("--seed", type=float, default=0.42, help="duckdb setseed value")
    args = ap.parse_args()

    import duckdb
    c = duckdb.connect()
    c.sql("SET memory_limit='24GB'; SET threads=8")
    c.sql(f"select setseed({args.seed})")
    c.sql(f"""create table C as select centroid_orf_id, centroid_origin, centroid_length, centroid_has_stop,
              centroid_afdb_s3_key, plasticome_ids, centroid_hmm_ids, n_members
              from '{args.parquet}'""")
    case_bin = "case " + " ".join(
        f"when centroid_length >= {lo} and centroid_length < {hi} then '{lo}-{hi - 1}'" for lo, hi in LEN_BINS
    ) + " end"
    cls = ("case when centroid_origin='NR' then 'NR' when centroid_origin='Logan' and centroid_has_stop "
           "then 'Logan_complete' when centroid_origin='Logan' then 'Logan_truncated' else 'PAZy' end")
    c.sql(f"create table S as select *, {case_bin} len_bin, {cls} origin_class from C")

    strat = c.sql(f"""select * exclude(r), 'strat' arm from (
        select *, row_number() over (partition by len_bin, origin_class order by random()) r
        from S where len_bin is not null and origin_class != 'PAZy') where r <= {args.per_cell}""")
    c.sql("create table B as select * from strat")
    c.sql(f"""insert into B select * exclude(r), 'afdb' arm from (
        select *, row_number() over (order by random()) r from S
        where centroid_afdb_s3_key is not null and centroid_length <= 2048
          and centroid_orf_id not in (select centroid_orf_id from B)) where r <= {args.n_afdb}""")
    c.sql(f"""insert into B select * exclude(r), 'plasticome' arm from (
        select *, row_number() over (order by random()) r from S
        where plasticome_ids is not null and centroid_length <= 2048
          and centroid_orf_id not in (select centroid_orf_id from B)) where r <= {args.n_plasticome}""")

    # sequences only for the chosen rows (the full column is ~5 GB of text)
    c.sql(f"""create table B2 as select B.*, starts_with(q.centroid_sequence,'M') m_start,
              rtrim(upper(q.centroid_sequence),'*') seq
              from B join (select centroid_orf_id, centroid_sequence from '{args.parquet}'
                           where centroid_orf_id in (select centroid_orf_id from B)) q using (centroid_orf_id)""")
    c.sql("drop table B; alter table B2 rename to B")
    os.makedirs(args.outdir, exist_ok=True)
    rows = c.sql("select * from B order by centroid_length, centroid_orf_id").fetchall()
    cols = [d[0] for d in c.sql("select * from B limit 0").description]
    with open(os.path.join(args.outdir, "bench_B.fa"), "w") as fa:
        for r in rows:
            d = dict(zip(cols, r))
            fa.write(f">orf{d['centroid_orf_id']}\n{d['seq']}\n")
    c.sql(f"copy (select * exclude(seq) from B order by centroid_length) to '{args.outdir}/bench_B_map.tsv' (header, delimiter '\t')")
    print(c.sql("select arm, origin_class, count(*) n from B group by all order by 1,2"))
    print(c.sql("select len_bin, count(*) n from B group by 1 order by min(centroid_length)"))
    print(f"wrote {len(rows)} seqs -> {args.outdir}/bench_B.fa (+ bench_B_map.tsv)")


if __name__ == "__main__":
    main()
