#!/usr/bin/env python
"""Build the hybrid fold input for a centroid level straight from the PETadex HF corpora.

Centroid IDs come from petadex/petadex-corpus-full dims/clusters<NN>.parquet (centroid_orf_id);
sequences from petadex/petadex-corpus-centroids90 (centroid_sequence). Every 60% and 30%
centroid is also a 90% centroid, so centroids90 holds all sequences.

Tiers (Part II decision, notebook Step 15):
  short_full.fa   L <  --split (150)            -> ESMFold2 Full,  32 steps / 3 loops
  mid_fast.fa     --split <= L <= --l40s-cap    -> ESMFold2-Fast,  32 steps / 3 loops
  tail_h100.fa    L >  --l40s-cap (1740)        -> H100 (separate run)
Skipped (recorded in manifest.tsv, not written to any FASTA):
  afdb       centroid has its own AFDB structure (unless --include-afdb)
  duplicate  exact sequence (md5) already present under a lower orf_id

    python make_fold_input_hf.py --level 60 --outdir data/fold60pid_input
"""
import argparse
import os

CORPUS = "/project/6122149/purav/petadex-corpus"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, choices=[30, 60, 90], default=60)
    ap.add_argument("--corpus", default=CORPUS)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--split", type=int, default=150, help="Full below this length, Fast at or above")
    ap.add_argument("--l40s-cap", type=int, default=1740, help="longer sequences go to the H100 tier")
    ap.add_argument("--include-afdb", action="store_true",
                    help="also fold centroids that already have their own AFDB structure")
    args = ap.parse_args()

    import duckdb
    c = duckdb.connect()
    c.sql("SET threads=8; SET memory_limit='24GB'")
    c90 = f"{args.corpus}/centroids90/centroids90.parquet"
    if args.level == 90:
        ids = f"select centroid_orf_id from '{c90}'"
    else:
        ids = f"select centroid_orf_id from '{args.corpus}/full/dims/clusters{args.level}.parquet'"

    c.sql(f"""create table T as
        select q.centroid_orf_id orf_id, q.c90_id, q.c60_id, q.c30_id, q.centroid_origin origin,
               q.centroid_md5 md5, q.centroid_afdb_s3_key afdb_s3_key,
               rtrim(upper(q.centroid_sequence), '*') seq
        from '{c90}' q where q.centroid_orf_id in ({ids})""")
    n_ids = c.sql(f"select count(*) from ({ids})").fetchone()[0]
    n_found = c.sql("select count(*) from T").fetchone()[0]
    if n_found != n_ids:
        raise SystemExit(f"{n_ids - n_found} centroid ids have no sequence in centroids90")

    c.sql(f"""create table M as select *, length(seq) L,
        case when afdb_s3_key is not null and not {str(args.include_afdb).lower()} then 'afdb'
             when row_number() over (partition by md5 order by orf_id) > 1 then 'duplicate'
             else null end skip_reason,
        case when length(seq) < {args.split} then 'short_full'
             when length(seq) <= {args.l40s_cap} then 'mid_fast'
             else 'tail_h100' end tier
        from T""")

    os.makedirs(args.outdir, exist_ok=True)
    for tier in ("short_full", "mid_fast", "tail_h100"):
        rows = c.sql(f"select orf_id, seq from M where tier='{tier}' and skip_reason is null order by L, orf_id").fetchall()
        with open(os.path.join(args.outdir, f"{tier}.fa"), "w") as fa:
            for oid, seq in rows:
                fa.write(f">orf{oid}\n{seq}\n")
    # duplicates map to the copy that is folded, so their structure can be symlinked/joined later
    c.sql(f"""copy (select m.orf_id, m.c90_id, m.c60_id, m.c30_id, m.origin, m.L length, m.tier, m.skip_reason,
                   m.md5, m.afdb_s3_key,
                   case when m.skip_reason='duplicate' then first_value(m.orf_id) over (partition by m.md5 order by m.orf_id) end folded_as
              from M m order by m.L, m.orf_id) to '{args.outdir}/manifest.tsv' (header, delimiter '\t')""")

    print(f"level {args.level}: {n_ids:,} centroids")
    print(c.sql("""select tier, count(*) filter (where skip_reason is null) to_fold,
                   count(*) filter (where skip_reason='afdb') skip_afdb, count(*) filter (where skip_reason='duplicate') skip_dup,
                   sum(L) filter (where skip_reason is null) residues, min(L) min_len, max(L) max_len
                   from M group by tier order by min(L)"""))
    print(f"wrote {args.outdir}/{{short_full,mid_fast,tail_h100}}.fa and manifest.tsv")


if __name__ == "__main__":
    main()
