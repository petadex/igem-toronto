"""Ground-truth set, step 1: ORFs in petadex-corpus-full whose own accession is a PDB entry
(e.g. 7QJM or 5YNS_A), i.e. Alex's issue-110 "direct" class. Reads the HF orfs spine remotely
(only 7 columns, so the 28 GB sequence column is never fetched). Run from this directory.
"""
import duckdb
c = duckdb.connect()
c.sql("SET threads=8; SET memory_limit='16GB'; INSTALL httpfs; LOAD httpfs;")
B = "https://huggingface.co/datasets/petadex/petadex-corpus-full/resolve/main/orfs/part-"
files = [f"{B}{i}.parquet" for i in range(131)]
c.sql(f"""create table P as select orf_id, origin, accession, length, md5, c90_id, is_c90_centroid
  from read_parquet({files}) where origin != 'Logan'
  and regexp_matches(accession, '^[0-9][A-Za-z0-9]{{3}}(_[A-Za-z0-9]+)?$')""")
print(c.sql("select origin, count(*) n, count(distinct accession) acc, count(distinct md5) seqs, sum(is_c90_centroid::int) centroids from P group by 1"))
print(c.sql("select * from P limit 5"))
c.sql("copy P to 'data/gt/pdb_accession_orfs.tsv' (header, delimiter '\t')")
