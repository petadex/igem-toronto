from __future__ import annotations

import hashlib
import json
import os
import platform
import random
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import psutil
import requests


ACCESSIONS = [
    "DRR001152",
    "DRR000584",
    "DRR001026",
    "DRR001104",
    "DRR000185",
    "DRR001199",
]

ENA_FIELDS = [
    "run_accession",
    "study_accession",
    "sample_accession",
    "scientific_name",
    "tax_id",
    "instrument_platform",
    "instrument_model",
    "library_source",
    "library_strategy",
    "library_layout",
    "collection_date",
    "country",
    "host",
    "environment_biome",
    "environment_feature",
    "environment_material",
]

UPSTREAM_REVISIONS = {
    "IndexThePlanet/Logan": "b3da25311802544c4c0e53c44b98380acf4d04d6",
    "rchikhi_pasteur/logan-analysis": "5e92dafb87e40834ff2123d9b55e3bb27887fedb",
    "ababaian/petadex": "a6bc326eb187925591b3988f57b4fadb1dba557d",
}


def project_root() -> Path:
    cwd = Path.cwd().resolve()
    return cwd if cwd.name == "logan-database-benchmarks" else cwd / "logan-database-benchmarks"


def ensure_directories(root: Path) -> dict[str, Path]:
    paths = {
        "raw": root / "data" / "raw",
        "metadata": root / "data" / "metadata",
        "databases": root / "data" / "databases",
        "results": root / "results",
        "figures": root / "results" / "figures",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, path: Path) -> dict:
    if not path.exists():
        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with path.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        handle.write(chunk)
    head = requests.head(url, timeout=60)
    head.raise_for_status()
    return {
        "url": url,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "etag": head.headers.get("ETag", "").strip('"'),
        "last_modified": head.headers.get("Last-Modified", ""),
    }


def acquire_logan_sample(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = ensure_directories(root)
    downloads = []
    metadata_rows = []

    for accession in ACCESSIONS:
        compressed = paths["raw"] / f"{accession}.contigs.fa.zst"
        fasta = paths["raw"] / f"{accession}.contigs.fa"
        url = (
            f"https://s3.amazonaws.com/logan-pub/c/{accession}/"
            f"{accession}.contigs.fa.zst"
        )
        item = download_file(url, compressed)
        if not fasta.exists():
            with fasta.open("wb") as output:
                subprocess.run(
                    ["zstd", "-dc", str(compressed)],
                    check=True,
                    stdout=output,
                    stderr=subprocess.PIPE,
                )
        item["decompressed_bytes"] = fasta.stat().st_size
        item["decompressed_sha256"] = sha256_file(fasta)
        downloads.append(item)

        response = requests.get(
            "https://www.ebi.ac.uk/ena/portal/api/filereport",
            params={
                "accession": accession,
                "result": "read_run",
                "fields": ",".join(ENA_FIELDS),
                "format": "json",
            },
            timeout=60,
        )
        response.raise_for_status()
        rows = response.json()
        if len(rows) != 1:
            raise RuntimeError(
                f"Expected one ENA metadata row for {accession}, received {len(rows)}"
            )
        metadata_rows.append({field: rows[0].get(field, "") for field in ENA_FIELDS})

    downloads_df = pd.DataFrame(downloads)
    metadata_df = pd.DataFrame(metadata_rows)
    downloads_df.to_csv(paths["results"] / "source_manifest.csv", index=False)
    metadata_df.to_csv(paths["metadata"] / "ena_run_metadata.csv", index=False)
    return downloads_df, metadata_df


def parse_fasta(path: Path, accession: str):
    sequence_id = None
    description = ""
    chunks = []
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    yield make_record(sequence_id, description, chunks, accession)
                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                sequence_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line.upper())
    if sequence_id is not None:
        yield make_record(sequence_id, description, chunks, accession)


def make_record(sequence_id: str, description: str, chunks: list[str], accession: str):
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"Empty sequence for {sequence_id}")
    gc = (sequence.count("G") + sequence.count("C")) / len(sequence)
    abundance = None
    for token in description.split():
        if token.startswith("ka:f:"):
            try:
                abundance = float(token.split(":", 2)[2])
            except ValueError:
                pass
    return sequence_id, sequence, len(sequence), gc, accession, abundance, description


def load_records(root: Path) -> list[tuple]:
    paths = ensure_directories(root)
    records = []
    for accession in ACCESSIONS:
        fasta = paths["raw"] / f"{accession}.contigs.fa"
        records.extend(parse_fasta(fasta, accession))
    if not records:
        raise RuntimeError("No Logan records were parsed")
    if len({record[0] for record in records}) != len(records):
        raise RuntimeError("Logan sequence identifiers are not unique")
    return records


def configure_connection(connection: sqlite3.Connection):
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("PRAGMA cache_size=-32768")


def create_schema(connection: sqlite3.Connection, variant: str):
    if variant == "no_metadata":
        connection.executescript(
            """
            CREATE TABLE sequences (
                sequence_id TEXT PRIMARY KEY,
                sequence TEXT NOT NULL,
                length INTEGER NOT NULL,
                gc_fraction REAL NOT NULL
            ) WITHOUT ROWID;
            CREATE INDEX idx_sequences_length ON sequences(length);
            CREATE INDEX idx_sequences_gc ON sequences(gc_fraction);
            """
        )
    elif variant == "metadata":
        connection.executescript(
            """
            CREATE TABLE runs (
                run_accession TEXT PRIMARY KEY,
                study_accession TEXT,
                sample_accession TEXT,
                scientific_name TEXT,
                tax_id TEXT,
                instrument_platform TEXT,
                instrument_model TEXT,
                library_source TEXT,
                library_strategy TEXT,
                library_layout TEXT,
                collection_date TEXT,
                country TEXT,
                host TEXT,
                environment_biome TEXT,
                environment_feature TEXT,
                environment_material TEXT
            ) WITHOUT ROWID;
            CREATE TABLE sequences (
                sequence_id TEXT PRIMARY KEY,
                sequence TEXT NOT NULL,
                length INTEGER NOT NULL,
                gc_fraction REAL NOT NULL,
                run_accession TEXT NOT NULL,
                abundance REAL,
                header_metadata TEXT,
                FOREIGN KEY (run_accession) REFERENCES runs(run_accession)
            ) WITHOUT ROWID;
            CREATE INDEX idx_sequences_length ON sequences(length);
            CREATE INDEX idx_sequences_gc ON sequences(gc_fraction);
            CREATE INDEX idx_sequences_run ON sequences(run_accession);
            CREATE INDEX idx_runs_scientific_name ON runs(scientific_name);
            CREATE INDEX idx_runs_library_source ON runs(library_source);
            """
        )
    else:
        raise ValueError(f"Unknown variant: {variant}")


def build_database(
    path: Path, variant: str, records: list[tuple], metadata_df: pd.DataFrame
) -> dict:
    path.unlink(missing_ok=True)
    start = time.perf_counter()
    connection = sqlite3.connect(path)
    configure_connection(connection)
    create_schema(connection, variant)
    if variant == "metadata":
        connection.executemany(
            f"INSERT INTO runs ({','.join(ENA_FIELDS)}) "
            f"VALUES ({','.join('?' for _ in ENA_FIELDS)})",
            metadata_df[ENA_FIELDS].itertuples(index=False, name=None),
        )
        connection.executemany(
            """
            INSERT INTO sequences
            (sequence_id, sequence, length, gc_fraction, run_accession,
             abundance, header_metadata)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            records,
        )
    else:
        connection.executemany(
            """
            INSERT INTO sequences
            (sequence_id, sequence, length, gc_fraction)
            VALUES (?, ?, ?, ?)
            """,
            (record[:4] for record in records),
        )
    connection.commit()
    integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    page_count = connection.execute("PRAGMA page_count").fetchone()[0]
    page_size = connection.execute("PRAGMA page_size").fetchone()[0]
    connection.close()
    elapsed = time.perf_counter() - start
    return {
        "variant": variant,
        "build_seconds": elapsed,
        "records_per_second": len(records) / elapsed,
        "database_bytes": path.stat().st_size,
        "page_count": page_count,
        "page_size": page_size,
        "integrity_check": integrity,
    }


def sequence_digest(path: Path) -> tuple[int, int, str]:
    connection = sqlite3.connect(path)
    digest = hashlib.sha256()
    count = 0
    bases = 0
    for sequence_id, sequence in connection.execute(
        "SELECT sequence_id, sequence FROM sequences ORDER BY sequence_id"
    ):
        digest.update(sequence_id.encode())
        digest.update(b"\0")
        digest.update(sequence.encode())
        digest.update(b"\n")
        count += 1
        bases += len(sequence)
    connection.close()
    return count, bases, digest.hexdigest()


def build_variants(root: Path, repetitions: int = 3):
    paths = ensure_directories(root)
    metadata_df = pd.read_csv(
        paths["metadata"] / "ena_run_metadata.csv", dtype=str, keep_default_na=False
    )
    records = load_records(root)
    build_rows = []
    variants = ["no_metadata", "metadata"]
    for repetition in range(1, repetitions + 1):
        order = variants if repetition % 2 else list(reversed(variants))
        for variant in order:
            suffix = "" if repetition == repetitions else f".rep{repetition}"
            path = paths["databases"] / f"logan_{variant}{suffix}.sqlite"
            row = build_database(path, variant, records, metadata_df)
            row["repetition"] = repetition
            build_rows.append(row)
            if repetition != repetitions:
                path.unlink()

    no_metadata_path = paths["databases"] / "logan_no_metadata.sqlite"
    metadata_path = paths["databases"] / "logan_metadata.sqlite"
    no_metadata_digest = sequence_digest(no_metadata_path)
    metadata_digest = sequence_digest(metadata_path)
    if no_metadata_digest != metadata_digest:
        raise AssertionError("Sequence payload differs between database variants")

    builds = pd.DataFrame(build_rows)
    builds.to_csv(paths["results"] / "build_measurements.csv", index=False)
    validation = {
        "sequence_count": no_metadata_digest[0],
        "total_bases": no_metadata_digest[1],
        "sequence_payload_sha256": no_metadata_digest[2],
        "variants_identical": True,
        "accessions": ACCESSIONS,
        "upstream_revisions": UPSTREAM_REVISIONS,
        "python": sys.version,
        "platform": platform.platform(),
        "sqlite": sqlite3.sqlite_version,
    }
    (paths["results"] / "validation.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    return records, builds, validation


def prepare_and_build(root: Path):
    downloads, metadata = acquire_logan_sample(root)
    records, builds, validation = build_variants(root)
    summary = pd.DataFrame(
        [
            {
                "accessions": len(ACCESSIONS),
                "sequences": len(records),
                "total_bases": sum(record[2] for record in records),
                "compressed_mib": downloads["bytes"].sum() / (1024 * 1024),
                "decompressed_mib": downloads["decompressed_bytes"].sum()
                / (1024 * 1024),
                "metadata_rows": len(metadata),
            }
        ]
    )
    return downloads, metadata, summary, builds, validation


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, int(round((len(ordered) - 1) * fraction))))
    return ordered[index]


def timed_query(
    connection: sqlite3.Connection,
    sql: str,
    params: tuple,
    repetitions: int,
    warmups: int,
):
    for _ in range(warmups):
        connection.execute(sql, params).fetchall()
    samples = []
    result = None
    for _ in range(repetitions):
        start = time.perf_counter_ns()
        result = connection.execute(sql, params).fetchall()
        samples.append((time.perf_counter_ns() - start) / 1_000_000)
    return samples, result


def benchmark_variants(root: Path, repetitions: int = 30, warmups: int = 5):
    paths = ensure_directories(root)
    db_paths = {
        "no_metadata": paths["databases"] / "logan_no_metadata.sqlite",
        "metadata": paths["databases"] / "logan_metadata.sqlite",
    }
    for path in db_paths.values():
        if not path.exists():
            raise FileNotFoundError(f"Run notebook 01 first: {path}")

    ids_connection = sqlite3.connect(db_paths["no_metadata"])
    all_ids = [
        row[0]
        for row in ids_connection.execute(
            "SELECT sequence_id FROM sequences ORDER BY sequence_id"
        )
    ]
    ids_connection.close()
    rng = random.Random(20260908)
    lookup_ids = rng.sample(all_ids, min(200, len(all_ids)))
    length_low = 500
    length_high = 2000

    common_queries = {
        "count_all": ("SELECT COUNT(*) FROM sequences", ()),
        "length_range_count": (
            "SELECT COUNT(*) FROM sequences WHERE length BETWEEN ? AND ?",
            (length_low, length_high),
        ),
        "average_length": ("SELECT AVG(length) FROM sequences", ()),
        "gc_range_count": (
            "SELECT COUNT(*) FROM sequences WHERE gc_fraction BETWEEN ? AND ?",
            (0.45, 0.55),
        ),
        "top_100_longest": (
            "SELECT sequence_id, length FROM sequences ORDER BY length DESC LIMIT 100",
            (),
        ),
    }

    raw_rows = []
    result_proof = {}
    for variant, path in db_paths.items():
        connection = sqlite3.connect(path)
        connection.execute("PRAGMA cache_size=-32768")
        connection.execute("PRAGMA temp_store=MEMORY")
        for query_name, (sql, params) in common_queries.items():
            samples, result = timed_query(
                connection, sql, params, repetitions, warmups
            )
            result_proof.setdefault(query_name, {})[variant] = result
            raw_rows.extend(
                {
                    "variant": variant,
                    "query": query_name,
                    "repetition": index + 1,
                    "latency_ms": sample,
                }
                for index, sample in enumerate(samples)
            )

        for _ in range(warmups):
            for sequence_id in lookup_ids:
                connection.execute(
                    "SELECT sequence, length FROM sequences WHERE sequence_id = ?",
                    (sequence_id,),
                ).fetchone()
        for repetition in range(repetitions):
            start = time.perf_counter_ns()
            rows = [
                connection.execute(
                    "SELECT sequence, length FROM sequences WHERE sequence_id = ?",
                    (sequence_id,),
                ).fetchone()
                for sequence_id in lookup_ids
            ]
            elapsed = (time.perf_counter_ns() - start) / 1_000_000
            raw_rows.append(
                {
                    "variant": variant,
                    "query": "indexed_lookup_batch_200",
                    "repetition": repetition + 1,
                    "latency_ms": elapsed,
                }
            )
            result_proof.setdefault("indexed_lookup_batch_200", {})[variant] = rows

        accession = ACCESSIONS[0]
        if variant == "metadata":
            accession_sql = (
                "SELECT COUNT(*) FROM sequences WHERE run_accession = ?",
                (accession,),
            )
        else:
            accession_sql = (
                "SELECT COUNT(*) FROM sequences WHERE sequence_id GLOB ?",
                (f"{accession}_*",),
            )
        samples, result = timed_query(
            connection, accession_sql[0], accession_sql[1], repetitions, warmups
        )
        result_proof.setdefault("accession_filter", {})[variant] = result
        raw_rows.extend(
            {
                "variant": variant,
                "query": "accession_filter",
                "repetition": index + 1,
                "latency_ms": sample,
            }
            for index, sample in enumerate(samples)
        )
        connection.close()

    for query_name, variants in result_proof.items():
        if variants["no_metadata"] != variants["metadata"]:
            raise AssertionError(f"Result mismatch for {query_name}")

    persisted_proof = {}
    for query_name, variants in result_proof.items():
        digests = {}
        for variant, result in variants.items():
            encoded = json.dumps(result, separators=(",", ":"), ensure_ascii=True).encode()
            digests[variant] = hashlib.sha256(encoded).hexdigest()
        persisted_proof[query_name] = {
            "equivalent": digests["no_metadata"] == digests["metadata"],
            "no_metadata_sha256": digests["no_metadata"],
            "metadata_sha256": digests["metadata"],
            "result_rows": len(variants["metadata"]),
        }
    (paths["results"] / "shared_query_proofs.json").write_text(
        json.dumps(persisted_proof, indent=2), encoding="utf-8"
    )

    raw = pd.DataFrame(raw_rows)
    summary = (
        raw.groupby(["variant", "query"])["latency_ms"]
        .agg(
            median_ms="median",
            mean_ms="mean",
            stddev_ms="std",
            minimum_ms="min",
            maximum_ms="max",
        )
        .reset_index()
    )
    p95 = (
        raw.groupby(["variant", "query"])["latency_ms"]
        .apply(lambda values: percentile(values.tolist(), 0.95))
        .rename("p95_ms")
        .reset_index()
    )
    summary = summary.merge(p95, on=["variant", "query"])
    summary["workloads_per_second"] = 1000 / summary["median_ms"]
    summary["logical_operations_per_second"] = summary.apply(
        lambda row: row["workloads_per_second"] * 200
        if row["query"] == "indexed_lookup_batch_200"
        else row["workloads_per_second"],
        axis=1,
    )

    raw.to_csv(paths["results"] / "query_measurements.csv", index=False)
    summary.to_csv(paths["results"] / "query_summary.csv", index=False)
    return raw, summary, result_proof


def storage_and_plans(root: Path):
    paths = ensure_directories(root)
    rows = []
    plans = {}
    for variant in ["no_metadata", "metadata"]:
        path = paths["databases"] / f"logan_{variant}.sqlite"
        connection = sqlite3.connect(path)
        page_count = connection.execute("PRAGMA page_count").fetchone()[0]
        page_size = connection.execute("PRAGMA page_size").fetchone()[0]
        integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
        sequence_count = connection.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
        first_id = connection.execute(
            "SELECT sequence_id FROM sequences ORDER BY sequence_id LIMIT 1"
        ).fetchone()[0]
        if variant == "metadata":
            run_count = connection.execute("SELECT COUNT(*) FROM runs").fetchone()[0]
            accession_plan = (
                "SELECT COUNT(*) FROM sequences WHERE run_accession = ?",
                (ACCESSIONS[0],),
            )
        else:
            run_count = 0
            accession_plan = (
                "SELECT COUNT(*) FROM sequences WHERE sequence_id GLOB ?",
                (f"{ACCESSIONS[0]}_*",),
            )
        plan_queries = {
            "count_all": ("SELECT COUNT(*) FROM sequences", ()),
            "length_range_count": (
                "SELECT COUNT(*) FROM sequences WHERE length BETWEEN ? AND ?",
                (500, 2000),
            ),
            "average_length": ("SELECT AVG(length) FROM sequences", ()),
            "gc_range_count": (
                "SELECT COUNT(*) FROM sequences WHERE gc_fraction BETWEEN ? AND ?",
                (0.45, 0.55),
            ),
            "top_100_longest": (
                "SELECT sequence_id, length FROM sequences ORDER BY length DESC LIMIT 100",
                (),
            ),
            "indexed_lookup_batch_200": (
                "SELECT sequence, length FROM sequences WHERE sequence_id = ?",
                (first_id,),
            ),
            "accession_filter": accession_plan,
        }
        plans[variant] = {
            name: connection.execute(f"EXPLAIN QUERY PLAN {sql}", params).fetchall()
            for name, (sql, params) in plan_queries.items()
        }
        rows.append(
            {
                "variant": variant,
                "database_mib": path.stat().st_size / (1024 * 1024),
                "page_count": page_count,
                "page_size": page_size,
                "sequence_count": sequence_count,
                "metadata_rows": run_count,
                "integrity_check": integrity,
            }
        )
        connection.close()
    frame = pd.DataFrame(rows)
    frame.to_csv(paths["results"] / "storage_summary.csv", index=False)
    (paths["results"] / "query_plans.json").write_text(
        json.dumps(plans, indent=2), encoding="utf-8"
    )
    return frame, plans


def metadata_capabilities(root: Path, repetitions: int = 30, warmups: int = 5):
    paths = ensure_directories(root)
    connection = sqlite3.connect(paths["databases"] / "logan_metadata.sqlite")
    checks = {
        "sequences_by_library_source": """
            SELECT r.library_source, COUNT(*)
            FROM sequences s JOIN runs r USING (run_accession)
            GROUP BY r.library_source ORDER BY COUNT(*) DESC
        """,
        "sequences_by_run_scientific_name": """
            SELECT r.scientific_name, COUNT(*)
            FROM sequences s JOIN runs r USING (run_accession)
            GROUP BY r.scientific_name ORDER BY COUNT(*) DESC
        """,
        "high_abundance_sequences_by_run_scientific_name": """
            SELECT r.scientific_name, COUNT(*)
            FROM sequences s JOIN runs r USING (run_accession)
            WHERE s.abundance >= 10
            GROUP BY r.scientific_name ORDER BY COUNT(*) DESC
        """,
    }
    rows = []
    for name, sql in checks.items():
        samples, result = timed_query(connection, sql, (), repetitions, warmups)
        rows.append(
            {
                "capability": name,
                "median_ms": statistics_median(samples),
                "p95_ms": percentile(samples, 0.95),
                "minimum_ms": min(samples),
                "maximum_ms": max(samples),
                "repetitions": repetitions,
                "result_rows": len(result),
                "result_sha256": hashlib.sha256(
                    json.dumps(result, separators=(",", ":")).encode()
                ).hexdigest(),
                "result": json.dumps(result),
            }
        )
    connection.close()
    frame = pd.DataFrame(rows)
    frame.to_csv(paths["results"] / "metadata_capabilities.csv", index=False)
    return frame


def memory_benchmark(root: Path, repetitions: int = 5):
    paths = ensure_directories(root)
    rows = []
    worker = r'''
import json, os, sqlite3, sys, threading, time
import psutil
path, variant = sys.argv[1], sys.argv[2]
process = psutil.Process(os.getpid())
baseline = process.memory_info().rss
peak = baseline
stop = threading.Event()
def sample():
    global peak
    while not stop.is_set():
        peak = max(peak, process.memory_info().rss)
        time.sleep(0.002)
thread = threading.Thread(target=sample, daemon=True)
thread.start()
start = time.perf_counter()
connection = sqlite3.connect(path)
connection.execute("PRAGMA cache_size=-32768")
for _ in range(20):
    connection.execute("SELECT AVG(length), AVG(gc_fraction) FROM sequences").fetchone()
    connection.execute("SELECT sequence_id, sequence FROM sequences ORDER BY length DESC LIMIT 250").fetchall()
    if variant == "metadata":
        connection.execute("SELECT r.library_source, COUNT(*) FROM sequences s JOIN runs r USING (run_accession) GROUP BY r.library_source").fetchall()
connection.close()
stop.set()
thread.join()
peak = max(peak, process.memory_info().rss)
print(json.dumps({
    "baseline_rss_mib": baseline / (1024 * 1024),
    "peak_rss_mib": peak / (1024 * 1024),
    "rss_delta_mib": (peak - baseline) / (1024 * 1024),
    "workload_seconds": time.perf_counter() - start
}))
'''
    for repetition in range(1, repetitions + 1):
        order = ["no_metadata", "metadata"]
        if repetition % 2 == 0:
            order.reverse()
        for variant in order:
            path = paths["databases"] / f"logan_{variant}.sqlite"
            completed = subprocess.run(
                [sys.executable, "-c", worker, str(path), variant],
                check=True,
                text=True,
                capture_output=True,
            )
            row = json.loads(completed.stdout)
            row.update({"variant": variant, "repetition": repetition})
            rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(paths["results"] / "memory_summary.csv", index=False)
    return frame


def create_figures(root: Path, query_summary: pd.DataFrame, storage: pd.DataFrame):
    paths = ensure_directories(root)
    pivot = query_summary.pivot(index="query", columns="variant", values="median_ms")
    ax = pivot.plot(kind="bar", figsize=(12, 6), logy=True, color=["#E85D3F", "#3C82F6"])
    ax.set_ylabel("Median latency, milliseconds, log scale")
    ax.set_xlabel("Query workload")
    ax.set_title("Logan SQLite benchmark: metadata impact on query latency")
    ax.grid(axis="y", alpha=0.25)
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    latency_path = paths["figures"] / "query_latency.png"
    plt.savefig(latency_path, dpi=180)
    plt.close()

    ax = storage.set_index("variant")["database_mib"].plot(
        kind="bar", figsize=(7, 5), color=["#E85D3F", "#3C82F6"]
    )
    ax.set_ylabel("Database size, MiB")
    ax.set_xlabel("Variant")
    ax.set_title("Storage cost of the metadata-enriched schema")
    ax.grid(axis="y", alpha=0.25)
    plt.xticks(rotation=0)
    plt.tight_layout()
    storage_path = paths["figures"] / "storage_size.png"
    plt.savefig(storage_path, dpi=180)
    plt.close()
    return latency_path, storage_path


def benchmark_environment() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "timestamp_utc": pd.Timestamp.now(tz="UTC").isoformat(),
                "python": sys.version.split()[0],
                "sqlite": sqlite3.sqlite_version,
                "platform": platform.platform(),
                "cpu_logical": psutil.cpu_count(logical=True),
                "cpu_physical": psutil.cpu_count(logical=False),
                "memory_gib": psutil.virtual_memory().total / (1024**3),
                "logan_revision": UPSTREAM_REVISIONS["IndexThePlanet/Logan"],
                "logan_analysis_revision": UPSTREAM_REVISIONS[
                    "rchikhi_pasteur/logan-analysis"
                ],
                "petadex_revision": UPSTREAM_REVISIONS["ababaian/petadex"],
            }
        ]
    )


def statistics_median(values):
    return float(pd.Series(values).median())