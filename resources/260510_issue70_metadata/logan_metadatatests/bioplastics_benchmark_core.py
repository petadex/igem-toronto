from __future__ import annotations

import hashlib
import json
import random
import shutil
import sqlite3
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from benchmark_core import (
    ACCESSIONS as BACKGROUND_ACCESSIONS,
    configure_connection,
    download_file,
    ensure_directories,
    parse_fasta,
    percentile,
    sha256_file,
    timed_query,
)


PLASTIC_CONTEXT_ACCESSIONS = {
    "DRR821254": {
        "study_accession": "PRJDB39034",
        "study_title": (
            "Multi-omics reveals the response of anaerobic microorganisms "
            "to microplastic pressure with different aging states"
        ),
        "context_class": "plastic_context",
        "selection_term": "microplastic",
    },
    "DRR821257": {
        "study_accession": "PRJDB39034",
        "study_title": (
            "Multi-omics reveals the response of anaerobic microorganisms "
            "to microplastic pressure with different aging states"
        ),
        "context_class": "plastic_context",
        "selection_term": "microplastic",
    },
    "DRR821260": {
        "study_accession": "PRJDB39034",
        "study_title": (
            "Multi-omics reveals the response of anaerobic microorganisms "
            "to microplastic pressure with different aging states"
        ),
        "context_class": "plastic_context",
        "selection_term": "microplastic",
    },
}
ALL_ACCESSIONS = BACKGROUND_ACCESSIONS + list(PLASTIC_CONTEXT_ACCESSIONS)

ROOT = Path(__file__).resolve().parent
WORKSPACE = ROOT.parent
PLASTICDB_RAW = (
    WORKSPACE
    / "plastic-biodegradation-analysis"
    / "data"
    / "plasticdb_microorganisms.tsv"
)
PLASTICDB_SCORED = (
    WORKSPACE
    / "plastic-biodegradation-analysis"
    / "outputs"
    / "reports"
    / "plasticdb_scored.csv"
)
PAZY_CSV = (
    WORKSPACE
    / "plastic-biodegradation-analysis"
    / "outputs"
    / "reports"
    / "pazy_proteins.csv"
)

BIOPLASTIC_OR_BIODEGRADABLE_PREFIXES = (
    "PHA", "PHB", "PHO", "PHV", "P3H", "P4H",
    "PLA", "PCL", "PBS", "PBAT", "PES", "PEC", "PVA", "OPVA",
    "PEF", "PEA", "PBSET", "PTS", "PTC",
)
COMMERCIAL_BIODEGRADABLE_FORMULATIONS = {"IMPRANIL", "ECOVIOFT"}
CLASSIFICATION_EXPECTATIONS = {
    "PHO": "bioplastic_or_biodegradable_polymer",
    "PHV": "bioplastic_or_biodegradable_polymer",
    "P3HP": "bioplastic_or_biodegradable_polymer",
    "P3HV": "bioplastic_or_biodegradable_polymer",
    "P(3HB-co-3HV)": "bioplastic_or_biodegradable_polymer",
    "P(3HB-co-3MP)": "bioplastic_or_biodegradable_polymer",
    "P(3HB-co-4HB)": "bioplastic_or_biodegradable_polymer",
    "PEF": "bioplastic_or_biodegradable_polymer",
    "PEA": "bioplastic_or_biodegradable_polymer",
    "PBSeT": "bioplastic_or_biodegradable_polymer",
    "PTS": "bioplastic_or_biodegradable_polymer",
    "PTC": "bioplastic_or_biodegradable_polymer",
    "Impranil": "commercial_biodegradable_formulation",
    "Ecovio-FT": "commercial_biodegradable_formulation",
    "NR": "natural_biopolymer",
}
ENZYME_FAMILIES = {
    "PETase": ["petase", "pet hydrolase", "pet-hydrolas"],
    "Cutinase": ["cutin", "cutinase", "lcc"],
    "Lipase/esterase": ["lipase", "esterase", "carboxylesterase"],
    "Laccase": ["laccase", "multicopper oxidase"],
    "Peroxidase": ["peroxidase"],
    "Depolymerase": ["depolymerase"],
    "Protease": ["protease", "subtilisin"],
    "Alkane hydroxylase": ["alkane hydroxylase", "monooxygenase", "alkb"],
    "Oxidase": ["oxidase", "catalase"],
}
AA = set("ACDEFGHIKLMNPQRSTVWY")
CODONS = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}
HOMOLOGY_THRESHOLDS = {
    "evalue_max": 1e-5,
    "identity_min_pct": 30.0,
    "aligned_aa_min": 50,
    "subject_coverage_min": 0.50,
}
RULE_SET_VERSION = "bioplastics-metadata-rules-1"
SOURCE_URLS = {
    "PlasticDB": "https://plasticdb.org/static/degraders_list.tsv",
    "PAZy": "https://www.pazy.eu/proteins",
    "UniProt": "https://rest.uniprot.org/uniprotkb/{accession}.fasta",
    "Logan": "https://s3.amazonaws.com/logan-pub/c/{accession}/{accession}.contigs.fa.zst",
    "ENA": "https://www.ebi.ac.uk/ena/portal/api/filereport",
}


def classify_plastic(value: str) -> str:
    normalized = "".join(
        character for character in str(value).strip().upper()
        if character.isalnum()
    )
    if normalized in COMMERCIAL_BIODEGRADABLE_FORMULATIONS:
        return "commercial_biodegradable_formulation"
    if normalized == "NR":
        return "natural_biopolymer"
    if normalized.startswith(BIOPLASTIC_OR_BIODEGRADABLE_PREFIXES):
        return "bioplastic_or_biodegradable_polymer"
    return "other_plastic"


def classify_enzyme(value: str) -> str:
    lowered = str(value).lower()
    for family, terms in ENZYME_FAMILIES.items():
        if any(term in lowered for term in terms):
            return family
    return "Other enzyme"


def clean_protein(value: str) -> str:
    return "".join(character for character in str(value).upper() if character in AA)


def protein_properties(sequence: str) -> dict:
    analysis = ProteinAnalysis(sequence)
    return {
        "protein_length": len(sequence),
        "molecular_weight": analysis.molecular_weight(),
        "isoelectric_point": analysis.isoelectric_point(),
        "instability_index": analysis.instability_index(),
        "gravy": analysis.gravy(),
        "aromaticity": analysis.aromaticity(),
    }


def fetch_uniprot(accession: str, cache_dir: Path) -> tuple[str, dict]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{accession}.fasta"
    url = SOURCE_URLS["UniProt"].format(accession=accession)
    retrieved_at = ""
    if not path.exists():
        response = requests.get(url, timeout=60)
        if response.status_code != 200:
            return "", {
                "uniprot": accession,
                "status": "not_resolved",
                "sequence_length": 0,
                "url": url,
                "cached_fasta": "",
                "sha256": "",
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            }
        path.write_text(response.text, encoding="utf-8")
        retrieved_at = datetime.now(timezone.utc).isoformat()
    else:
        retrieved_at = datetime.fromtimestamp(
            path.stat().st_mtime, tz=timezone.utc
        ).isoformat()
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if not line.startswith(">")]
    sequence = clean_protein("".join(lines))
    return sequence, {
        "uniprot": accession,
        "status": "resolved" if sequence else "not_resolved",
        "sequence_length": len(sequence),
        "url": url,
        "cached_fasta": str(path),
        "sha256": sha256_file(path),
        "retrieved_at_utc": retrieved_at,
    }


def acquire_corrected_sample(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = ensure_directories(root)
    downloads = []
    study_rows = []
    ena_fields = [
        "run_accession", "study_accession", "study_title", "scientific_name",
        "tax_id", "library_source", "library_strategy", "library_layout",
    ]
    for accession in ALL_ACCESSIONS:
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
        item["context_class"] = (
            "plastic_context"
            if accession in PLASTIC_CONTEXT_ACCESSIONS
            else "background_control"
        )
        downloads.append(item)

        response = requests.get(
            "https://www.ebi.ac.uk/ena/portal/api/filereport",
            params={
                "accession": accession,
                "result": "read_run",
                "fields": ",".join(ena_fields),
                "format": "json",
            },
            timeout=60,
        )
        response.raise_for_status()
        rows = response.json()
        if len(rows) != 1:
            raise RuntimeError(f"Expected one ENA row for {accession}")
        row = {field: rows[0].get(field, "") for field in ena_fields}
        row["context_class"] = item["context_class"]
        row["selection_term"] = PLASTIC_CONTEXT_ACCESSIONS.get(
            accession, {}
        ).get("selection_term", "")
        study_rows.append(row)

    downloads_df = pd.DataFrame(downloads)
    studies_df = pd.DataFrame(study_rows)
    downloads_df.to_csv(
        paths["results"] / "bioplastics_source_manifest.csv", index=False
    )
    studies_df.to_csv(paths["metadata"] / "logan_study_context.csv", index=False)
    return downloads_df, studies_df


def load_corrected_records(root: Path) -> list[tuple]:
    paths = ensure_directories(root)
    records = []
    for accession in ALL_ACCESSIONS:
        for record in parse_fasta(paths["raw"] / f"{accession}.contigs.fa", accession):
            records.append(record[:5])
    if len({record[0] for record in records}) != len(records):
        raise RuntimeError("Sequence identifiers are not unique")
    return records


def build_reference_library(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = ensure_directories(root)
    raw = pd.read_csv(PLASTICDB_RAW, sep="\t", dtype=str, keep_default_na=False)
    scored = pd.read_csv(PLASTICDB_SCORED, dtype=str, keep_default_na=False)
    if len(raw) != len(scored):
        raise RuntimeError("PlasticDB raw and scored row counts differ")

    plastic_audit = (
        raw.groupby("Plastic", dropna=False)
        .size()
        .rename("source_rows")
        .reset_index()
        .rename(columns={"Plastic": "source_label"})
    )
    plastic_audit["normalized_label"] = plastic_audit["source_label"].map(
        lambda value: "".join(
            character for character in str(value).strip().upper()
            if character.isalnum()
        )
    )
    plastic_audit["assigned_class"] = plastic_audit["source_label"].map(
        classify_plastic
    )
    plastic_audit["classification_rule"] = plastic_audit.apply(
        lambda row: (
            "explicit commercial formulation"
            if row["assigned_class"] == "commercial_biodegradable_formulation"
            else "explicit natural-biopolymer separation"
            if row["assigned_class"] == "natural_biopolymer"
            else "normalized polymer-family prefix"
            if row["assigned_class"] == "bioplastic_or_biodegradable_polymer"
            else "outside declared target families"
        ),
        axis=1,
    )
    observed_classes = dict(
        zip(plastic_audit["source_label"], plastic_audit["assigned_class"])
    )
    for label, expected_class in CLASSIFICATION_EXPECTATIONS.items():
        if label in observed_classes and observed_classes[label] != expected_class:
            raise AssertionError(
                f"Plastic classification mismatch for {label}: "
                f"{observed_classes[label]} != {expected_class}"
            )
    plastic_audit.to_csv(
        paths["results"] / "plastic_classification_audit.csv", index=False
    )

    annotations = []
    sequences = {}
    for index, row in raw.iterrows():
        sequence = clean_protein(row["Sequence"])
        if len(sequence) < 30:
            continue
        reference_id = "ref_" + hashlib.sha256(sequence.encode()).hexdigest()[:16]
        sequences[reference_id] = sequence
        scored_row = scored.iloc[index]
        annotations.append(
            {
                "reference_id": reference_id,
                "source": "PlasticDB",
                "organism": row["Microorganism"],
                "tax_id": row["Tax ID"],
                "plastic": row["Plastic"],
                "plastic_class": classify_plastic(row["Plastic"]),
                "enzyme_name": row["Enzyme"],
                "enzyme_family": classify_enzyme(row["Enzyme"]),
                "genbank_id": row["GenbankID"],
                "evidence": row["Evidence"],
                "evidence_score": scored_row.get("evidence_score", ""),
                "evidence_tier": scored_row.get("evidence_tier", ""),
                "year": row["Year"],
                "doi": row["DOI"],
            }
        )

    pazy = pd.read_csv(PAZY_CSV, dtype=str, keep_default_na=False)
    pazy_resolution = []
    for _, row in pazy.iterrows():
        accession = row.get("uniprot", "").strip()
        if not accession:
            pazy_resolution.append(
                {
                    "uniprot": "", "status": "missing_accession",
                    "sequence_length": 0, "url": "", "cached_fasta": "",
                    "sha256": "", "retrieved_at_utc": "",
                }
            )
            continue
        sequence, resolution = fetch_uniprot(
            accession, paths["metadata"] / "uniprot"
        )
        pazy_resolution.append(resolution)
        if len(sequence) < 30:
            continue
        reference_id = "ref_" + hashlib.sha256(sequence.encode()).hexdigest()[:16]
        sequences[reference_id] = sequence
        annotations.append(
            {
                "reference_id": reference_id,
                "source": "PAZy",
                "organism": row.get("organism", ""),
                "tax_id": "",
                "plastic": row.get("plastic", ""),
                "plastic_class": classify_plastic(row.get("plastic", "")),
                "enzyme_name": row.get("enzyme_name", ""),
                "enzyme_family": classify_enzyme(row.get("enzyme_name", "")),
                "genbank_id": accession,
                "evidence": "thoroughly characterised PAZy protein",
                "evidence_score": "",
                "evidence_tier": "PAZy",
                "year": row.get("year", ""),
                "doi": "",
            }
        )

    reference_rows = []
    for reference_id, sequence in sorted(sequences.items()):
        row = {
            "reference_id": reference_id,
            "sequence": sequence,
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
        row.update(protein_properties(sequence))
        reference_rows.append(row)
    references = pd.DataFrame(reference_rows)
    annotations_df = pd.DataFrame(annotations).drop_duplicates()
    references.to_csv(paths["metadata"] / "biodegradation_references.csv", index=False)
    annotations_df.to_csv(
        paths["metadata"] / "biodegradation_reference_annotations.csv", index=False
    )
    pd.DataFrame(pazy_resolution).to_csv(
        paths["results"] / "pazy_uniprot_resolution.csv", index=False
    )

    fasta = paths["metadata"] / "biodegradation_references.faa"
    with fasta.open("w", encoding="utf-8") as handle:
        for row in references.itertuples(index=False):
            handle.write(f">{row.reference_id}\n{row.sequence}\n")
    manifest = {
        "manifest_created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rule_set_version": RULE_SET_VERSION,
        "source_urls": SOURCE_URLS,
        "plasticdb_raw_path": str(PLASTICDB_RAW),
        "plasticdb_raw_sha256": sha256_file(PLASTICDB_RAW),
        "plasticdb_snapshot_mtime_utc": datetime.fromtimestamp(
            PLASTICDB_RAW.stat().st_mtime, tz=timezone.utc
        ).isoformat(),
        "plasticdb_scored_path": str(PLASTICDB_SCORED),
        "plasticdb_scored_sha256": sha256_file(PLASTICDB_SCORED),
        "plasticdb_scored_snapshot_mtime_utc": datetime.fromtimestamp(
            PLASTICDB_SCORED.stat().st_mtime, tz=timezone.utc
        ).isoformat(),
        "pazy_path": str(PAZY_CSV),
        "pazy_sha256": sha256_file(PAZY_CSV),
        "pazy_snapshot_mtime_utc": datetime.fromtimestamp(
            PAZY_CSV.stat().st_mtime, tz=timezone.utc
        ).isoformat(),
        "unique_reference_sequences": len(references),
        "unique_reference_sequence_hashes": references["sequence_sha256"].nunique(),
        "annotation_rows": len(annotations_df),
        "reference_fasta_sha256": sha256_file(fasta),
        "thresholds": HOMOLOGY_THRESHOLDS,
        "threshold_interpretation": (
            "Heuristic candidate-similarity retrieval filter only; not calibrated "
            "to infer plastic-degradation function."
        ),
        "threshold_rationale": {
            "evalue_max": "Limits chance alignments in this reference search.",
            "identity_min_pct": "Permissive remote-similarity screen.",
            "aligned_aa_min": "Rejects very short local motif matches.",
            "subject_coverage_min": "Requires at least half of the reference protein span.",
        },
        "plastic_classification_scope": (
            "Broad study category combining bio-based and/or biodegradable "
            "polymers; it does not assert that every material is both."
        ),
        "bioplastic_or_biodegradable_normalized_prefixes": (
            BIOPLASTIC_OR_BIODEGRADABLE_PREFIXES
        ),
        "commercial_biodegradable_formulations": sorted(
            COMMERCIAL_BIODEGRADABLE_FORMULATIONS
        ),
        "classification_expectations": CLASSIFICATION_EXPECTATIONS,
        "classification_audit_path": str(
            paths["results"] / "plastic_classification_audit.csv"
        ),
        "classification_audit_sha256": sha256_file(
            paths["results"] / "plastic_classification_audit.csv"
        ),
        "classification_normalization": (
            "Uppercase and remove non-alphanumeric characters before prefix matching."
        ),
        "deliberate_separation": (
            "Impranil and Ecovio-FT are reported as commercial biodegradable "
            "formulations rather than polymer identities. Natural rubber (NR) "
            "is reported as a natural biopolymer rather than a plastic."
        ),
        "enzyme_family_keyword_rules": ENZYME_FAMILIES,
        "petadex_derivations": [
            "bioplastic classification",
            "enzyme-family classification",
            "ProtParam biochemical properties",
        ],
    }
    (paths["results"] / "bioplastics_reference_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return references, annotations_df


def _write_combined_queries(root: Path, records: list[tuple]) -> Path:
    paths = ensure_directories(root)
    path = paths["metadata"] / "logan_queries.fna"
    with path.open("w", encoding="utf-8") as handle:
        for sequence_id, sequence, *_ in records:
            handle.write(f">{sequence_id}\n{sequence}\n")
    return path


def _parse_diamond(path: Path) -> pd.DataFrame:
    columns = [
        "sequence_id", "reference_id", "identity_pct", "aligned_aa",
        "query_length_nt", "subject_length_aa", "subject_start", "subject_end",
        "evalue", "bitscore",
    ]
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns + ["subject_coverage"])
    frame = pd.read_csv(path, sep="\t", names=columns)
    frame["subject_coverage"] = (
        (frame["subject_end"] - frame["subject_start"]).abs() + 1
    ) / frame["subject_length_aa"]
    return frame


def _qualified(frame: pd.DataFrame) -> pd.DataFrame:
    t = HOMOLOGY_THRESHOLDS
    return frame[
        (frame["evalue"] <= t["evalue_max"])
        & (frame["identity_pct"] >= t["identity_min_pct"])
        & (frame["aligned_aa"] >= t["aligned_aa_min"])
        & (frame["subject_coverage"] >= t["subject_coverage_min"])
    ].copy()


def run_homology(root: Path, records: list[tuple], references: pd.DataFrame):
    paths = ensure_directories(root)
    reference_fasta = paths["metadata"] / "biodegradation_references.faa"
    database = paths["metadata"] / "biodegradation_references"
    query_fasta = _write_combined_queries(root, records)
    raw_output = paths["results"] / "logan_biodegradation_homology_raw.tsv"

    start = time.perf_counter()
    subprocess.run(
        ["diamond", "makedb", "--in", str(reference_fasta), "--db", str(database)],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [
            "diamond", "blastx", "--db", str(database), "--query", str(query_fasta),
            "--out", str(raw_output), "--outfmt", "6", "qseqid", "sseqid",
            "pident", "length", "qlen", "slen", "sstart", "send", "evalue", "bitscore",
            "--max-target-seqs", "10", "--evalue", "1e-5", "--threads", "2",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - start
    raw = _parse_diamond(raw_output)
    qualified = _qualified(raw)
    if not qualified.empty:
        qualified = qualified.sort_values(
            ["sequence_id", "bitscore"], ascending=[True, False]
        )
        qualified["hit_rank"] = qualified.groupby("sequence_id").cumcount() + 1
    else:
        qualified["hit_rank"] = pd.Series(dtype=int)
    qualified.to_csv(
        paths["results"] / "logan_biodegradation_homology_candidates.csv", index=False
    )

    eligible_controls = references[references["protein_length"] >= 50]
    controls = eligible_controls.sort_values("reference_id").head(
        min(25, len(eligible_controls))
    )
    rng = random.Random(20260908)
    control_fasta = paths["metadata"] / "homology_controls.fna"
    expected = {}
    reference_sequences = dict(
        zip(references["reference_id"], references["sequence"])
    )
    with control_fasta.open("w", encoding="utf-8") as handle:
        for row in controls.itertuples(index=False):
            nucleotide = "".join(CODONS[aa] for aa in row.sequence)
            positive_id = f"positive__{row.reference_id}"
            handle.write(f">{positive_id}\n{nucleotide}\n")
            expected[positive_id] = row.sequence
            shuffled = list(nucleotide)
            rng.shuffle(shuffled)
            handle.write(f">negative__{row.reference_id}\n{''.join(shuffled)}\n")
    control_output = paths["results"] / "homology_controls_raw.tsv"
    subprocess.run(
        [
            "diamond", "blastx", "--db", str(database), "--query", str(control_fasta),
            "--out", str(control_output), "--outfmt", "6", "qseqid", "sseqid",
            "pident", "length", "qlen", "slen", "sstart", "send", "evalue", "bitscore",
            "--max-target-seqs", "10", "--evalue", "1e-5", "--threads", "2",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    control_hits = _qualified(_parse_diamond(control_output))
    observed_queries = set(control_hits["sequence_id"])
    true_positive = 0
    for query_id, expected_sequence in expected.items():
        query_hits = control_hits[control_hits["sequence_id"] == query_id]
        if any(
            reference_sequences.get(reference_id) == expected_sequence
            for reference_id in query_hits["reference_id"]
        ):
            true_positive += 1
    false_positive = sum(
        query_id.startswith("negative__") for query_id in observed_queries
    )
    control_validation = {
        "positive_controls": len(expected),
        "expected_reference_recovered": true_positive,
        "exact_sequence_pipeline_recovery_rate": (
            true_positive / len(expected) if expected else None
        ),
        "negative_controls": len(expected),
        "artificial_decoys_with_candidate_hit": false_positive,
        "artificial_decoy_hit_rate": (
            false_positive / len(expected) if expected else None
        ),
        "control_scope": (
            "Pipeline plumbing check using exact codon-encoded references and "
            "shuffled artificial decoys; not biological sensitivity or specificity."
        ),
        "thresholds": HOMOLOGY_THRESHOLDS,
    }
    (paths["results"] / "homology_control_validation.json").write_text(
        json.dumps(control_validation, indent=2), encoding="utf-8"
    )
    run_summary = {
        "diamond_version": subprocess.run(
            ["diamond", "version"], capture_output=True, text=True, check=True
        ).stdout.strip(),
        "query_sequences": len(records),
        "raw_alignment_rows": len(raw),
        "candidate_alignment_rows": len(qualified),
        "sequences_with_candidate_similarity": int(qualified["sequence_id"].nunique()),
        "runtime_seconds": elapsed,
        "query_fasta_sha256": sha256_file(query_fasta),
        "reference_fasta_sha256": sha256_file(reference_fasta),
        "raw_output_sha256": sha256_file(raw_output),
        "thresholds": HOMOLOGY_THRESHOLDS,
    }
    (paths["results"] / "homology_run_summary.json").write_text(
        json.dumps(run_summary, indent=2), encoding="utf-8"
    )
    return qualified, control_validation, run_summary


def create_corrected_schema(connection: sqlite3.Connection, variant: str):
    connection.executescript(
        """
        CREATE TABLE sequences (
            sequence_id TEXT PRIMARY KEY,
            sequence TEXT NOT NULL,
            length INTEGER NOT NULL,
            gc_fraction REAL NOT NULL,
            run_accession TEXT NOT NULL
        ) WITHOUT ROWID;
        CREATE INDEX idx_sequences_length ON sequences(length);
        CREATE INDEX idx_sequences_gc ON sequences(gc_fraction);
        CREATE INDEX idx_sequences_run ON sequences(run_accession);
        """
    )
    if variant == "bioplastics_metadata":
        connection.executescript(
            """
            CREATE TABLE studies (
                run_accession TEXT PRIMARY KEY,
                study_accession TEXT,
                study_title TEXT,
                scientific_name TEXT,
                tax_id TEXT,
                library_source TEXT,
                library_strategy TEXT,
                library_layout TEXT,
                context_class TEXT NOT NULL,
                selection_term TEXT
            ) WITHOUT ROWID;
            CREATE TABLE biodegradation_references (
                reference_id TEXT PRIMARY KEY,
                sequence TEXT NOT NULL,
                sequence_sha256 TEXT NOT NULL,
                protein_length INTEGER,
                molecular_weight REAL,
                isoelectric_point REAL,
                instability_index REAL,
                gravy REAL,
                aromaticity REAL
            ) WITHOUT ROWID;
            CREATE TABLE reference_annotations (
                annotation_id INTEGER PRIMARY KEY,
                reference_id TEXT NOT NULL,
                source TEXT NOT NULL,
                organism TEXT,
                tax_id TEXT,
                plastic TEXT,
                plastic_class TEXT,
                enzyme_name TEXT,
                enzyme_family TEXT,
                genbank_id TEXT,
                evidence TEXT,
                evidence_score REAL,
                evidence_tier TEXT,
                year TEXT,
                doi TEXT,
                FOREIGN KEY(reference_id) REFERENCES biodegradation_references(reference_id)
            );
            CREATE TABLE homology_hits (
                sequence_id TEXT NOT NULL,
                reference_id TEXT NOT NULL,
                identity_pct REAL NOT NULL,
                aligned_aa INTEGER NOT NULL,
                query_length_nt INTEGER NOT NULL,
                subject_length_aa INTEGER NOT NULL,
                subject_start INTEGER NOT NULL,
                subject_end INTEGER NOT NULL,
                evalue REAL NOT NULL,
                bitscore REAL NOT NULL,
                subject_coverage REAL NOT NULL,
                hit_rank INTEGER NOT NULL,
                PRIMARY KEY(sequence_id, reference_id),
                FOREIGN KEY(sequence_id) REFERENCES sequences(sequence_id),
                FOREIGN KEY(reference_id) REFERENCES biodegradation_references(reference_id)
            ) WITHOUT ROWID;
            CREATE TABLE provenance (
                resource TEXT PRIMARY KEY,
                sha256 TEXT,
                role TEXT NOT NULL
            ) WITHOUT ROWID;
            CREATE INDEX idx_studies_context ON studies(context_class);
            CREATE INDEX idx_annotations_plastic ON reference_annotations(plastic);
            CREATE INDEX idx_annotations_class ON reference_annotations(plastic_class);
            CREATE INDEX idx_annotations_family ON reference_annotations(enzyme_family);
            CREATE INDEX idx_hits_sequence ON homology_hits(sequence_id);
            CREATE INDEX idx_hits_reference ON homology_hits(reference_id);
            CREATE INDEX idx_hits_evalue ON homology_hits(evalue);
            """
        )
    elif variant != "no_bioplastics_metadata":
        raise ValueError(variant)


def _sequence_digest(path: Path) -> tuple[int, int, str]:
    connection = sqlite3.connect(path)
    digest = hashlib.sha256()
    count = bases = 0
    for row in connection.execute(
        "SELECT sequence_id,sequence,length,gc_fraction,run_accession "
        "FROM sequences ORDER BY sequence_id"
    ):
        digest.update(json.dumps(row, separators=(",", ":")).encode())
        count += 1
        bases += row[2]
    connection.close()
    return count, bases, digest.hexdigest()


def build_corrected_databases(
    root: Path,
    records: list[tuple],
    studies: pd.DataFrame,
    references: pd.DataFrame,
    annotations: pd.DataFrame,
    hits: pd.DataFrame,
    repetitions: int = 3,
):
    paths = ensure_directories(root)
    measurements = []
    variants = ["no_bioplastics_metadata", "bioplastics_metadata"]
    for repetition in range(1, repetitions + 1):
        order = variants if repetition % 2 else list(reversed(variants))
        for variant in order:
            suffix = "" if repetition == repetitions else f".rep{repetition}"
            path = paths["databases"] / f"logan_{variant}{suffix}.sqlite"
            path.unlink(missing_ok=True)
            start = time.perf_counter()
            connection = sqlite3.connect(path)
            configure_connection(connection)
            connection.execute("PRAGMA foreign_keys=ON")
            create_corrected_schema(connection, variant)
            connection.executemany(
                "INSERT INTO sequences VALUES (?,?,?,?,?)", records
            )
            if variant == "bioplastics_metadata":
                study_columns = [
                    "run_accession", "study_accession", "study_title",
                    "scientific_name", "tax_id", "library_source",
                    "library_strategy", "library_layout", "context_class",
                    "selection_term",
                ]
                connection.executemany(
                    f"INSERT INTO studies ({','.join(study_columns)}) "
                    f"VALUES ({','.join('?' for _ in study_columns)})",
                    studies[study_columns].itertuples(index=False, name=None),
                )
                connection.executemany(
                    "INSERT INTO biodegradation_references VALUES (?,?,?,?,?,?,?,?,?)",
                    references[
                        [
                            "reference_id", "sequence", "sequence_sha256",
                            "protein_length", "molecular_weight", "isoelectric_point",
                            "instability_index", "gravy", "aromaticity",
                        ]
                    ].itertuples(index=False, name=None),
                )
                annotation_columns = [
                    "reference_id", "source", "organism", "tax_id", "plastic",
                    "plastic_class", "enzyme_name", "enzyme_family", "genbank_id",
                    "evidence", "evidence_score", "evidence_tier", "year", "doi",
                ]
                connection.executemany(
                    f"INSERT INTO reference_annotations ({','.join(annotation_columns)}) "
                    f"VALUES ({','.join('?' for _ in annotation_columns)})",
                    annotations[annotation_columns].replace("", None).itertuples(
                        index=False, name=None
                    ),
                )
                hit_columns = [
                    "sequence_id", "reference_id", "identity_pct", "aligned_aa",
                    "query_length_nt", "subject_length_aa", "subject_start",
                    "subject_end", "evalue", "bitscore", "subject_coverage",
                    "hit_rank",
                ]
                if not hits.empty:
                    connection.executemany(
                        f"INSERT INTO homology_hits ({','.join(hit_columns)}) "
                        f"VALUES ({','.join('?' for _ in hit_columns)})",
                        hits[hit_columns].itertuples(index=False, name=None),
                    )
                provenance = [
                    ("PlasticDB raw", sha256_file(PLASTICDB_RAW), "reference metadata"),
                    ("PlasticDB scored", sha256_file(PLASTICDB_SCORED), "evidence tiers"),
                    ("PAZy", sha256_file(PAZY_CSV), "curated proteins"),
                    (
                        "DIAMOND raw output",
                        sha256_file(paths["results"] / "logan_biodegradation_homology_raw.tsv"),
                        "sequence-level homology",
                    ),
                ]
                connection.executemany("INSERT INTO provenance VALUES (?,?,?)", provenance)
            connection.commit()
            integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
            page_count = connection.execute("PRAGMA page_count").fetchone()[0]
            page_size = connection.execute("PRAGMA page_size").fetchone()[0]
            connection.close()
            elapsed = time.perf_counter() - start
            measurements.append(
                {
                    "variant": variant,
                    "repetition": repetition,
                    "build_seconds": elapsed,
                    "records_per_second": len(records) / elapsed,
                    "database_mib": path.stat().st_size / (1024 * 1024),
                    "page_count": page_count,
                    "page_size": page_size,
                    "integrity_check": integrity,
                }
            )
            if repetition != repetitions:
                path.unlink()

    no_path = paths["databases"] / "logan_no_bioplastics_metadata.sqlite"
    yes_path = paths["databases"] / "logan_bioplastics_metadata.sqlite"
    no_digest = _sequence_digest(no_path)
    yes_digest = _sequence_digest(yes_path)
    if no_digest != yes_digest:
        raise AssertionError("Canonical sequence rows differ")
    validation = {
        "sequence_count": no_digest[0],
        "total_bases": no_digest[1],
        "canonical_sequence_rows_sha256": no_digest[2],
        "variants_identical": True,
        "accessions": ALL_ACCESSIONS,
        "plastic_context_accessions": list(PLASTIC_CONTEXT_ACCESSIONS),
        "background_control_accessions": BACKGROUND_ACCESSIONS,
        "reference_sequences": len(references),
        "reference_annotations": len(annotations),
        "candidate_similarity_rows": len(hits),
        "sequences_with_candidate_similarity": int(hits["sequence_id"].nunique()),
    }
    frame = pd.DataFrame(measurements)
    frame.to_csv(paths["results"] / "bioplastics_build_measurements.csv", index=False)
    (paths["results"] / "bioplastics_validation.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    return frame, validation


def prepare_bioplastics_benchmark(root: Path):
    downloads, studies = acquire_corrected_sample(root)
    records = load_corrected_records(root)
    references, annotations = build_reference_library(root)
    hits, controls, homology = run_homology(root, records, references)
    builds, validation = build_corrected_databases(
        root, records, studies, references, annotations, hits
    )
    summary = pd.DataFrame(
        [
            {
                "logan_accessions": len(ALL_ACCESSIONS),
                "plastic_context_accessions": len(PLASTIC_CONTEXT_ACCESSIONS),
                "sequence_records": len(records),
                "total_bases": sum(record[2] for record in records),
                "reference_sequences": len(references),
                "reference_annotations": len(annotations),
                "candidate_similarity_rows": len(hits),
            }
        ]
    )
    return (
        downloads, studies, references, annotations, hits, controls,
        homology, builds, validation, summary,
    )


def benchmark_bioplastics_databases(
    root: Path, repetitions: int = 30, warmups: int = 5
):
    paths = ensure_directories(root)
    dbs = {
        "no_bioplastics_metadata": (
            paths["databases"] / "logan_no_bioplastics_metadata.sqlite"
        ),
        "bioplastics_metadata": (
            paths["databases"] / "logan_bioplastics_metadata.sqlite"
        ),
    }
    common_queries = {
        "count_all_sequences": ("SELECT COUNT(*) FROM sequences", ()),
        "average_sequence_length": ("SELECT AVG(length) FROM sequences", ()),
        "gc_range_count": (
            "SELECT COUNT(*) FROM sequences WHERE gc_fraction BETWEEN ? AND ?",
            (0.45, 0.55),
        ),
        "top_100_longest": (
            "SELECT sequence_id,length FROM sequences ORDER BY length DESC LIMIT 100",
            (),
        ),
        "plastic_context_sequence_count": (
            "SELECT COUNT(*) FROM sequences WHERE run_accession IN (?,?,?)",
            tuple(PLASTIC_CONTEXT_ACCESSIONS),
        ),
    }
    raw_rows = []
    proofs = {}
    plans = {}
    for variant, path in dbs.items():
        connection = sqlite3.connect(path)
        connection.execute("PRAGMA cache_size=-32768")
        plans[variant] = {}
        for name, (sql, params) in common_queries.items():
            samples, result = timed_query(
                connection, sql, params, repetitions, warmups
            )
            encoded = json.dumps(result, separators=(",", ":")).encode()
            proofs.setdefault(name, {})[variant] = hashlib.sha256(encoded).hexdigest()
            plans[variant][name] = connection.execute(
                f"EXPLAIN QUERY PLAN {sql}", params
            ).fetchall()
            raw_rows.extend(
                {
                    "variant": variant,
                    "query": name,
                    "repetition": index + 1,
                    "latency_ms": sample,
                }
                for index, sample in enumerate(samples)
            )
        connection.close()
    persisted_proofs = {}
    for name, variants in proofs.items():
        equivalent = (
            variants["no_bioplastics_metadata"]
            == variants["bioplastics_metadata"]
        )
        if not equivalent:
            raise AssertionError(f"Shared query mismatch: {name}")
        persisted_proofs[name] = dict(variants, equivalent=True)

    raw = pd.DataFrame(raw_rows)
    summary = (
        raw.groupby(["variant", "query"])["latency_ms"]
        .agg(
            median_ms="median", mean_ms="mean", stddev_ms="std",
            minimum_ms="min", maximum_ms="max",
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
    raw.to_csv(paths["results"] / "bioplastics_query_measurements.csv", index=False)
    summary.to_csv(paths["results"] / "bioplastics_query_summary.csv", index=False)
    (paths["results"] / "bioplastics_shared_query_proofs.json").write_text(
        json.dumps(persisted_proofs, indent=2), encoding="utf-8"
    )
    (paths["results"] / "bioplastics_query_plans.json").write_text(
        json.dumps(plans, indent=2), encoding="utf-8"
    )
    return raw, summary, persisted_proofs, plans


def benchmark_bioplastics_capabilities(
    root: Path, repetitions: int = 30, warmups: int = 5
):
    paths = ensure_directories(root)
    connection = sqlite3.connect(
        paths["databases"] / "logan_bioplastics_metadata.sqlite"
    )
    queries = {
        "candidate_sequences_by_reference_plastic": """
            SELECT a.plastic, COUNT(DISTINCT h.sequence_id)
            FROM homology_hits h
            JOIN reference_annotations a USING(reference_id)
            GROUP BY a.plastic ORDER BY COUNT(DISTINCT h.sequence_id) DESC
        """,
        "candidate_sequences_by_reference_enzyme_family": """
            SELECT a.enzyme_family, COUNT(DISTINCT h.sequence_id)
            FROM homology_hits h
            JOIN reference_annotations a USING(reference_id)
            GROUP BY a.enzyme_family ORDER BY COUNT(DISTINCT h.sequence_id) DESC
        """,
        "candidate_sequences_matching_bioplastic_or_biodegradable_references": """
            SELECT COUNT(DISTINCT h.sequence_id)
            FROM homology_hits h
            JOIN reference_annotations a USING(reference_id)
            WHERE a.plastic_class='bioplastic_or_biodegradable_polymer'
        """,
        "hits_by_study_context": """
            SELECT s.context_class, COUNT(DISTINCT h.sequence_id)
            FROM homology_hits h
            JOIN sequences q USING(sequence_id)
            JOIN studies s USING(run_accession)
            GROUP BY s.context_class ORDER BY s.context_class
        """,
        "reference_biochemistry_by_plastic_class": """
            SELECT a.plastic_class, COUNT(DISTINCT r.reference_id),
                   AVG(r.molecular_weight), AVG(r.isoelectric_point),
                   AVG(r.instability_index), AVG(r.gravy)
            FROM biodegradation_references r
            JOIN reference_annotations a USING(reference_id)
            GROUP BY a.plastic_class ORDER BY a.plastic_class
        """,
        "reference_coverage_by_source": """
            SELECT a.source, COUNT(DISTINCT a.reference_id), COUNT(*)
            FROM reference_annotations a GROUP BY a.source ORDER BY a.source
        """,
    }
    rows = []
    plans = {}
    for name, sql in queries.items():
        samples, result = timed_query(
            connection, sql, (), repetitions, warmups
        )
        plans[name] = connection.execute(f"EXPLAIN QUERY PLAN {sql}").fetchall()
        rows.append(
            {
                "capability": name,
                "median_ms": pd.Series(samples).median(),
                "p95_ms": percentile(samples, 0.95),
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
    frame.to_csv(
        paths["results"] / "bioplastics_metadata_capabilities.csv", index=False
    )
    (paths["results"] / "bioplastics_capability_query_plans.json").write_text(
        json.dumps(plans, indent=2), encoding="utf-8"
    )
    return frame, plans


def bioplastics_storage_summary(root: Path):
    paths = ensure_directories(root)
    rows = []
    for variant in ["no_bioplastics_metadata", "bioplastics_metadata"]:
        path = paths["databases"] / f"logan_{variant}.sqlite"
        connection = sqlite3.connect(path)
        row = {
            "variant": variant,
            "database_mib": path.stat().st_size / (1024 * 1024),
            "page_count": connection.execute("PRAGMA page_count").fetchone()[0],
            "page_size": connection.execute("PRAGMA page_size").fetchone()[0],
            "integrity_check": connection.execute(
                "PRAGMA integrity_check"
            ).fetchone()[0],
            "sequence_rows": connection.execute(
                "SELECT COUNT(*) FROM sequences"
            ).fetchone()[0],
            "reference_rows": 0,
            "annotation_rows": 0,
            "homology_hit_rows": 0,
        }
        if variant == "bioplastics_metadata":
            row["reference_rows"] = connection.execute(
                "SELECT COUNT(*) FROM biodegradation_references"
            ).fetchone()[0]
            row["annotation_rows"] = connection.execute(
                "SELECT COUNT(*) FROM reference_annotations"
            ).fetchone()[0]
            row["homology_hit_rows"] = connection.execute(
                "SELECT COUNT(*) FROM homology_hits"
            ).fetchone()[0]
        rows.append(row)
        connection.close()
    frame = pd.DataFrame(rows)
    frame.to_csv(
        paths["results"] / "bioplastics_storage_summary.csv", index=False
    )
    return frame