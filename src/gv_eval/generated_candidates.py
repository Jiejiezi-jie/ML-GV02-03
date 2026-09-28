from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import yaml

from .io import FastaRecord, read_fasta, sequence_sha256, sha256_file, write_fasta, write_json, write_tsv
from .quality import assess_candidate_files


def _read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _truth(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _load(root: Path, config_path: str | Path) -> tuple[Path, dict[str, Any]]:
    path = Path(config_path)
    if not path.is_absolute():
        path = root / path
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError("Generated-candidate config must be a mapping")
    return path, config


def prepare_generated_qc(root: str | Path, config_path: str | Path) -> dict[str, Any]:
    """Audit generation metadata and write an immutable pre-family QC split."""

    root = Path(root).resolve()
    config_path, config = _load(root, config_path)
    inputs = config["inputs"]
    outputs = config["outputs"]
    candidate_path = root / inputs["candidates"]
    metadata_path = root / inputs["generation_metadata"]
    training_path = root / inputs["training_fasta"]
    reference_path = root / inputs["reference_fasta"]
    records = list(read_fasta(candidate_path))
    metadata = _read_table(metadata_path)
    if not records or len({record.identifier for record in records}) != len(records):
        raise ValueError("Generated candidates must be nonempty with unique identifiers")
    if [record.identifier for record in records] != [row.get("sequence_id") for row in metadata]:
        raise ValueError("Generation metadata IDs/order do not match candidate FASTA")

    for record, row in zip(records, metadata, strict=True):
        if int(row["sequence_length"]) != len(record.sequence):
            raise ValueError(f"Generation metadata length mismatch: {record.identifier}")
        if row["sequence_sha256"] != sequence_sha256(record.sequence):
            raise ValueError(f"Generation metadata sequence hash mismatch: {record.identifier}")
        if _truth(row["terminated_by_eos"]) == _truth(row["hit_generation_cap"]):
            raise ValueError(f"Generation stop flags are not mutually exclusive: {record.identifier}")

    quality = assess_candidate_files(
        records,
        training_fasta=training_path,
        reference_fasta=reference_path,
        **config["quality"],
    )
    rows: list[dict[str, object]] = []
    for record, generation in zip(records, metadata, strict=True):
        result = dict(quality[record.identifier])
        cap_warning = _truth(generation["hit_generation_cap"])
        reasons = [reason for reason in str(result["qc_reasons"]).split(";") if reason]
        if cap_warning and "generation_cap_warning" not in reasons:
            reasons.append("generation_cap_warning")
        result["generation_cap_warning"] = cap_warning
        result["qc_reasons"] = ";".join(reasons)
        rows.append(
            {
                "sequence_id": record.identifier,
                "sequence_sha256": sequence_sha256(record.sequence),
                "terminated_by_eos": generation["terminated_by_eos"],
                "hit_generation_cap": generation["hit_generation_cap"],
                "stop_reason": generation["stop_reason"],
                **result,
            }
        )

    output_dir = root / outputs["processed_dir"]
    qc_path = output_dir / "candidate_qc.tsv"
    passing_path = output_dir / "qc_pass.fasta"
    summary_path = output_dir / "qc_summary.json"
    manifest_path = output_dir / "qc_manifest.json"
    write_tsv(qc_path, rows, list(rows[0]))
    passing_ids = {str(row["sequence_id"]) for row in rows if bool(row["qc_pass"])}
    passing = [record for record in records if record.identifier in passing_ids]
    write_fasta(passing, passing_path)
    reason_counts = Counter(
        reason
        for row in rows
        for reason in str(row["qc_reasons"]).split(";")
        if reason
    )
    summary = {
        "input_count": len(records),
        "pass_count": len(passing),
        "failure_count": len(records) - len(passing),
        "natural_eos_count": sum(_truth(row["terminated_by_eos"]) for row in rows),
        "generation_cap_count": sum(_truth(row["hit_generation_cap"]) for row in rows),
        "reason_counts": dict(sorted(reason_counts.items())),
    }
    write_json(summary_path, summary)
    write_json(
        manifest_path,
        {
            "version": config.get("version", 1),
            "seed": config["seed"],
            "config": config,
            "inputs": {
                path.relative_to(root).as_posix(): sha256_file(path)
                for path in (config_path, candidate_path, metadata_path, training_path, reference_path)
            },
            "outputs": {
                path.relative_to(root).as_posix(): sha256_file(path)
                for path in (qc_path, passing_path, summary_path)
            },
        },
    )
    return summary


def build_generated_candidate_pools(root: str | Path, config_path: str | Path) -> dict[str, Any]:
    """Join QC and competitive family evidence into disjoint candidate pools."""

    root = Path(root).resolve()
    config_path, config = _load(root, config_path)
    inputs = config["inputs"]
    output_dir = root / config["outputs"]["processed_dir"]
    candidate_path = root / inputs["candidates"]
    qc_path = output_dir / "candidate_qc.tsv"
    family_path = root / inputs["family_classification"]
    records = list(read_fasta(candidate_path))
    qc_rows = _read_table(qc_path)
    family_rows = _read_table(family_path)
    if [record.identifier for record in records] != [row["sequence_id"] for row in qc_rows]:
        raise ValueError("QC rows do not match generated candidate IDs/order")
    qc_by_id = {row["sequence_id"]: row for row in qc_rows}
    expected_family_ids = [row["sequence_id"] for row in qc_rows if _truth(row["qc_pass"])]
    if expected_family_ids != [row.get("sequence_id") for row in family_rows]:
        raise ValueError("Family rows do not match QC-passing candidate IDs/order")
    family_by_id = {row["sequence_id"]: row for row in family_rows}

    pool_records: dict[str, list[FastaRecord]] = {
        "main_supported_gvpa": [],
        "ambiguous_exploration": [],
        "excluded": [],
    }
    audit_rows: list[dict[str, object]] = []
    for record in records:
        qc = qc_by_id[record.identifier]
        family = family_by_id.get(record.identifier, {})
        if not _truth(qc["qc_pass"]):
            pool = "excluded"
            pool_reason = "qc_failed"
        elif family.get("family_status") == "supported_gvpa" and family.get("pf00741_status") == "pass":
            pool = "main_supported_gvpa"
            pool_reason = "qc_and_competitive_gvpa_and_pf00741_pass"
        elif family.get("family_status") == "gvpa_gvpj_ambiguous":
            pool = "ambiguous_exploration"
            pool_reason = "qc_pass_but_gvpa_gvpj_ambiguous"
        else:
            pool = "excluded"
            pool_reason = "qc_pass_but_not_supported_gvpa"
        if family and family.get("sequence_sha256") != sequence_sha256(record.sequence):
            raise ValueError(f"Family classification sequence hash mismatch: {record.identifier}")
        pool_records[pool].append(record)
        audit_rows.append(
            {
                **qc,
                **{f"family_{key}": value for key, value in family.items() if key != "sequence_id"},
                "pool": pool,
                "pool_reason": pool_reason,
            }
        )

    output_paths = {
        name: output_dir / f"{name}.fasta" for name in pool_records
    }
    for name, path in output_paths.items():
        write_fasta(pool_records[name], path)
    audit_path = output_dir / "candidate_pool_audit.tsv"
    write_tsv(audit_path, audit_rows, list(audit_rows[0]))
    summary = {
        "input_count": len(records),
        "pool_counts": {name: len(rows) for name, rows in pool_records.items()},
        "family_status_counts": dict(
            sorted(Counter(row.get("family_status", "qc_failed") for row in family_rows).items())
        ),
    }
    summary_path = output_dir / "candidate_pool_summary.json"
    manifest_path = output_dir / "candidate_pool_manifest.json"
    write_json(summary_path, summary)
    write_json(
        manifest_path,
        {
            "version": config.get("version", 1),
            "seed": config["seed"],
            "config": config,
            "inputs": {
                path.relative_to(root).as_posix(): sha256_file(path)
                for path in (config_path, candidate_path, qc_path, family_path)
            },
            "outputs": {
                path.relative_to(root).as_posix(): sha256_file(path)
                for path in (*output_paths.values(), audit_path, summary_path)
            },
        },
    )
    return summary
