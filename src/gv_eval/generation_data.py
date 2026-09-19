from __future__ import annotations

import csv
import random
import re
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import yaml

from .io import FastaRecord, read_fasta, sha256_file, write_fasta, write_json


_CDHIT_MEMBER = re.compile(r">(.+?)\.\.\.")
_SPLIT_NAMES = ("train", "validation", "test")


def parse_cdhit_clusters(path: str | Path) -> list[list[str]]:
    """Parse CD-HIT's .clstr format while preserving cluster membership."""

    clusters: list[list[str]] = []
    current: list[str] | None = None
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if line.startswith(">Cluster "):
            current = []
            clusters.append(current)
            continue
        if not line:
            continue
        if current is None:
            raise ValueError(f"{path}:{line_number}: member before cluster header")
        match = _CDHIT_MEMBER.search(line)
        if not match:
            raise ValueError(f"{path}:{line_number}: malformed CD-HIT member row")
        current.append(match.group(1))

    if not clusters or any(not cluster for cluster in clusters):
        raise ValueError(f"CD-HIT cluster file is empty or contains an empty cluster: {path}")
    members = [member for cluster in clusters for member in cluster]
    duplicates = sorted(member for member, count in Counter(members).items() if count > 1)
    if duplicates:
        raise ValueError(f"CD-HIT members occur in multiple clusters: {duplicates[:5]}")
    return clusters


def cluster_aware_split(
    clusters: Sequence[Sequence[str]],
    ratios: Mapping[str, float],
    *,
    seed: int,
) -> tuple[dict[str, str], dict[str, list[int]]]:
    """Assign whole homology clusters to splits with near-target sequence counts."""

    if tuple(ratios) != _SPLIT_NAMES:
        raise ValueError(f"Split names and order must be {_SPLIT_NAMES}")
    if any(value <= 0 for value in ratios.values()) or abs(sum(ratios.values()) - 1.0) > 1e-9:
        raise ValueError("Split ratios must be positive and sum to 1")
    if len(clusters) < len(_SPLIT_NAMES):
        raise ValueError("At least three homology clusters are required")

    total = sum(len(cluster) for cluster in clusters)
    targets = {name: ratios[name] * total for name in _SPLIT_NAMES}
    counts = {name: 0 for name in _SPLIT_NAMES}
    cluster_ids = {name: [] for name in _SPLIT_NAMES}
    member_to_split: dict[str, str] = {}

    rng = random.Random(seed)
    decorated = [(rng.random(), index, list(cluster)) for index, cluster in enumerate(clusters)]
    decorated.sort(key=lambda item: (-len(item[2]), item[0], item[1]))

    for _, cluster_index, members in decorated:
        split = max(
            _SPLIT_NAMES,
            key=lambda name: (
                (targets[name] - counts[name]) / targets[name],
                -_SPLIT_NAMES.index(name),
            ),
        )
        cluster_ids[split].append(cluster_index)
        counts[split] += len(members)
        for member in members:
            if member in member_to_split:
                raise ValueError(f"Duplicate cluster member: {member}")
            member_to_split[member] = split

    if any(not cluster_ids[name] for name in _SPLIT_NAMES):
        raise ValueError("Cluster assignment produced an empty split")
    return member_to_split, cluster_ids


def _read_reference_mapping(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    mapping: dict[str, dict[str, str]] = {}
    for row in rows:
        reference_id = row.get("reference_id", "")
        if reference_id and row.get("unique_representative") == "True":
            if reference_id in mapping:
                raise ValueError(f"Duplicate unique mapping row for {reference_id}")
            mapping[reference_id] = row
    return mapping


def prepare_generator_data(root: str | Path, config_path: str | Path) -> dict:
    """Build deterministic, cluster-isolated data files for generator training."""

    root = Path(root).resolve()
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = root / config_path
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    data_config = config["data"]
    split_config = data_config["split"]
    if split_config.get("algorithm") != "seeded_largest_cluster_first":
        raise ValueError("Unsupported split algorithm")
    ratios = {name: float(split_config[name]) for name in _SPLIT_NAMES}

    reference_path = root / data_config["reference_fasta"]
    mapping_path = root / data_config["reference_mapping"]
    cluster_path = root / data_config["cluster_file"]
    records_list = list(read_fasta(reference_path))
    records = {record.identifier: record for record in records_list}
    if not records or len(records) != len(records_list):
        raise ValueError("Reference FASTA is empty or contains duplicate identifiers")

    clusters = parse_cdhit_clusters(cluster_path)
    clustered_ids = {member for cluster in clusters for member in cluster}
    if clustered_ids != set(records):
        missing = sorted(set(records) - clustered_ids)
        unexpected = sorted(clustered_ids - set(records))
        raise ValueError(
            f"Reference/cluster membership mismatch; missing={missing[:5]}, "
            f"unexpected={unexpected[:5]}"
        )

    member_to_split, split_cluster_ids = cluster_aware_split(
        clusters, ratios, seed=int(config["seed"])
    )
    source_mapping = _read_reference_mapping(mapping_path)
    if set(source_mapping) != set(records):
        raise ValueError("Reference mapping does not cover every cleaned reference exactly once")

    output_dir = root / config["outputs"]["data_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_paths: dict[str, Path] = {}
    for split in _SPLIT_NAMES:
        fasta_paths[split] = output_dir / f"{split}.fasta"
        split_records = []
        for record in records_list:
            if member_to_split[record.identifier] != split:
                continue
            description = record.description
            identifier_prefix = f"{record.identifier} "
            if description.startswith(identifier_prefix):
                description = description[len(identifier_prefix) :]
            split_records.append(FastaRecord(record.identifier, description, record.sequence))
        write_fasta(
            split_records,
            fasta_paths[split],
        )

    cluster_by_member = {
        member: cluster_index
        for cluster_index, cluster in enumerate(clusters)
        for member in cluster
    }
    audit_rows = []
    for record in records_list:
        source = source_mapping[record.identifier]
        source_path = Path(source["source_path"])
        try:
            source_path_value = source_path.resolve().relative_to(root).as_posix()
        except ValueError:
            source_path_value = source_path.as_posix()
        audit_rows.append(
            {
                "sequence_id": record.identifier,
                "split": member_to_split[record.identifier],
                "cluster_id": f"cluster_{cluster_by_member[record.identifier]:06d}",
                "length": len(record.sequence),
                "sequence_sha256": source["sequence_sha256"],
                "source_path": source_path_value,
                "source_id": source["source_id"],
                "source_header": source["source_header"],
                "annotation_support": "explicit_gvpa_source_label",
                "family_verification_status": data_config["verification_status"],
            }
        )
    audit_path = output_dir / "sequence_manifest.tsv"
    with audit_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(audit_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(audit_rows)

    special_tokens = list(config["vocabulary"]["special_tokens"])
    amino_acids = list(config["vocabulary"]["amino_acids"])
    tokens = special_tokens + amino_acids
    if len(tokens) != len(set(tokens)):
        raise ValueError("Vocabulary contains duplicate tokens")
    vocabulary = {
        "tokens": tokens,
        "token_to_id": {token: index for index, token in enumerate(tokens)},
        "special_tokens": special_tokens,
        "amino_acids": amino_acids,
    }
    vocabulary_path = root / config["outputs"]["vocabulary"]
    write_json(vocabulary_path, vocabulary)

    split_summary = {}
    for split in _SPLIT_NAMES:
        lengths = [len(record.sequence) for record in records_list if member_to_split[record.identifier] == split]
        split_summary[split] = {
            "sequences": len(lengths),
            "clusters": len(split_cluster_ids[split]),
            "fraction": len(lengths) / len(records_list),
            "minimum_length": min(lengths),
            "maximum_length": max(lengths),
        }

    output_hashes = {
        fasta_paths[name].relative_to(root).as_posix(): sha256_file(fasta_paths[name])
        for name in _SPLIT_NAMES
    }
    output_hashes[audit_path.relative_to(root).as_posix()] = sha256_file(audit_path)
    output_hashes[vocabulary_path.relative_to(root).as_posix()] = sha256_file(vocabulary_path)
    manifest = {
        "version": config["version"],
        "seed": int(config["seed"]),
        "dataset_status": "provisional_nominal_gvpa",
        "family_verification_status": data_config["verification_status"],
        "warning": (
            "The source annotation, length and character filters are complete, but "
            "competitive GvpA/GvpJ verification is still pending. These files are not "
            "the final high-confidence GvpA training set."
        ),
        "split_algorithm": split_config["algorithm"],
        "cluster_identity": float(data_config["cluster_identity"]),
        "ratios": ratios,
        "total_sequences": len(records_list),
        "total_clusters": len(clusters),
        "splits": split_summary,
        "input_sha256": {
            config_path.relative_to(root).as_posix(): sha256_file(config_path),
            reference_path.relative_to(root).as_posix(): sha256_file(reference_path),
            mapping_path.relative_to(root).as_posix(): sha256_file(mapping_path),
            cluster_path.relative_to(root).as_posix(): sha256_file(cluster_path),
        },
        "output_sha256": dict(sorted(output_hashes.items())),
    }
    manifest_path = output_dir / "split_manifest.json"
    write_json(manifest_path, manifest)
    return manifest
