"""Create an isolated development split from the pinned conservative release."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.generation_data import constrained_cluster_split
from gv_eval.io import read_fasta, write_fasta, write_json, write_tsv, sha256_file
from prepare_member_b_handoff import read_table


def prepare(allow_provisional=False):
    directory = ROOT / "data/processed/gv02_03_v2/member_b_handoff"
    manifest_path = directory / "release_manifest.json"
    release = json.loads(manifest_path.read_text(encoding="utf-8"))
    if release["status"] != "final" and not allow_provisional:
        raise ValueError("Reference is not final. Development only: pass --allow-provisional and label all results provisional.")
    for section in ("inputs", "outputs"):
        for path, digest in release[section].items():
            if sha256_file(ROOT / path) != digest:
                raise ValueError(f"Release hash mismatch: {path}")
    records = list(read_fasta(directory / "gvpa_training_eligible.fasta"))
    audit = {r["sequence_id"]:r for r in read_table(directory / "release_audit.tsv")}
    groups = defaultdict(list)
    for record in records:
        groups[audit[record.identifier]["cluster_id"]].append(record.identifier)
    # Preserve original mixed-family 70% clusters after filtering. Reassign
    # groups afresh; do not recycle generator or classification fold labels.
    clusters = [sorted(groups[k]) for k in sorted(groups)]
    member_lengths = {record.identifier: len(record.sequence) for record in records}
    target_sequences = {"train": 277, "validation": 33, "test": 36}
    target_clusters = {"train": 16, "validation": 9, "test": 6}
    target_total_lengths = {"train": 30151, "validation": 3590, "test": 4021}
    assigned, split_cluster_ids = constrained_cluster_split(
        clusters,
        member_lengths,
        target_sequences,
        target_clusters,
        target_total_lengths,
        seed=42,
    )
    output = directory / "development_split"
    paths = []
    for split in ("train", "validation", "test"):
        path = output / (split + ".fasta")
        write_fasta([r for r in records if assigned[r.identifier] == split], path)
        paths.append(path)
    rows = [dict(**audit[r.identifier], split=assigned[r.identifier]) for r in records]
    write_tsv(output / "sequence_manifest.tsv", rows, list(rows[0]))
    special_tokens = ["<PAD>", "<BOS>", "<EOS>", "<UNK>"]
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    tokens = special_tokens + amino_acids
    write_json(
        output / "vocabulary.json",
        dict(
            tokens=tokens,
            token_to_id={symbol: index for index, symbol in enumerate(tokens)},
            special_tokens=special_tokens,
            amino_acids=amino_acids,
        ),
    )
    paths += [output / "sequence_manifest.tsv", output / "vocabulary.json"]
    split_statistics = {}
    for split in ("train", "validation", "test"):
        lengths = [member_lengths[identifier] for identifier, value in assigned.items() if value == split]
        split_statistics[split] = dict(
            sequences=len(lengths),
            clusters=len(split_cluster_ids[split]),
            fraction=len(lengths) / len(records),
            mean_length=sum(lengths) / len(lengths),
            minimum_length=min(lengths),
            maximum_length=max(lengths),
            total_length=sum(lengths),
        )
    summary = dict(dataset_status="provisional_development_only", seed=42,
        release_id=release["release_id"], release_sha256=sha256_file(manifest_path),
        total_sequences=len(records), total_clusters=len(clusters),
        splits=split_statistics, cluster_identity=0.7, cluster_bidirectional_coverage=0.8,
        targets=dict(sequences=target_sequences, clusters=target_clusters, total_lengths=target_total_lengths),
        algorithm="constrained_cluster_stratification_by_count_and_length",
        warning="Do not call this the final training set. Rebuild when A release changes; test set must not tune model or sampling.",
        outputs={p.relative_to(ROOT).as_posix():sha256_file(p) for p in paths})
    write_json(output / "split_manifest.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-provisional", action="store_true", help="Explicitly produce a development-only split")
    args = parser.parse_args()
    prepare(args.allow_provisional)
