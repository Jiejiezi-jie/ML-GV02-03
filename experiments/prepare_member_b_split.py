"""Create an isolated development split from the pinned conservative release."""
from __future__ import annotations

import argparse
from collections import defaultdict, Counter
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.generation_data import cluster_aware_split
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
    assigned, _ = cluster_aware_split(clusters, {"train":0.8,"validation":0.1,"test":0.1}, seed=42)
    output = directory / "development_split"
    paths = []
    for split in ("train", "validation", "test"):
        path = output / (split + ".fasta")
        write_fasta([r for r in records if assigned[r.identifier] == split], path)
        paths.append(path)
    rows = [dict(**audit[r.identifier], split=assigned[r.identifier]) for r in records]
    write_tsv(output / "sequence_manifest.tsv", rows, list(rows[0]))
    tokens = ["<PAD>","<BOS>","<EOS>","<UNK>"] + list("ACDEFGHIKLMNPQRSTVWY")
    write_json(output / "vocabulary.json", dict(tokens=tokens, token_to_id={s:i for i,s in enumerate(tokens)}))
    paths += [output / "sequence_manifest.tsv", output / "vocabulary.json"]
    summary = dict(dataset_status="provisional_development_only", seed=42,
        release_id=release["release_id"], release_sha256=sha256_file(manifest_path),
        total_sequences=len(records), total_clusters=len(clusters),
        splits=dict(Counter(assigned.values())), cluster_identity=0.7, cluster_bidirectional_coverage=0.8,
        algorithm="seeded_largest_cluster_first_on_filtered_original_70_percent_clusters",
        warning="Do not call this the final training set. Rebuild when A release changes; test set must not tune model or sampling.",
        outputs={p.relative_to(ROOT).as_posix():sha256_file(p) for p in paths})
    write_json(output / "split_manifest.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--allow-provisional", action="store_true", help="Explicitly produce a development-only split")
    args = parser.parse_args()
    prepare(args.allow_provisional)
