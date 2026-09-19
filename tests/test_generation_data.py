from __future__ import annotations

import csv
import json
from pathlib import Path

from gv_eval.generation_data import cluster_aware_split, parse_cdhit_clusters
from gv_eval.io import read_fasta, sha256_file


ROOT = Path(__file__).resolve().parents[1]


def test_parse_cdhit_clusters(tmp_path):
    path = tmp_path / "toy.clstr"
    path.write_text(
        ">Cluster 0\n"
        "0 70aa, >a... *\n"
        "1 69aa, >b... at 95.00%\n"
        ">Cluster 1\n"
        "0 80aa, >c... *\n",
        encoding="utf-8",
    )
    assert parse_cdhit_clusters(path) == [["a", "b"], ["c"]]


def test_cluster_split_is_deterministic_complete_and_leak_free():
    clusters = [
        ["a", "b", "c"],
        ["d", "e"],
        ["f"],
        ["g"],
        ["h"],
        ["i"],
        ["j"],
    ]
    ratios = {"train": 0.6, "validation": 0.2, "test": 0.2}
    first, first_cluster_ids = cluster_aware_split(clusters, ratios, seed=42)
    second, second_cluster_ids = cluster_aware_split(clusters, ratios, seed=42)

    assert first == second
    assert first_cluster_ids == second_cluster_ids
    assert set(first) == set("abcdefghij")
    assert set(first.values()) == {"train", "validation", "test"}
    assert all(len({first[member] for member in cluster}) == 1 for cluster in clusters)


def test_cluster_split_rejects_invalid_ratios():
    clusters = [["a"], ["b"], ["c"]]
    try:
        cluster_aware_split(
            clusters,
            {"train": 0.8, "validation": 0.1, "test": 0.2},
            seed=42,
        )
    except ValueError as error:
        assert "sum to 1" in str(error)
    else:
        raise AssertionError("invalid ratios were accepted")


def test_generator_data_artifacts_are_complete_isolated_and_current():
    output_dir = ROOT / "data/processed/gv02_03_v2/generator_data"
    manifest = json.loads((output_dir / "split_manifest.json").read_text(encoding="utf-8"))
    assert manifest["total_sequences"] == 1224
    assert manifest["total_clusters"] == 418
    assert manifest["family_verification_status"] == (
        "pending_competitive_gvpa_gvpj_verification"
    )

    with (output_dir / "sequence_manifest.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 1224
    assert len({row["sequence_id"] for row in rows}) == 1224
    cluster_splits: dict[str, set[str]] = {}
    for row in rows:
        cluster_splits.setdefault(row["cluster_id"], set()).add(row["split"])
    assert len(cluster_splits) == 418
    assert all(len(splits) == 1 for splits in cluster_splits.values())

    fasta_ids = {
        split: {record.identifier for record in read_fasta(output_dir / f"{split}.fasta")}
        for split in ("train", "validation", "test")
    }
    assert [len(fasta_ids[name]) for name in fasta_ids] == [978, 123, 123]
    assert not (fasta_ids["train"] & fasta_ids["validation"])
    assert not (fasta_ids["train"] & fasta_ids["test"])
    assert not (fasta_ids["validation"] & fasta_ids["test"])

    for relative, expected in manifest["input_sha256"].items():
        assert sha256_file(ROOT / relative) == expected
    for relative, expected in manifest["output_sha256"].items():
        assert sha256_file(ROOT / relative) == expected
