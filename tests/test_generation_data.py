from __future__ import annotations

import csv
import json
from pathlib import Path

from gv_eval.generation_data import (
    cluster_aware_split,
    constrained_cluster_split,
    parse_cdhit_clusters,
)
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


def test_constrained_cluster_split_matches_counts_lengths_and_is_deterministic():
    clusters = [["a", "b"], ["c"], ["d"], ["e"], ["f"], ["g"]]
    lengths = {"a": 4, "b": 6, "c": 8, "d": 9, "e": 10, "f": 11, "g": 12}
    sequence_targets = {"train": 3, "validation": 2, "test": 2}
    cluster_targets = {"train": 2, "validation": 2, "test": 2}
    length_targets = {"train": 22, "validation": 17, "test": 21}

    first, first_clusters = constrained_cluster_split(
        clusters,
        lengths,
        sequence_targets,
        cluster_targets,
        length_targets,
        seed=42,
    )
    second, second_clusters = constrained_cluster_split(
        clusters,
        lengths,
        sequence_targets,
        cluster_targets,
        length_targets,
        seed=42,
    )

    assert first == second
    assert first_clusters == second_clusters
    assert all(len({first[member] for member in cluster}) == 1 for cluster in clusters)
    for split in ("train", "validation", "test"):
        members = [member for member, assigned in first.items() if assigned == split]
        assert len(members) == sequence_targets[split]
        assert len(first_clusters[split]) == cluster_targets[split]
        assert sum(lengths[member] for member in members) == length_targets[split]
