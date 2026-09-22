from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import yaml

from gv_eval import similarity
from gv_eval.io import FastaRecord, sha256_file


def rec(identifier, sequence):
    return FastaRecord(identifier, identifier, sequence)


def settings(**kwargs):
    return similarity.SimilarityConfig(minimum_aligned_residues=1, **kwargs)


def test_identical_and_single_substitution_have_hand_checked_distances():
    config = settings()
    same = similarity.align_pair(rec("q", "ACDEFGHIK"), rec("r", "acdefghik"), config)
    assert same["distance_status"] == "resolved"
    assert same["distance"] == 0
    assert same["identity"] == 1
    changed = similarity.align_pair(rec("q", "ACDEFGHIK"), rec("r", "ACDEYGHIK"), config)
    assert changed["identical_residues"] == 8
    assert changed["aligned_pairs"] == changed["alignment_length"] == 9
    assert changed["query_coverage"] == changed["target_coverage"] == 1
    assert changed["identity"] == pytest.approx(8 / 9)
    assert changed["distance"] == pytest.approx(1 / 9)


def test_terminal_overhang_has_partial_coverage_and_unresolved_distance():
    q, t = rec("q", "ACDEFGHIK"), rec("t", "ACDEFGHIKWWWWWWWWWWWW")
    result = similarity.align_pair(q, t, settings())
    assert result["identical_residues"] == result["aligned_pairs"] == 9
    assert result["alignment_length"] == 21
    assert result["query_coverage"] == 1
    assert result["target_coverage"] == pytest.approx(9 / 21)
    assert result["distance"] is None
    assert "low_target_coverage" in result["distance_reason"]
    assert result["distance_status"] == "unresolved"


def test_reverse_order_preserves_distance_and_swaps_coverages():
    q, t = rec("q", "ACDEFGHIK"), rec("t", "ACDEFGHIKWWW")
    forward = similarity.align_pair(q, t, settings(minimum_coverage=0.5))
    reverse = similarity.align_pair(t, q, settings(minimum_coverage=0.5))
    assert forward["distance"] == reverse["distance"] == pytest.approx(0.25)
    assert forward["query_coverage"] == reverse["target_coverage"]
    assert forward["target_coverage"] == reverse["query_coverage"]


def test_unrelated_sequences_are_not_assigned_maximum_diversity():
    result = similarity.align_pair(rec("a", "AAAAAAAA"), rec("w", "WWWWWWWW"), settings())
    assert result["identity"] == 0
    assert result["distance_status"] == "unresolved"
    assert result["distance"] is None
    assert "low_identity" in result["distance_reason"]


def test_internal_gap_counts_only_paired_residues_for_coverage():
    result = similarity.align_pair(rec("q", "ACDEFGHIK"), rec("t", "ACDEWFGHIK"), settings())
    assert result["alignment_length"] == 10
    assert result["aligned_pairs"] == result["identical_residues"] == 9
    assert result["identity"] == result["query_coverage"] == 1
    assert result["target_coverage"] == 0.9
    assert result["distance"] == pytest.approx(0.1)


def test_coverage_boundary_is_inclusive_and_score_gate_is_strict():
    q, t = rec("q", "ACDEFGHIK"), rec("t", "ACDEFGHIKWWW")
    assert similarity.align_pair(q, t, settings(minimum_coverage=0.75))["distance_status"] == "resolved"
    assert similarity.align_pair(q, t, settings(minimum_coverage=0.751))["distance_status"] == "unresolved"
    score = similarity.align_pair(q, q, settings())["alignment_score"]
    assert similarity.align_pair(q, q, settings(minimum_score=score))["distance"] is None


def test_short_alignment_is_unresolved_even_when_identical():
    result = similarity.align_pair(rec("a", "ACDE"), rec("b", "ACDE"))
    assert result["distance"] is None
    assert result["distance_reason"] == "too_few_aligned_residues"


def test_matrix_order_symmetry_nan_and_singleton():
    rows = [rec("z", "ACDEFGHIK"), rec("a", "ACDEYGHIK"), rec("w", "WWWWWWWWWWWWWWWWWWWWW")]
    ids, matrix, pairs = similarity.build_distance_matrix(rows, settings())
    assert ids == ["z", "a", "w"]
    assert len(pairs) == 3
    np.testing.assert_allclose(matrix, matrix.T, equal_nan=True)
    np.testing.assert_array_equal(np.diag(matrix), [0, 0, 0])
    assert matrix[0, 1] == pytest.approx(1 / 9)
    assert np.isnan(matrix[0, 2])
    _, single, pairs = similarity.build_distance_matrix(rows[:1], settings())
    assert single.tolist() == [[0.0]] and pairs == []


@pytest.mark.parametrize("kwargs", [
    {"minimum_coverage": -0.1}, {"minimum_identity": 1.1},
    {"minimum_score": float("nan")}, {"gap_open_score": 1},
    {"gap_extend_score": float("inf")}, {"minimum_coverage": True},
    {"maximum_pairs": 0}, {"maximum_sequence_length": 1.5},
])
def test_invalid_settings_rejected(kwargs):
    with pytest.raises(ValueError):
        settings(**kwargs)


@pytest.mark.parametrize("rows", [[], [rec("a", "")], [rec("a", "ACDX")],
    [rec("same", "ACDE"), rec("same", "ACDF")]])
def test_invalid_fasta_records_rejected(rows):
    with pytest.raises(ValueError):
        similarity.build_distance_matrix(rows, settings())


def test_resource_limits_fail_before_alignment():
    with pytest.raises(ValueError, match="length"):
        similarity.build_distance_matrix([rec("a", "ACDE")], settings(maximum_sequence_length=3))
    with pytest.raises(ValueError, match="pairs"):
        similarity.build_distance_matrix([rec(str(i), "ACDE") for i in range(3)], settings(maximum_pairs=1))


def test_real_cli_outputs_repeatable_auditable_artifacts(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">z\nACDEFGHIK\n>a\nACDEYGHIK\n>w\nWWWWWWWWWWWWWWWWWWWWW\n", encoding="utf-8")
    input_hash = sha256_file(fasta)
    config_file = tmp_path / "config.yaml"
    hashes = []
    for run in ("one", "two"):
        output = tmp_path / run
        config_file.write_text(yaml.safe_dump({
            "input_fasta": str(fasta), "output_dir": str(output),
            "pool_label": "synthetic_test", "alignment": {"minimum_aligned_residues": 1},
        }), encoding="utf-8")
        script = Path(__file__).resolve().parents[1] / "experiments/run_similarity.py"
        process = subprocess.run([sys.executable, str(script), "--config", str(config_file)],
                                 capture_output=True, text=True, timeout=30)
        assert process.returncode == 0, process.stderr
        summary = json.loads(process.stdout)
        assert summary["pair_count"] == 3
        assert summary["resolved_pairs"] == 1
        assert summary["unresolved_pairs"] == 2
        assert json.loads((output / "candidate_distance_ids.json").read_text()) == ["z", "a", "w"]
        matrix = np.load(output / "candidate_distance.npy", allow_pickle=False)
        assert np.isnan(matrix[0, 2])
        with (output / "pairwise_similarity.tsv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        assert rows[1]["distance"] == ""
        assert rows[1]["distance_status"] == "unresolved"
        manifest = json.loads((output / "manifest.json").read_text())
        assert manifest["inputs"]["fasta"]["sha256"] == input_hash
        for name, digest in manifest["outputs"].items():
            assert sha256_file(output / name) == digest
        hashes.append(manifest["outputs"])
        again = subprocess.run([sys.executable, str(script), "--config", str(config_file)],
                               capture_output=True, text=True, timeout=30)
        assert again.returncode != 0  # Never silently overwrite a previous run.
    assert hashes[0] == hashes[1]
    assert sha256_file(fasta) == input_hash


def test_committed_development_artifacts_match_manifest_and_matrix():
    folder = Path(__file__).resolve().parents[1] / "results/gv02_03_similarity_dev"
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    for name, digest in manifest["outputs"].items():
        assert sha256_file(folder / name) == digest
    summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
    matrix = np.load(folder / "candidate_distance.npy", allow_pickle=False)
    ids = json.loads((folder / "candidate_distance_ids.json").read_text(encoding="utf-8"))
    assert matrix.shape == (len(ids), len(ids))
    np.testing.assert_allclose(matrix, matrix.T, equal_nan=True)
    np.testing.assert_array_equal(np.diag(matrix), np.zeros(len(ids)))
    upper = matrix[np.triu_indices(len(ids), 1)]
    assert len(upper) == summary["pair_count"]
    assert np.isfinite(upper).sum() == summary["resolved_pairs"]
    assert np.isnan(upper).sum() == summary["unresolved_pairs"]
    check = json.loads((folder / "reproducibility_check.json").read_text(encoding="utf-8"))
    assert check["run_1"] == check["run_2"] == manifest["outputs"]
