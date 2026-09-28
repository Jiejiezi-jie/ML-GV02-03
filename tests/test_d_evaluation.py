from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gv_eval.d_evaluation import distance_statistics, prepare_scored_pool, run_evaluation, validate_partial_distance
from gv_eval.io import read_fasta, sha256_file

ROOT = Path(__file__).resolve().parents[1]


def test_missing_distances_are_not_rewarded_as_maximum_diversity():
    stats = distance_statistics(np.array([0.2, 0.4, np.nan]))
    assert stats["resolved_fraction"] == pytest.approx(2 / 3)
    assert stats["observed_mean"] == pytest.approx(0.3)
    assert stats["lower_bound"] == pytest.approx(0.2)
    assert stats["upper_bound"] == pytest.approx(1.6 / 3)
    empty = distance_statistics(np.array([np.nan, np.nan]))
    assert empty["observed_mean"] is None
    assert empty["resolved_fraction"] == 0
    assert empty["lower_bound"] == 0
    assert empty["upper_bound"] == 1


@pytest.mark.parametrize("matrix", [
    np.array([[0, np.nan], [0.2, 0]]),
    np.array([[0, np.inf], [np.inf, 0]]),
    np.array([[0, 1.1], [1.1, 0]]),
    np.eye(2),
])
def test_partial_distance_rejects_invalid_evidence(matrix):
    with pytest.raises(ValueError):
        validate_partial_distance(matrix, ["a", "b"])


def test_hard_gates_before_uniqueness_and_id_alignment():
    frame = pd.DataFrame({
        "sequence_id": ["a", "b", "c", "other"], "qc_pass": [True] * 4,
        "family_status": ["supported_gvpa"] * 3 + ["other_gvp"],
        "domain_pass": [True] * 4, "global_distance_status": ["resolved"] * 4,
        "constraint_score": [0.8] * 4, "conservation_score": [0.8] * 4,
        "novelty_score": [0.2] * 4,
    })
    # Order differs from frame. c has no reliable partners and must not rank.
    distance = np.array([[0, np.nan, 0.2], [np.nan, 0, np.nan], [0.2, np.nan, 0]])
    ranking, aligned, pool = prepare_scored_pool(frame, distance, ["b", "c", "a"], 0.5)
    assert ranking.sequence_id.tolist() == ["a", "b"]
    assert ranking.uniqueness_score.tolist() == pytest.approx([0.2, 0.2])
    assert ranking.uniqueness_reference_pool_size.tolist() == [3, 3]
    np.testing.assert_allclose(aligned, [[0, 0.2], [0.2, 0]])
    assert not pool.loc[pool.sequence_id.eq("c"), "score_ready"].item()


def test_frozen_handoff_end_to_end_and_reproducibility(tmp_path):
    config_path = ROOT / "configs/d_evaluation.json"
    first, second = tmp_path / "first", tmp_path / "second"
    summary = run_evaluation(ROOT, config_path, first)
    assert summary["retrained"] is False
    assert summary["funnel"]["input_candidate_count"] == 1000
    assert summary["funnel"]["eligible_candidate_count"] == 136
    assert summary["ranking_count"] == 114
    ranking = pd.read_csv(first / "ranking_pool_scores.tsv", sep="\t")
    for row in summary["strategy_summary"]:
        stem = f'{row["strategy"]}_top{row["budget"]}'
        selected = pd.read_csv(first / f"{stem}.tsv", sep="\t")
        assert len(selected) == row["budget"]
        assert selected.sequence_id.is_unique
        assert set(selected.sequence_id) <= set(ranking.sequence_id)
        assert selected.family_status.eq("supported_gvpa").all()
        assert selected.domain_pass.all() and selected.qc_pass.all()
        assert [r.identifier for r in read_fasta(first / f"{stem}.fasta")] == selected.sequence_id.tolist()
    audit = pd.read_csv(first / "candidate_audit.tsv", sep="\t")
    assert len(audit) == 1000
    assert audit.ranking_eligible.sum() == 114
    manifest = json.loads((first / "manifest.json").read_text())
    for name, digest in manifest["output_sha256"].items():
        assert sha256_file(first / name) == digest
    run_evaluation(ROOT, config_path, second)
    repeated = json.loads((second / "manifest.json").read_text())
    assert manifest == repeated
    with pytest.raises(ValueError, match="new or empty"):
        run_evaluation(ROOT, config_path, first)


def test_insufficient_budget_is_audited_without_padding(tmp_path):
    config = json.loads((ROOT / "configs/d_evaluation.json").read_text())
    config.update(budgets=[200], primary_k=200, weight_samples=1, random_samples=1)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    destination = tmp_path / "results"
    summary = run_evaluation(ROOT, config_path, destination)
    assert summary["primary_budget_status"] == "insufficient"
    assert summary["strategy_summary"] == []
    budgets = pd.read_csv(destination / "budget_audit.tsv", sep="\t")
    assert budgets.status.tolist() == ["insufficient"]
    assert not list(destination.glob("*_top200.fasta"))
