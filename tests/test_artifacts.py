from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from gv_eval.io import sha256_file


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/gv02_03_v2/d_evaluation_member_a_v1"


def test_candidate_score_artifact_is_complete():
    frame = pd.read_csv(RESULTS / "ranking_pool_scores.tsv", sep="\t")
    assert len(frame) == 114
    assert frame["sequence_id"].is_unique
    metrics = ["constraint_score", "conservation_score", "novelty_score", "uniqueness_score"]
    assert frame[metrics].notna().all().all()
    assert ((frame[metrics] >= 0) & (frame[metrics] <= 1)).all().all()


def test_distance_artifact_matches_candidate_order():
    matrix = np.load(RESULTS / "ranking_distance.npy", allow_pickle=False)
    identifiers = json.loads(
        (RESULTS / "ranking_distance_ids.json").read_text(encoding="utf-8")
    )
    assert matrix.shape == (len(identifiers), len(identifiers)) == (114, 114)
    np.testing.assert_allclose(matrix, matrix.T, equal_nan=True)
    np.testing.assert_allclose(np.diag(matrix), 0.0)


def test_all_strategy_budgets_are_exact_and_unique():
    for strategy in ("weighted_sum", "pareto", "dimension_round_robin"):
        for budget in (10, 20, 50):
            frame = pd.read_csv(
                RESULTS / f"{strategy}_top{budget}.tsv", sep="\t"
            )
            assert len(frame) == budget
            assert frame["sequence_id"].is_unique
            assert frame["selection_rank"].tolist() == list(range(1, budget + 1))


def test_manifest_input_and_output_hashes_match_files():
    from gv_eval.d_evaluation import FrozenInputs
    manifest = json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8"))
    inputs = FrozenInputs(ROOT, manifest["source_commits"])
    for relative, expected in manifest["input_sha256"].items():
        owner, path = relative.split(":", 1)
        inputs.read(owner, path, expected)
    for relative, expected in manifest["output_sha256"].items():
        assert sha256_file(RESULTS / relative) == expected


def test_recorded_reproducibility_hashes_are_equal_and_current():
    record = json.loads(
        (RESULTS.parent / "d_evaluation_member_a_v1_reproducibility.json").read_text(encoding="utf-8")
    )
    assert record["passed"] is True
    assert record["compared_artifact_count"] == 44
    assert record["manifest_sha256"] == sha256_file(RESULTS / "manifest.json")

