from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from gv_eval.io import sha256_file


ROOT = Path(__file__).resolve().parents[1]


def test_candidate_score_artifact_is_complete():
    frame = pd.read_csv(ROOT / "results/gv02_03/tables/candidate_scores.tsv", sep="\t")
    assert len(frame) == 200
    assert frame["sequence_id"].is_unique
    metrics = ["constraint_score", "conservation_score", "novelty_score", "uniqueness_score"]
    assert frame[metrics].notna().all().all()
    assert ((frame[metrics] >= 0) & (frame[metrics] <= 1)).all().all()


def test_distance_artifact_matches_candidate_order():
    matrix = np.load(ROOT / "data/processed/gv02_03/candidate_distance.npy", allow_pickle=False)
    identifiers = json.loads(
        (ROOT / "data/processed/gv02_03/candidate_distance_ids.json").read_text(encoding="utf-8")
    )
    assert matrix.shape == (len(identifiers), len(identifiers)) == (200, 200)
    np.testing.assert_allclose(matrix, matrix.T)
    np.testing.assert_allclose(np.diag(matrix), 0.0)


def test_all_strategy_budgets_are_exact_and_unique():
    for strategy in ("weighted_sum", "pareto", "dimension_round_robin"):
        for budget in (10, 20, 50):
            frame = pd.read_csv(
                ROOT / f"results/gv02_03/selections/{strategy}_top{budget}.tsv", sep="\t"
            )
            assert len(frame) == budget
            assert frame["sequence_id"].is_unique
            assert frame["selection_rank"].tolist() == list(range(1, budget + 1))


def test_manifest_input_and_output_hashes_match_files():
    manifest = json.loads((ROOT / "results/gv02_03/manifest.json").read_text(encoding="utf-8"))
    for section in ("inputs", "outputs"):
        for relative, expected in manifest[section].items():
            assert sha256_file(ROOT / relative) == expected


def test_recorded_reproducibility_hashes_are_equal_and_current():
    record = json.loads(
        (ROOT / "results/gv02_03/reproducibility_check.json").read_text(encoding="utf-8")
    )
    assert record["all_equal"] is True
    assert record["run_1"] == record["run_2"]
    for relative, expected in record["run_2"].items():
        assert sha256_file(ROOT / relative) == expected

