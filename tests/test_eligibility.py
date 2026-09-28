from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gv_eval.analysis import (
    build_eligible_strategy_results,
    eligibility_summary,
    prepare_eligible_analysis_pool,
)
from gv_eval.selection import (
    filter_eligible_candidates,
    round_robin_top_k,
    validate_selection_budget,
)


def _candidate_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sequence_id": ["a", "b", "c", "d", "e", "f"],
            "qc_pass": [True, False, True, True, True, True],
            "family_status": [
                "supported_gvpa",
                "supported_gvpa",
                "other_gvp",
                "gvpa_gvpj_ambiguous",
                "supported_gvpa",
                "supported_gvpa",
            ],
            "domain_pass": [True, True, True, True, False, True],
            "m1": [0.9, 1.0, 1.0, 1.0, 1.0, 0.7],
            "m2": [0.6, 1.0, 1.0, 1.0, 1.0, 0.8],
        }
    )


def _distance() -> np.ndarray:
    positions = np.arange(6, dtype=float)
    return np.abs(positions[:, None] - positions[None, :]) / 10.0


def test_filter_eligible_candidates_enforces_all_three_gates():
    eligible = filter_eligible_candidates(_candidate_frame())
    assert eligible["sequence_id"].tolist() == ["a", "f"]
    assert eligible["qc_pass"].all()
    assert eligible["domain_pass"].all()
    assert eligible["family_status"].eq("supported_gvpa").all()


def test_filter_rejects_missing_or_invalid_evidence():
    frame = _candidate_frame().drop(columns="family_status")
    with pytest.raises(ValueError, match="Missing eligibility columns"):
        filter_eligible_candidates(frame)

    frame = _candidate_frame()
    frame["qc_pass"] = frame["qc_pass"].astype(object)
    frame.loc[0, "qc_pass"] = "yes"
    with pytest.raises(ValueError, match="must contain only booleans"):
        filter_eligible_candidates(frame)

    frame = _candidate_frame()
    frame.loc[0, "family_status"] = " "
    with pytest.raises(ValueError, match="family_status requires nonempty strings"):
        filter_eligible_candidates(frame)

    frame = _candidate_frame()
    frame.loc[0, "sequence_id"] = " "
    with pytest.raises(ValueError, match="nonempty string sequence identifiers"):
        filter_eligible_candidates(frame)


def test_filter_rejects_duplicate_identifiers():
    frame = _candidate_frame()
    frame.loc[1, "sequence_id"] = "a"
    with pytest.raises(ValueError, match="duplicate sequence identifiers"):
        filter_eligible_candidates(frame)


def test_eligibility_summary_is_a_sequential_audit_funnel():
    summary = eligibility_summary(_candidate_frame())
    assert summary == {
        "required_family_status": "supported_gvpa",
        "input_candidate_count": 6,
        "qc_passing_count": 5,
        "supported_gvpa_count": 3,
        "domain_passing_count": 2,
        "eligible_candidate_count": 2,
        "excluded_qc_count": 1,
        "excluded_family_count": 2,
        "excluded_domain_count": 1,
    }
    assert (
        summary["excluded_qc_count"]
        + summary["excluded_family_count"]
        + summary["excluded_domain_count"]
        + summary["eligible_candidate_count"]
        == summary["input_candidate_count"]
    )


def test_prepare_pool_slices_distance_matrix_in_candidate_order():
    eligible, distance, summary = prepare_eligible_analysis_pool(
        _candidate_frame(),
        _distance(),
        distance_ids=_candidate_frame()["sequence_id"],
    )
    assert eligible["sequence_id"].tolist() == ["a", "f"]
    np.testing.assert_allclose(distance, [[0.0, 0.5], [0.5, 0.0]])
    assert summary["eligible_candidate_count"] == 2


def test_prepare_pool_rejects_misaligned_distance_matrix():
    with pytest.raises(ValueError, match="does not match"):
        prepare_eligible_analysis_pool(
            _candidate_frame(),
            np.zeros((5, 5)),
            distance_ids=_candidate_frame()["sequence_id"],
        )


def test_prepare_pool_aligns_distance_matrix_by_sequence_id():
    frame = _candidate_frame()
    order = [5, 0, 4, 3, 2, 1]
    shuffled_ids = frame.iloc[order]["sequence_id"].tolist()
    shuffled_distance = _distance()[np.ix_(order, order)]

    eligible, distance, _ = prepare_eligible_analysis_pool(
        frame,
        shuffled_distance,
        distance_ids=shuffled_ids,
    )

    assert eligible["sequence_id"].tolist() == ["a", "f"]
    np.testing.assert_allclose(distance, [[0.0, 0.5], [0.5, 0.0]])


def test_prepare_pool_rejects_mismatched_distance_ids():
    ids = _candidate_frame()["sequence_id"].tolist()
    ids[-1] = "unexpected"
    with pytest.raises(ValueError, match="IDs do not match candidate IDs"):
        prepare_eligible_analysis_pool(
            _candidate_frame(),
            _distance(),
            distance_ids=ids,
        )


@pytest.mark.parametrize(
    ("distance", "message"),
    [
        (np.full((6, 6), np.nan), "non-finite"),
        (np.eye(6) * 2.0, "must be in"),
        (np.triu(np.ones((6, 6)), k=1) / 10.0, "symmetric"),
        (np.eye(6) / 10.0, "diagonal must be zero"),
    ],
)
def test_prepare_pool_rejects_invalid_distance_values(distance, message):
    with pytest.raises(ValueError, match=message):
        prepare_eligible_analysis_pool(
            _candidate_frame(),
            distance,
            distance_ids=_candidate_frame()["sequence_id"],
        )


def test_all_strategies_only_rank_eligible_candidates():
    selections, summaries, overlap, audit = build_eligible_strategy_results(
        _candidate_frame(),
        ["m1", "m2"],
        {"m1": 0.5, "m2": 0.5},
        [2],
        _distance(),
        distance_ids=_candidate_frame()["sequence_id"],
    )
    assert audit["eligible_candidate_count"] == 2
    assert len(summaries) == 3
    assert len(overlap) == 6
    for selected in selections.values():
        assert set(selected["sequence_id"]) == {"a", "f"}
        assert selected["selection_rank"].tolist() == [1, 2]
        assert selected["qc_pass"].all()
        assert selected["domain_pass"].all()
        assert selected["family_status"].eq("supported_gvpa").all()


def test_infeasible_budgets_raise_instead_of_silently_filling():
    with pytest.raises(ValueError, match="only 2 available"):
        validate_selection_budget(2, 3)
    with pytest.raises(ValueError, match="only 2 available"):
        build_eligible_strategy_results(
            _candidate_frame(),
            ["m1", "m2"],
            {"m1": 0.5, "m2": 0.5},
            [3],
            _distance(),
            distance_ids=_candidate_frame()["sequence_id"],
        )
    with pytest.raises(ValueError, match="only 2 available"):
        round_robin_top_k(_candidate_frame().head(2), ["m1", "m2"], 3)
