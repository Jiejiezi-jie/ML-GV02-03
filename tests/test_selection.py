from __future__ import annotations

import numpy as np
import pandas as pd

from gv_eval.selection import jaccard, non_dominated_sort, pareto_ranking, round_robin_top_k, weighted_ranking


def _frame():
    return pd.DataFrame(
        {
            "sequence_id": ["a", "b", "c", "d"],
            "m1": [1.0, 0.8, 0.2, 0.1],
            "m2": [0.1, 0.4, 0.9, 0.0],
        }
    )


def test_weighted_ranking_is_deterministic():
    one = weighted_ranking(_frame(), ["m1", "m2"], {"m1": 0.5, "m2": 0.5})
    two = weighted_ranking(_frame(), ["m1", "m2"], {"m1": 0.5, "m2": 0.5})
    assert one["sequence_id"].tolist() == two["sequence_id"].tolist()
    assert one.iloc[0]["sequence_id"] == "b"


def test_pareto_front_and_round_robin_budget():
    frame = _frame()
    ranks = non_dominated_sort(frame[["m1", "m2"]].to_numpy())
    assert ranks.tolist() == [0, 0, 0, 1]
    pareto = pareto_ranking(frame, ["m1", "m2"])
    assert set(pareto.head(3)["sequence_id"]) == {"a", "b", "c"}
    selected = round_robin_top_k(frame, ["m1", "m2"], 3)
    assert selected["sequence_id"].tolist() == ["a", "c", "b"]
    assert len(selected) == 3


def test_jaccard():
    assert jaccard(["a", "b"], ["b", "c"]) == 1 / 3
    assert jaccard([], []) == 1.0

