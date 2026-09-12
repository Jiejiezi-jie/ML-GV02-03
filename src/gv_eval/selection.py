from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


def validate_metrics(frame: pd.DataFrame, metrics: Iterable[str]) -> list[str]:
    metrics = list(metrics)
    missing = [column for column in metrics if column not in frame]
    if missing:
        raise ValueError(f"Missing selection metrics: {missing}")
    values = frame[metrics].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Selection metrics contain non-finite values")
    if ((values < 0) | (values > 1)).any():
        raise ValueError("Selection metrics must be in [0, 1]")
    return metrics


def weighted_ranking(
    frame: pd.DataFrame, metrics: Iterable[str], weights: Mapping[str, float]
) -> pd.DataFrame:
    metrics = validate_metrics(frame, metrics)
    weight = np.array([float(weights[name]) for name in metrics], dtype=float)
    if np.any(weight < 0) or not np.isfinite(weight).all() or weight.sum() <= 0:
        raise ValueError("Weights must be finite, nonnegative, and have positive sum")
    weight /= weight.sum()
    result = frame.copy()
    result["aggregate_score"] = result[metrics].to_numpy(dtype=float) @ weight
    return result.sort_values(
        ["aggregate_score", "sequence_id"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)


def non_dominated_sort(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    count = values.shape[0]
    dominates: list[list[int]] = [[] for _ in range(count)]
    dominated_by = np.zeros(count, dtype=int)
    ranks = np.full(count, -1, dtype=int)
    fronts: list[list[int]] = [[]]
    for left in range(count):
        for right in range(left + 1, count):
            left_dominates = np.all(values[left] >= values[right]) and np.any(
                values[left] > values[right]
            )
            right_dominates = np.all(values[right] >= values[left]) and np.any(
                values[right] > values[left]
            )
            if left_dominates:
                dominates[left].append(right)
                dominated_by[right] += 1
            elif right_dominates:
                dominates[right].append(left)
                dominated_by[left] += 1
    fronts[0] = [index for index in range(count) if dominated_by[index] == 0]
    for index in fronts[0]:
        ranks[index] = 0
    front_index = 0
    while fronts[front_index]:
        next_front: list[int] = []
        for left in fronts[front_index]:
            for right in dominates[left]:
                dominated_by[right] -= 1
                if dominated_by[right] == 0:
                    ranks[right] = front_index + 1
                    next_front.append(right)
        front_index += 1
        fronts.append(next_front)
    return ranks


def crowding_distance(values: np.ndarray, ranks: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    crowding = np.zeros(values.shape[0], dtype=float)
    for rank in sorted(set(ranks.tolist())):
        members = np.where(ranks == rank)[0]
        if len(members) <= 2:
            crowding[members] = np.inf
            continue
        for dimension in range(values.shape[1]):
            ordered = members[np.argsort(values[members, dimension], kind="mergesort")]
            crowding[ordered[0]] = crowding[ordered[-1]] = np.inf
            span = values[ordered[-1], dimension] - values[ordered[0], dimension]
            if span == 0:
                continue
            for position in range(1, len(ordered) - 1):
                current = ordered[position]
                if np.isfinite(crowding[current]):
                    crowding[current] += (
                        values[ordered[position + 1], dimension]
                        - values[ordered[position - 1], dimension]
                    ) / span
    return crowding


def pareto_ranking(frame: pd.DataFrame, metrics: Iterable[str]) -> pd.DataFrame:
    metrics = validate_metrics(frame, metrics)
    values = frame[metrics].to_numpy(dtype=float)
    ranks = non_dominated_sort(values)
    crowding = crowding_distance(values, ranks)
    result = frame.copy()
    result["pareto_rank"] = ranks
    result["crowding_distance"] = crowding
    return result.sort_values(
        ["pareto_rank", "crowding_distance", "sequence_id"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def round_robin_top_k(
    frame: pd.DataFrame, metrics: Iterable[str], k: int
) -> pd.DataFrame:
    metrics = validate_metrics(frame, metrics)
    rankings = {
        metric: frame.sort_values(
            [metric, "sequence_id"], ascending=[False, True], kind="mergesort"
        )["sequence_id"].tolist()
        for metric in metrics
    }
    lookup = frame.set_index("sequence_id", drop=False)
    positions = defaultdict(int)
    selected: list[tuple[str, str, int]] = []
    seen: set[str] = set()
    while len(selected) < min(k, len(frame)):
        progress = False
        for metric in metrics:
            ranking = rankings[metric]
            while positions[metric] < len(ranking) and ranking[positions[metric]] in seen:
                positions[metric] += 1
            if positions[metric] >= len(ranking):
                continue
            identifier = ranking[positions[metric]]
            positions[metric] += 1
            seen.add(identifier)
            selected.append((identifier, metric, positions[metric]))
            progress = True
            if len(selected) == min(k, len(frame)):
                break
        if not progress:
            break
    result = lookup.loc[[item[0] for item in selected]].copy()
    result["selected_by_metric"] = [item[1] for item in selected]
    result["source_metric_rank"] = [item[2] for item in selected]
    return result.reset_index(drop=True)


def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    a, b = set(left), set(right)
    union = a | b
    return len(a & b) / len(union) if union else 1.0

