"""Descriptive sequence analysis of a fixed, four-objective Pareto front.

The feature list is fixed in advance.  No significance tests or functional
claims are made: the front is selected using scores from this same pool.
"""
from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Iterable, Mapping
from itertools import groupby
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from .selection import non_dominated_sort


METRICS = (
    "constraint_score", "conservation_score", "novelty_score", "uniqueness_score"
)
FEATURES = (
    "length", "shannon_entropy_bits", "max_residue_fraction", "charged_fraction",
    "gly_pro_fraction", "longest_homopolymer",
)
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_COLUMNS = tuple(f"aa_{aa}" for aa in AMINO_ACIDS)
_FLAG_COLUMNS = ("manual_review", "c_manual_review")
_COVERAGE_COLUMN = "uniqueness_resolved_fraction"


def _validated_sequences(
    sequences: Mapping[str, str] | Iterable[tuple[str, str]],
) -> dict[str, str]:
    records = sequences.items() if isinstance(sequences, Mapping) else sequences
    result: dict[str, str] = {}
    seen_sequences: dict[str, str] = {}
    for identifier, sequence in records:
        if not isinstance(identifier, str) or not identifier or any(
            character.isspace() for character in identifier
        ):
            raise ValueError("Sequence identifiers must be nonempty strings without whitespace")
        if identifier in result:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError(f"Empty or non-string sequence: {identifier}")
        illegal = sorted(set(sequence) - set(AMINO_ACIDS))
        if illegal:
            raise ValueError(f"Nonstandard amino acids in {identifier}: {illegal}")
        if sequence in seen_sequences:
            raise ValueError(
                f"Duplicate FASTA sequence: {identifier} and {seen_sequences[sequence]}"
            )
        result[identifier] = sequence
        seen_sequences[sequence] = identifier
    if not result:
        raise ValueError("FASTA must contain at least one sequence")
    return result


def read_fasta_strict(path: str | Path) -> dict[str, str]:
    """Read standard uppercase protein FASTA and reject ambiguous input.

    The first whitespace-delimited header token is the ID.  Wrapped sequences
    are supported. Duplicate IDs and identical sequences both raise errors.
    """
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    fragments: list[str] = []
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, "".join(fragments)))
            tokens = line[1:].split()
            if not tokens:
                raise ValueError(f"Empty FASTA header at line {line_number}")
            current_id, fragments = tokens[0], []
        elif current_id is None:
            raise ValueError(f"FASTA sequence before first header at line {line_number}")
        else:
            fragments.append(line)
    if current_id is not None:
        records.append((current_id, "".join(fragments)))
    return _validated_sequences(records)


def sequence_features(sequence: str) -> dict[str, float | int]:
    """Return six prespecified features and all 20 residue frequencies."""
    if not isinstance(sequence, str) or not sequence or set(sequence) - set(AMINO_ACIDS):
        raise ValueError("Features require a nonempty, standard uppercase amino acid sequence")
    counts = Counter(sequence)
    length = len(sequence)
    probabilities = np.array(list(counts.values()), dtype=float) / length
    return {
        "length": length,
        "shannon_entropy_bits": float(-np.sum(probabilities * np.log2(probabilities))),
        "max_residue_fraction": max(counts.values()) / length,
        "charged_fraction": sum(counts[aa] for aa in "DEKR") / length,
        "gly_pro_fraction": sum(counts[aa] for aa in "GP") / length,
        "longest_homopolymer": max(sum(1 for _ in run) for _, run in groupby(sequence)),
        **{f"aa_{aa}": counts[aa] / length for aa in AMINO_ACIDS},
    }


def cliffs_delta(front: Iterable[float], control: Iterable[float]) -> float:
    """P(front > control) - P(front < control); ties contribute zero."""
    left = np.asarray(list(front), dtype=float)
    right = np.asarray(list(control), dtype=float)
    if not left.size or not right.size or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("Cliff's delta requires two nonempty finite groups")
    difference = left[:, None] - right[None, :]
    return float((np.count_nonzero(difference > 0) - np.count_nonzero(difference < 0)) / difference.size)


def match_by_length(front: pd.DataFrame, control: pd.DataFrame) -> pd.DataFrame:
    """Match each front member to one unique control with minimum total gap.

    Both groups are sorted by ID before the one-to-one assignment, making the
    result independent of input row order. Multiple optimal solutions may
    exist; the selected one follows SciPy's deterministic assignment order.
    """
    for name, frame in (("front", front), ("control", control)):
        if not {"sequence_id", "length"}.issubset(frame.columns):
            raise ValueError(f"{name} matching input requires sequence_id and length")
        if frame.empty or not frame["sequence_id"].is_unique:
            raise ValueError(f"{name} matching input must be nonempty with unique IDs")
        if not frame["sequence_id"].map(lambda value: isinstance(value, str) and bool(value)).all():
            raise ValueError(f"{name} matching input requires string IDs")
        lengths = frame["length"].to_numpy(dtype=float)
        if not np.isfinite(lengths).all() or (lengths <= 0).any() or not (lengths == np.floor(lengths)).all():
            raise ValueError(f"{name} matching input requires positive integer lengths")
    if set(front["sequence_id"]) & set(control["sequence_id"]):
        raise ValueError("Front and control IDs must be disjoint")
    if len(control) < len(front):
        raise ValueError("Insufficient nonfront controls for unique length matching")
    left = front.sort_values("sequence_id").reset_index(drop=True)
    right = control.sort_values("sequence_id").reset_index(drop=True)
    cost = np.abs(
        left["length"].to_numpy(dtype=float)[:, None]
        - right["length"].to_numpy(dtype=float)[None, :]
    )
    rows, columns = linear_sum_assignment(cost)
    return pd.DataFrame({
        "front_id": left.iloc[rows]["sequence_id"].to_numpy(),
        "control_id": right.iloc[columns]["sequence_id"].to_numpy(),
        "front_length": left.iloc[rows]["length"].to_numpy(dtype=int),
        "control_length": right.iloc[columns]["length"].to_numpy(dtype=int),
        "length_gap": cost[rows, columns].astype(int),
    })


def length_match_sensitivity(
    front: pd.DataFrame,
    control: pd.DataFrame,
    feature_columns: Iterable[str] = (*METRICS, *FEATURES, *AA_COLUMNS),
    *,
    runs: int = 100,
    seed: int = 42,
) -> dict[str, Any]:
    """Describe sensitivity to choosing among equally optimal length controls.

    This post-hoc check perturbs only the tie-breaking cost, never the integer
    total length gap: each edge costs gap * (front_count + 1) + U[0, 1).
    The total perturbation is less than front_count, so it cannot overcome
    one residue of total-gap difference. Every run checks that it achieves
    the original unperturbed optimum. All generated matches are retained.

    Runs are not bootstrap samples, independent experiments, or a uniform
    enumeration of optimal matches. Sign counts and ranges describe the
    generated alternatives only, without confidence intervals or p-values.
    """
    if isinstance(runs, bool) or not isinstance(runs, (int, np.integer)) or runs < 1:
        raise ValueError("Sensitivity runs must be a positive integer")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Sensitivity seed must be a nonnegative integer")
    baseline = match_by_length(front, control)
    optimum = int(baseline["length_gap"].sum())
    left = front.sort_values("sequence_id").reset_index(drop=True)
    right = control.sort_values("sequence_id").reset_index(drop=True)
    features = list(feature_columns)
    if not features or len(set(features)) != len(features):
        raise ValueError("Sensitivity feature names must be nonempty and unique")
    for name, frame in (("front", left), ("control", right)):
        missing = set(features) - set(frame.columns)
        if missing:
            raise ValueError(f"Missing {name} sensitivity features: {sorted(missing)}")
        if not np.isfinite(frame[features].to_numpy(dtype=float)).all():
            raise ValueError(f"{name} sensitivity features must be finite")
    left_ids = left["sequence_id"].to_numpy()
    right_ids = right["sequence_id"].to_numpy()
    left_lengths = left["length"].to_numpy(dtype=float)
    right_lengths = right["length"].to_numpy(dtype=float)
    gaps = np.abs(left_lengths[:, None] - right_lengths[None, :])
    left_features = left[features].to_numpy(dtype=float)
    right_features = right[features].to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    pair_rows: list[dict[str, Any]] = []
    difference_rows: list[dict[str, Any]] = []
    assignments: set[tuple[str, ...]] = set()
    control_sets: set[tuple[str, ...]] = set()
    for run in range(1, runs + 1):
        perturbed = gaps * (len(left) + 1) + rng.uniform(0.0, 1.0, size=gaps.shape)
        rows, columns = linear_sum_assignment(perturbed)
        total_gap = float(gaps[rows, columns].sum())
        if total_gap != optimum:
            raise RuntimeError(
                f"Tie perturbation changed primary length optimum in run {run}: "
                f"{total_gap} instead of {optimum}"
            )
        assignments.add(tuple(right_ids[columns]))
        control_sets.add(tuple(sorted(right_ids[columns])))
        for row, column in zip(rows, columns):
            pair_rows.append({
                "run": run, "front_id": left_ids[row], "control_id": right_ids[column],
                "front_length": int(left_lengths[row]), "control_length": int(right_lengths[column]),
                "length_gap": int(gaps[row, column]),
            })
        differences = (left_features[rows] - right_features[columns]).mean(axis=0)
        for feature, difference in zip(features, differences):
            difference_rows.append({"run": run, "feature": feature, "mean_difference": float(difference)})
    difference_frame = pd.DataFrame(difference_rows)
    zero_tolerance = 1e-12
    summary_rows: list[dict[str, Any]] = []
    for feature in features:
        values = difference_frame.loc[difference_frame["feature"].eq(feature), "mean_difference"]
        summary_rows.append({
            "feature": feature, "runs": int(runs),
            "min_mean_difference": float(values.min()),
            "median_mean_difference": float(values.median()),
            "max_mean_difference": float(values.max()),
            "negative_count": int((values < -zero_tolerance).sum()),
            "zero_count": int((values.abs() <= zero_tolerance).sum()),
            "positive_count": int((values > zero_tolerance).sum()),
        })
    return {
        "tie_match_pairs": pd.DataFrame(pair_rows),
        "tie_match_differences": difference_frame,
        "tie_match_summary": pd.DataFrame(summary_rows),
        "summary": {
            "runs": int(runs), "seed": int(seed), "baseline_total_gap": optimum,
            "unique_assignments": len(assignments), "unique_control_sets": len(control_sets),
            "sign_zero_tolerance": zero_tolerance,
            "method": "Edge cost = absolute_length_gap * (front_count + 1) + Uniform[0,1); each run must retain the original minimum total gap; ascending ID order before perturbation",
            "scope": "Post-hoc sensitivity to equivalent optimal control selection. All runs retained without selection by result. These are not bootstrap samples, independent experiments, significance tests, or uniform enumeration; ranges and sign counts apply only to sampled optimal assignments.",
        },
    }


def _stats(values: pd.Series) -> dict[str, float | int]:
    return {
        "n": int(len(values)), "mean": float(values.mean()),
        "median": float(values.median()), "q1": float(values.quantile(0.25)),
        "q3": float(values.quantile(0.75)), "min": float(values.min()),
        "max": float(values.max()),
    }


def _group_comparison(front: pd.DataFrame, control: pd.DataFrame, name: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for feature in (*METRICS, *FEATURES, *AA_COLUMNS):
        left, right = front[feature], control[feature]
        rows.append({
            "comparison": name, "feature": feature,
            **{f"front_{key}": value for key, value in _stats(left).items()},
            **{f"control_{key}": value for key, value in _stats(right).items()},
            "mean_difference": float(left.mean() - right.mean()),
            "cliffs_delta": cliffs_delta(left, right),
        })
    return rows


def _evidence_summary(frame: pd.DataFrame) -> dict[str, Any]:
    return {
        "n": int(len(frame)),
        "manual_review_counts": {
            column: int(frame[column].sum()) if column in frame else None
            for column in _FLAG_COLUMNS
        },
        "uniqueness_resolved_fraction": _stats(frame[_COVERAGE_COLUMN])
        if _COVERAGE_COLUMN in frame else None,
    }


def analyze_pareto_features(
    scores: pd.DataFrame,
    sequences: Mapping[str, str] | Iterable[tuple[str, str]],
) -> dict[str, Any]:
    """Analyze a frozen ranking pool without changing its scores or membership.

    Returns candidate_features, front_details, group_comparison, matched_pairs,
    matched_feature_differences, matched_difference_summary, layer_summary,
    tie_match_pairs, tie_match_differences, tie_match_summary DataFrames and
    a JSON-serializable summary dictionary.
    Objective percentiles are average ascending ranks divided by pool size,
    multiplied by 100. All strongest/weakest percentile ties are retained.
    Length matching reduces length imbalance only; it does not remove other
    confounders or establish functional effects.
    """
    sequence_map = _validated_sequences(sequences)
    required = {"sequence_id", *METRICS}
    if not required.issubset(scores.columns):
        raise ValueError(f"Missing required score columns: {sorted(required - set(scores.columns))}")
    if not scores.columns.is_unique:
        raise ValueError("Score column names must be unique")
    if scores.empty or not scores["sequence_id"].is_unique:
        raise ValueError("Scores must be nonempty with unique sequence IDs")
    if not scores["sequence_id"].map(
        lambda value: isinstance(value, str) and bool(value) and not any(c.isspace() for c in value)
    ).all():
        raise ValueError("Score identifiers must be nonempty strings without whitespace")
    score_ids = set(scores["sequence_id"])
    if score_ids != set(sequence_map):
        raise ValueError(
            f"Score/FASTA ID mismatch: missing={sorted(score_ids - set(sequence_map))}, "
            f"extra={sorted(set(sequence_map) - score_ids)}"
        )
    candidates = scores.sort_values("sequence_id").copy().reset_index(drop=True)
    for column in (*METRICS, _COVERAGE_COLUMN):
        if column not in candidates:
            continue
        try:
            values = pd.to_numeric(candidates[column], errors="raise").to_numpy(dtype=float)
        except (ValueError, TypeError) as error:
            raise ValueError(f"{column} requires numeric values") from error
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"{column} requires finite values in [0, 1]")
        candidates[column] = values
    for column in _FLAG_COLUMNS:
        if column in candidates and not candidates[column].map(
            lambda value: isinstance(value, (bool, np.bool_))
        ).all():
            raise ValueError(f"{column} must contain actual booleans without missing values")
    features = pd.DataFrame([
        sequence_features(sequence_map[identifier]) for identifier in candidates["sequence_id"]
    ])
    for column in ("length", "sequence_length", "global_query_length"):
        if column in candidates:
            try:
                lengths = pd.to_numeric(candidates[column], errors="raise").to_numpy(dtype=float)
            except (ValueError, TypeError) as error:
                raise ValueError(f"Invalid metadata length column: {column}") from error
            if not np.array_equal(lengths, features["length"].to_numpy()):
                raise ValueError(f"FASTA/metadata length mismatch: {column}")
    if "sequence_sha256" in candidates:
        actual_hashes = candidates["sequence_id"].map(
            lambda identifier: hashlib.sha256(sequence_map[identifier].encode("ascii")).hexdigest()
        )
        if not actual_hashes.equals(candidates["sequence_sha256"]):
            raise ValueError("FASTA/metadata sequence_sha256 mismatch")
    for column in features:
        candidates[column] = features[column]
    candidates["pareto_rank"] = non_dominated_sort(candidates[list(METRICS)].to_numpy())
    candidates["is_front"] = candidates["pareto_rank"].eq(0)
    percentile_columns = []
    for metric in METRICS:
        column = f"{metric}_percentile"
        percentile_columns.append(column)
        candidates[column] = candidates[metric].rank(method="average", ascending=True, pct=True) * 100
    percentiles = candidates[percentile_columns].to_numpy()
    candidates["strongest_metrics"] = [
        ";".join(metric for metric, value in zip(METRICS, row) if value == row.max())
        for row in percentiles
    ]
    candidates["weakest_metrics"] = [
        ";".join(metric for metric, value in zip(METRICS, row) if value == row.min())
        for row in percentiles
    ]
    front = candidates.loc[candidates["is_front"]].copy().reset_index(drop=True)
    rest = candidates.loc[~candidates["is_front"]].copy().reset_index(drop=True)
    if front.empty or rest.empty:
        raise ValueError("Pareto comparison requires both a nonempty front and nonfront controls")
    pairs = match_by_length(front, rest)
    sensitivity = length_match_sensitivity(front, rest)
    matched = candidates.set_index("sequence_id", drop=False).loc[pairs["control_id"]].reset_index(drop=True)
    comparisons = pd.DataFrame(
        _group_comparison(front, rest, "front_vs_rest")
        + _group_comparison(front, matched, "front_vs_length_matched")
    )
    paired_rows: list[dict[str, Any]] = []
    candidate_lookup = candidates.set_index("sequence_id")
    for pair in pairs.itertuples(index=False):
        for feature in (*METRICS, *FEATURES, *AA_COLUMNS):
            left_value = float(candidate_lookup.at[pair.front_id, feature])
            right_value = float(candidate_lookup.at[pair.control_id, feature])
            paired_rows.append({
                "front_id": pair.front_id, "control_id": pair.control_id,
                "feature": feature, "front_value": left_value,
                "control_value": right_value, "difference": left_value - right_value,
            })
    paired_differences = pd.DataFrame(paired_rows)
    paired_summary = pd.DataFrame([
        {"feature": feature, **_stats(paired_differences.loc[
            paired_differences["feature"].eq(feature), "difference"
        ])} for feature in (*METRICS, *FEATURES, *AA_COLUMNS)
    ])
    layer_rows: list[dict[str, Any]] = []
    for rank, members in candidates.groupby("pareto_rank", sort=True):
        for feature in (*METRICS, *FEATURES, *AA_COLUMNS):
            layer_rows.append({"pareto_rank": int(rank), "feature": feature, **_stats(members[feature])})
    trends: dict[str, float | None] = {}
    for feature in (*METRICS, *FEATURES, *AA_COLUMNS):
        trends[feature] = float(candidates[feature].corr(candidates["pareto_rank"], method="spearman")) \
            if candidates[feature].nunique() > 1 else None
    summary = {
        "candidate_count": int(len(candidates)), "front_count": int(len(front)),
        "nonfront_count": int(len(rest)), "pareto_layer_count": int(candidates["pareto_rank"].nunique()),
        "front_ids": front["sequence_id"].tolist(), "metrics": list(METRICS),
        "predefined_features": list(FEATURES), "composition_features": list(AA_COLUMNS),
        "percentile_definition": "100 * average ascending rank / pool size; all ties retained",
        "cliffs_delta_definition": "P(front > control) - P(front < control); ties contribute zero",
        "quantile_definition": "pandas linear interpolation at probabilities 0.25, 0.5 and 0.75",
        "layer_trend_spearman": trends,
        "layer_trend_note": "Descriptive pooled Spearman with zero-based Pareto layer; positive means increasing in worse layers. No p-values.",
        "evidence": {
            "front": _evidence_summary(front), "rest": _evidence_summary(rest),
            "length_matched": _evidence_summary(matched),
        },
        "length_matching": {
            "method": "Minimum total absolute length gap in one-to-one matching; scipy.optimize.linear_sum_assignment after ascending ID sort; no control reuse",
            "pairs": int(len(pairs)), "total_gap": int(pairs["length_gap"].sum()),
            "mean_gap": float(pairs["length_gap"].mean()), "max_gap": int(pairs["length_gap"].max()),
            "exact_match_count": int(pairs["length_gap"].eq(0).sum()),
            "front_length_range": [int(front["length"].min()), int(front["length"].max())],
            "available_control_length_range": [int(rest["length"].min()), int(rest["length"].max())],
            "front_outside_control_length_range": int(
                ((front["length"] < rest["length"].min()) | (front["length"] > rest["length"].max())).sum()
            ),
        },
        "length_match_sensitivity": sensitivity["summary"],
        "scope": [
            "Descriptive analysis of a frozen, selected pool; no p-values or causal/functional claims.",
            "The first front is defined by the same objectives being summarized; score differences are expected by selection.",
            "Length matching addresses length balance only; inspect saved gaps and coverage before interpreting feature differences.",
            "Amino acid frequencies sum to one; the 20 residue differences are not independent effects.",
            "Manual-review flags and unresolved-pair coverage remain part of the candidate evidence.",
            "GP fraction and DEKR fraction are composition summaries, not predictions of protein function or structure.",
        ],
    }
    return {
        "candidate_features": candidates, "front_details": front,
        "group_comparison": comparisons, "matched_pairs": pairs,
        "matched_feature_differences": paired_differences,
        "matched_difference_summary": paired_summary,
        "tie_match_pairs": sensitivity["tie_match_pairs"],
        "tie_match_differences": sensitivity["tie_match_differences"],
        "tie_match_summary": sensitivity["tie_match_summary"],
        "layer_summary": pd.DataFrame(layer_rows), "summary": summary,
    }
