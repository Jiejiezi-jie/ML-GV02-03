from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .metrics import subset_diversity
from .selection import (
    eligibility_mask,
    jaccard,
    pareto_ranking,
    round_robin_top_k,
    validate_selection_budget,
    weighted_ranking,
)


def correlation_tables(frame: pd.DataFrame, metrics: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    return frame[metrics].corr(method="pearson"), frame[metrics].corr(method="spearman")


def eligibility_summary(
    frame: pd.DataFrame,
    *,
    required_family_status: str = "supported_gvpa",
) -> dict[str, int | str]:
    """Build a mutually exclusive audit funnel for the formal GvpA pool."""

    eligible = eligibility_mask(frame, required_family_status=required_family_status)
    qc_pass = frame["qc_pass"]
    family_pass = qc_pass & frame["family_status"].eq(required_family_status)
    domain_pass = family_pass & frame["domain_pass"]
    return {
        "required_family_status": required_family_status,
        "input_candidate_count": len(frame),
        "qc_passing_count": int(qc_pass.sum()),
        "supported_gvpa_count": int(family_pass.sum()),
        "domain_passing_count": int(domain_pass.sum()),
        "eligible_candidate_count": int(eligible.sum()),
        "excluded_qc_count": int((~qc_pass).sum()),
        "excluded_family_count": int(
            (qc_pass & ~frame["family_status"].eq(required_family_status)).sum()
        ),
        "excluded_domain_count": int((family_pass & ~frame["domain_pass"]).sum()),
    }


def prepare_eligible_analysis_pool(
    frame: pd.DataFrame,
    distance: np.ndarray,
    *,
    distance_ids: Iterable[str],
    required_family_status: str = "supported_gvpa",
) -> tuple[pd.DataFrame, np.ndarray, dict[str, int | str]]:
    """Align by sequence ID, then filter candidates and their distance matrix."""

    distance = np.asarray(distance, dtype=float)
    distance_ids = list(distance_ids)
    expected_shape = (len(distance_ids), len(distance_ids))
    if distance.shape != expected_shape:
        raise ValueError(
            f"Candidate distance matrix shape {distance.shape} does not match {expected_shape}"
        )
    if not np.isfinite(distance).all():
        raise ValueError("Candidate distance matrix contains non-finite values")
    if ((distance < 0.0) | (distance > 1.0)).any():
        raise ValueError("Candidate distances must be in [0, 1]")
    if not np.allclose(distance, distance.T, rtol=0.0, atol=1e-12):
        raise ValueError("Candidate distance matrix must be symmetric")
    if not np.allclose(np.diag(distance), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Candidate distance matrix diagonal must be zero")
    if not all(isinstance(value, str) and value.strip() for value in distance_ids):
        raise ValueError("Distance matrix IDs must be nonempty strings")
    if len(set(distance_ids)) != len(distance_ids):
        raise ValueError("Distance matrix IDs contain duplicates")

    frame_ids = frame["sequence_id"].tolist()
    frame_id_set = set(frame_ids)
    distance_id_set = set(distance_ids)
    if frame_id_set != distance_id_set:
        missing = sorted(frame_id_set - distance_id_set)
        unexpected = sorted(distance_id_set - frame_id_set)
        raise ValueError(
            "Distance matrix IDs do not match candidate IDs: "
            f"missing={missing}, unexpected={unexpected}"
        )

    distance_position = {identifier: index for index, identifier in enumerate(distance_ids)}
    frame_order = [distance_position[identifier] for identifier in frame_ids]
    aligned_distance = distance[np.ix_(frame_order, frame_order)]
    mask = eligibility_mask(frame, required_family_status=required_family_status)
    positions = np.flatnonzero(mask.to_numpy())
    eligible_frame = frame.iloc[positions].copy().reset_index(drop=True)
    eligible_distance = aligned_distance[np.ix_(positions, positions)].copy()
    summary = eligibility_summary(frame, required_family_status=required_family_status)
    return eligible_frame, eligible_distance, summary


def selected_quality(
    selected: pd.DataFrame,
    metrics: list[str],
    distance: np.ndarray,
    position_by_id: dict[str, int],
) -> dict:
    positions = [position_by_id[value] for value in selected["sequence_id"]]
    result = {
        "selected_count": len(selected),
        "subset_diversity": subset_diversity(distance, positions),
        "domain_pass_rate": float(selected["domain_pass"].mean()),
    }
    result.update({f"mean_{metric}": float(selected[metric].mean()) for metric in metrics})
    return result


def build_strategy_results(
    frame: pd.DataFrame,
    metrics: list[str],
    equal_weights: dict[str, float],
    budgets: Iterable[int],
    distance: np.ndarray,
) -> tuple[dict[tuple[str, int], pd.DataFrame], list[dict], list[dict]]:
    weighted = weighted_ranking(frame, metrics, equal_weights)
    pareto = pareto_ranking(frame, metrics)
    position = {value: index for index, value in enumerate(frame["sequence_id"])}
    selections: dict[tuple[str, int], pd.DataFrame] = {}
    summaries: list[dict] = []
    overlap: list[dict] = []
    for budget in budgets:
        validate_selection_budget(len(frame), budget)
        variants = {
            "weighted_sum": weighted.head(budget).copy(),
            "pareto": pareto.head(budget).copy(),
            "dimension_round_robin": round_robin_top_k(frame, metrics, budget),
        }
        for name, selected in variants.items():
            selected.insert(0, "selection_rank", np.arange(1, len(selected) + 1))
            selections[(name, budget)] = selected
            summaries.append(
                {
                    "strategy": name,
                    "budget": budget,
                    **selected_quality(selected, metrics, distance, position),
                }
            )
        names = list(variants)
        for left_index, left in enumerate(names):
            for right in names[left_index:]:
                overlap.append(
                    {
                        "budget": budget,
                        "strategy_a": left,
                        "strategy_b": right,
                        "jaccard": jaccard(variants[left]["sequence_id"], variants[right]["sequence_id"]),
                        "intersection": len(
                            set(variants[left]["sequence_id"]) & set(variants[right]["sequence_id"])
                        ),
                    }
                )
    return selections, summaries, overlap


def build_eligible_strategy_results(
    frame: pd.DataFrame,
    metrics: list[str],
    equal_weights: dict[str, float],
    budgets: Iterable[int],
    distance: np.ndarray,
    *,
    distance_ids: Iterable[str],
    required_family_status: str = "supported_gvpa",
) -> tuple[
    dict[tuple[str, int], pd.DataFrame],
    list[dict],
    list[dict],
    dict[str, int | str],
]:
    """Apply strict V2 gates, then run all strategies on the aligned eligible pool."""

    eligible, eligible_distance, summary = prepare_eligible_analysis_pool(
        frame,
        distance,
        distance_ids=distance_ids,
        required_family_status=required_family_status,
    )
    selections, summaries, overlap = build_strategy_results(
        eligible,
        metrics,
        equal_weights,
        budgets,
        eligible_distance,
    )
    return selections, summaries, overlap, summary


def ablation_results(
    frame: pd.DataFrame,
    metrics: list[str],
    weights: dict[str, float],
    budget: int,
) -> tuple[list[dict], dict[str, list[str]]]:
    baseline = weighted_ranking(frame, metrics, weights).head(budget)
    baseline_ids = baseline["sequence_id"].tolist()
    rows: list[dict] = []
    selections = {"all_dimensions": baseline_ids}
    variants: list[tuple[str, list[str]]] = [("all_dimensions", metrics)]
    variants.extend((f"without_{metric}", [x for x in metrics if x != metric]) for metric in metrics)
    for name, retained in variants:
        ranking = weighted_ranking(frame, retained, {x: weights[x] for x in retained})
        selected = ranking.head(budget)
        selections[name] = selected["sequence_id"].tolist()
        row = {
            "variant": name,
            "removed_metric": "" if name == "all_dimensions" else name.removeprefix("without_"),
            "jaccard_with_full": jaccard(baseline_ids, selected["sequence_id"]),
        }
        row.update({f"mean_{metric}": float(selected[metric].mean()) for metric in metrics})
        rows.append(row)
    return rows, selections


def weight_robustness(
    frame: pd.DataFrame,
    metrics: list[str],
    weights: dict[str, float],
    budget: int,
    samples: int,
    relative_sigma: float,
    seed: int,
) -> tuple[list[dict], dict[str, float]]:
    baseline = weighted_ranking(frame, metrics, weights)
    baseline_ids = baseline.head(budget)["sequence_id"].tolist()
    base_position = {identifier: rank for rank, identifier in enumerate(baseline["sequence_id"])}
    frequency = {identifier: 0 for identifier in frame["sequence_id"]}
    base_weight = np.asarray([weights[name] for name in metrics], dtype=float)
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for sample in range(samples):
        perturbation = rng.uniform(1.0 - relative_sigma, 1.0 + relative_sigma, len(metrics))
        current = base_weight * perturbation
        current /= current.sum()
        ranking = weighted_ranking(frame, metrics, dict(zip(metrics, current)))
        chosen = ranking.head(budget)["sequence_id"].tolist()
        for identifier in chosen:
            frequency[identifier] += 1
        aligned_current = [0] * len(ranking)
        for rank, identifier in enumerate(ranking["sequence_id"]):
            aligned_current[base_position[identifier]] = rank
        rho = float(spearmanr(np.arange(len(frame)), aligned_current).statistic)
        rows.append(
            {
                "sample": sample + 1,
                **{f"weight_{name}": float(current[index]) for index, name in enumerate(metrics)},
                "top_k_jaccard": jaccard(baseline_ids, chosen),
                "rank_spearman": rho,
            }
        )
    return rows, {key: value / samples for key, value in frequency.items()}


def threshold_robustness(
    weighted_order: pd.DataFrame,
    budget: int,
    score_cutoff: float,
    coverage_cutoff: float,
    cutoff_multipliers: Iterable[float],
    coverage_offsets: Iterable[float],
) -> list[dict]:
    baseline_mask = (weighted_order["domain_score"] >= score_cutoff) & (
        weighted_order["model_coverage"] >= coverage_cutoff
    )
    baseline = weighted_order.loc[baseline_mask].head(budget)["sequence_id"].tolist()
    rows: list[dict] = []
    for multiplier in cutoff_multipliers:
        for offset in coverage_offsets:
            current_score = score_cutoff * multiplier
            current_coverage = float(np.clip(coverage_cutoff + offset, 0.0, 1.0))
            eligible = weighted_order.loc[
                (weighted_order["domain_score"] >= current_score)
                & (weighted_order["model_coverage"] >= current_coverage)
            ]
            chosen = eligible.head(budget)["sequence_id"].tolist()
            rows.append(
                {
                    "domain_cutoff_multiplier": multiplier,
                    "domain_score_cutoff": current_score,
                    "coverage_offset": offset,
                    "model_coverage_cutoff": current_coverage,
                    "eligible_count": len(eligible),
                    "selected_count": len(chosen),
                    "jaccard_with_baseline_gate": jaccard(baseline, chosen),
                }
            )
    return rows


def random_baseline(
    frame: pd.DataFrame,
    strategy_summary: list[dict],
    metrics: list[str],
    distance: np.ndarray,
    budget: int,
    samples: int,
    seed: int,
) -> tuple[list[dict], list[dict]]:
    position = {value: index for index, value in enumerate(frame["sequence_id"])}
    rng = np.random.default_rng(seed)
    random_rows: list[dict] = []
    for sample in range(samples):
        chosen_positions = sorted(rng.choice(len(frame), size=budget, replace=False).tolist())
        selected = frame.iloc[chosen_positions]
        row = {"sample": sample + 1}
        row.update(selected_quality(selected, metrics, distance, position))
        random_rows.append(row)
    random_frame = pd.DataFrame(random_rows)
    comparisons: list[dict] = []
    for strategy in strategy_summary:
        if strategy["budget"] != budget:
            continue
        for key in [*(f"mean_{metric}" for metric in metrics), "subset_diversity", "domain_pass_rate"]:
            observed = float(strategy[key])
            values = random_frame[key].to_numpy(float)
            comparisons.append(
                {
                    "strategy": strategy["strategy"],
                    "metric": key,
                    "observed": observed,
                    "random_mean": float(values.mean()),
                    "random_p95": float(np.quantile(values, 0.95)),
                    "empirical_p_random_ge_observed": float((np.sum(values >= observed) + 1) / (samples + 1)),
                }
            )
    return random_rows, comparisons


def save_figures(
    frame: pd.DataFrame,
    metrics: list[str],
    pearson: pd.DataFrame,
    spearman: pd.DataFrame,
    strategy_summary: pd.DataFrame,
    overlap: pd.DataFrame,
    weight_rows: pd.DataFrame,
    threshold_rows: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    import seaborn as sns

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for metric, axis in zip(metrics, axes.flat):
        sns.histplot(frame[metric], bins=20, kde=True, ax=axis)
        axis.set_title(metric)
    fig.tight_layout()
    fig.savefig(output_dir / "metric_distributions.png", dpi=180)
    plt.close(fig)

    pair = sns.pairplot(frame[metrics], corner=False, plot_kws={"s": 18, "alpha": 0.65})
    pair.fig.suptitle("Four-objective scatter matrix", y=1.02)
    pair.savefig(output_dir / "metric_scatter_matrix.png", dpi=180)
    plt.close(pair.fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.heatmap(pearson, annot=True, fmt=".2f", vmin=-1, vmax=1, cmap="vlag", ax=axes[0])
    axes[0].set_title("Pearson")
    sns.heatmap(spearman, annot=True, fmt=".2f", vmin=-1, vmax=1, cmap="vlag", ax=axes[1])
    axes[1].set_title("Spearman")
    fig.tight_layout()
    fig.savefig(output_dir / "correlation_heatmaps.png", dpi=180)
    plt.close(fig)

    primary = strategy_summary[strategy_summary["budget"] == 20].copy()
    long = primary.melt(
        id_vars="strategy",
        value_vars=[*(f"mean_{metric}" for metric in metrics), "subset_diversity"],
        var_name="metric",
        value_name="value",
    )
    fig, axis = plt.subplots(figsize=(12, 6))
    sns.barplot(data=long, x="metric", y="value", hue="strategy", ax=axis)
    axis.tick_params(axis="x", rotation=25)
    axis.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(output_dir / "strategy_quality_k20.png", dpi=180)
    plt.close(fig)

    primary_overlap = overlap[overlap["budget"] == 20].pivot(
        index="strategy_a", columns="strategy_b", values="jaccard"
    )
    labels = sorted(set(primary_overlap.index) | set(primary_overlap.columns))
    matrix = pd.DataFrame(np.eye(len(labels)), index=labels, columns=labels)
    for left in primary_overlap.index:
        for right in primary_overlap.columns:
            if not math.isnan(primary_overlap.loc[left, right]):
                matrix.loc[left, right] = matrix.loc[right, left] = primary_overlap.loc[left, right]
    fig, axis = plt.subplots(figsize=(6, 5))
    sns.heatmap(matrix, annot=True, fmt=".2f", vmin=0, vmax=1, cmap="Blues", ax=axis)
    axis.set_title("Top-20 selection Jaccard")
    fig.tight_layout()
    fig.savefig(output_dir / "strategy_overlap_k20.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    sns.histplot(weight_rows["top_k_jaccard"], bins=15, ax=axes[0])
    axes[0].set_title("Weight perturbation: Top-K Jaccard")
    sns.histplot(weight_rows["rank_spearman"], bins=15, ax=axes[1])
    axes[1].set_title("Weight perturbation: rank Spearman")
    fig.tight_layout()
    fig.savefig(output_dir / "weight_robustness.png", dpi=180)
    plt.close(fig)

    pivot = threshold_rows.pivot(
        index="coverage_offset", columns="domain_cutoff_multiplier", values="eligible_count"
    )
    fig, axis = plt.subplots(figsize=(7, 5))
    sns.heatmap(pivot, annot=True, fmt=".0f", cmap="YlOrRd", ax=axis)
    axis.set_title("Eligible candidates under domain thresholds")
    fig.tight_layout()
    fig.savefig(output_dir / "threshold_robustness.png", dpi=180)
    plt.close(fig)
