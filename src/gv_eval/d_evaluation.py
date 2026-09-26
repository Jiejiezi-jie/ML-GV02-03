"""Evaluate frozen B/C deliveries without training, sampling or external tools."""
from __future__ import annotations

import hashlib
import io
import json
import subprocess
import platform
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .analysis import ablation_results, eligibility_summary, weight_robustness
from .io import read_fasta, sequence_sha256, sha256_file, write_fasta, write_json
from .selection import eligibility_mask, jaccard, pareto_ranking, round_robin_top_k, weighted_ranking

METRICS = ["constraint_score", "conservation_score", "novelty_score", "uniqueness_score"]
WEIGHTS = dict.fromkeys(METRICS, 0.25)
B_POOL = "data/processed/gv02_03_v2/vae_member_a_v1"
B_SCORE = "results/gv02_03_v2/vae_member_a_v1/scoring"
C_BASE = "results/gv02_03_c_handoff_member_a_v1"
C_FINAL = "results/gv02_03_c_final_handoff_member_a_v1"


class FrozenInputs:
    """Read exact committed bytes; never depend on the currently checked-out branch."""

    def __init__(self, root: Path, revisions: dict[str, str]):
        self.root = root
        self.revisions = {
            name: subprocess.check_output(
                ["git", "rev-parse", "--verify", f"{revision}^{{commit}}"], cwd=root,
                text=True,
            ).strip()
            for name, revision in revisions.items()
        }
        self.hashes: dict[str, str] = {}

    def read(self, owner: str, path: str, expected: str | None = None) -> bytes:
        data = subprocess.check_output(
            ["git", "show", f"{self.revisions[owner]}:{path}"], cwd=self.root,
        )
        digest = hashlib.sha256(data).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError(f"Frozen input hash mismatch: {owner}:{path}")
        self.hashes[f"{owner}:{path}"] = digest
        return data

    def json(self, owner: str, path: str) -> dict:
        return json.loads(self.read(owner, path))


def table(data: bytes) -> pd.DataFrame:
    frame = pd.read_csv(io.BytesIO(data), sep="\t")
    if "sequence_id" not in frame or frame.sequence_id.isna().any() or not frame.sequence_id.is_unique:
        raise ValueError("Input table requires unique nonmissing sequence IDs")
    return frame


def validate_partial_distance(distance: np.ndarray, ids: list[str]) -> None:
    if len(set(ids)) != len(ids) or not all(isinstance(x, str) and x.strip() for x in ids):
        raise ValueError("Distance IDs must be unique nonempty strings")
    if distance.shape != (len(ids), len(ids)):
        raise ValueError("Distance shape does not match IDs")
    if np.isinf(distance).any():
        raise ValueError("Infinite distance is invalid")
    finite = distance[np.isfinite(distance)]
    if ((finite < 0) | (finite > 1)).any():
        raise ValueError("Distances must be in [0, 1]")
    if not np.allclose(distance, distance.T, equal_nan=True, rtol=0, atol=1e-12):
        raise ValueError("Distances and missing entries must be symmetric")
    if not np.allclose(np.diag(distance), 0, rtol=0, atol=1e-12):
        raise ValueError("Distance diagonal must be zero")


def distance_statistics(values: np.ndarray) -> dict:
    """Observed mean and bounds; missing pairs never receive synthetic distances."""
    values = np.asarray(values, dtype=float)
    observed = values[np.isfinite(values)]
    total, count = len(values), len(observed)
    return {
        "possible_pairs": total,
        "resolved_pairs": count,
        "resolved_fraction": count / total if total else 0.0,
        "observed_mean": float(observed.mean()) if count else None,
        "lower_bound": float(observed.sum() / total) if total else None,
        "upper_bound": float((observed.sum() + total - count) / total) if total else None,
    }


def prepare_scored_pool(frame: pd.DataFrame, distance: np.ndarray, distance_ids: list[str],
                        minimum_resolved_fraction: float) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    if not 0 <= minimum_resolved_fraction <= 1:
        raise ValueError("Minimum resolved fraction must be in [0, 1]")
    mask = eligibility_mask(frame)
    validate_partial_distance(distance, distance_ids)
    pool = frame.loc[mask].copy().reset_index(drop=True)
    lookup = {value: i for i, value in enumerate(distance_ids)}
    if not set(pool.sequence_id) <= set(lookup):
        raise ValueError("Eligible candidates missing from distance matrix")
    positions = [lookup[value] for value in pool.sequence_id]
    pool_distance = distance[np.ix_(positions, positions)]
    statistics = [distance_statistics(np.delete(row, i)) for i, row in enumerate(pool_distance)]
    pool["uniqueness_score"] = [row["observed_mean"] for row in statistics]
    pool["uniqueness_resolved_fraction"] = [row["resolved_fraction"] for row in statistics]
    pool["uniqueness_lower_bound"] = [row["lower_bound"] for row in statistics]
    pool["uniqueness_upper_bound"] = [row["upper_bound"] for row in statistics]
    pool["uniqueness_reference_pool_size"] = len(pool)
    ready = (pool["global_distance_status"].eq("resolved")
             & np.isfinite(pool[METRICS].to_numpy(float)).all(axis=1)
             & pool.uniqueness_resolved_fraction.ge(minimum_resolved_fraction))
    pool["score_ready"] = ready
    pool["score_exclusion_reason"] = np.where(ready, "", "unresolved_metric_or_insufficient_pair_coverage")
    indices = np.flatnonzero(ready.to_numpy())
    return pool.iloc[indices].reset_index(drop=True), pool_distance[np.ix_(indices, indices)], pool


def quality(selected: pd.DataFrame, distance: np.ndarray, position: dict[str, int]) -> dict:
    indices = [position[x] for x in selected.sequence_id]
    sub = distance[np.ix_(indices, indices)]
    stats = distance_statistics(sub[np.triu_indices(len(sub), 1)])
    result = {"selected_count": len(selected), "manual_review_count": int(selected.c_manual_review.sum()),
              **{f"mean_{metric}": float(selected[metric].mean()) for metric in METRICS}}
    result.update({f"subset_diversity_{key}": value for key, value in stats.items()})
    return result


def rank_pool(frame: pd.DataFrame, budget: int) -> dict[str, pd.DataFrame]:
    return {
        "weighted_sum": weighted_ranking(frame, METRICS, WEIGHTS).head(budget),
        "pareto": pareto_ranking(frame, METRICS).head(budget),
        "dimension_round_robin": round_robin_top_k(frame, METRICS, budget),
    }


def run_evaluation(root: Path, config_path: Path, output: Path | None = None) -> dict:
    config = json.loads(config_path.read_text())
    for key in ("weight_samples", "random_samples", "primary_k"):
        if type(config[key]) is not int or config[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if not config["budgets"] or any(type(k) is not int or k < 1 for k in config["budgets"]):
        raise ValueError("budgets must contain positive integers")
    destination = output or root / config["output_dir"]
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("Output directory must be new or empty")
    inputs = FrozenInputs(root, {"b": config["b_revision"], "c": config["c_revision"]})
    cm = inputs.json("c", f"{C_BASE}/qc_similarity_manifest.json")
    bm = inputs.json("b", f"{B_SCORE}/manifest.json")
    fm = inputs.json("c", f"{C_FINAL}/manifest.json")
    if inputs.hashes[f"c:{C_BASE}/qc_similarity_manifest.json"] != fm["source_manifest_sha256"]["base"]:
        raise ValueError("C consolidated/base manifest mismatch")
    for name, folder in [("domains", "results/gv02_03_c_domains_member_a_v1"),
                         ("tm", "results/gv02_03_c_tm_member_a_v1")]:
        inputs.read("c", f"{folder}/manifest.json", fm["source_manifest_sha256"][name])
    # Verify the entire B package used by C, not just filenames or counts.
    b_data = {name: inputs.read("b", spec["path"], spec["sha256"])
              for name, spec in cm["inputs"].items()}
    for name in ("candidates", "reference_fasta", "pool_audit"):
        if hashlib.sha256(b_data[name]).hexdigest() != fm["frozen_input_sha256"][name]:
            raise ValueError(f"C review frozen input mismatch: {name}")
    for path, digest in bm["inputs"].items():
        inputs.read("b", path, digest)
    score_path = f"{B_SCORE}/tables/candidate_scores.tsv"
    scores = table(inputs.read("b", score_path, bm["outputs"][score_path]))
    review = table(inputs.read("c", f"{C_FINAL}/candidate_c_review.tsv",
                               fm["output_sha256"]["candidate_c_review.tsv"]))
    pool_audit = table(b_data["pool_audit"])
    metadata = table(b_data["generation_metadata"])
    if set(review.sequence_id) != set(pool_audit.sequence_id) or set(review.sequence_id) != set(metadata.sequence_id):
        raise ValueError("B/C full candidate ID sets differ")
    b_index = pool_audit.set_index("sequence_id").loc[review.sequence_id]
    for left, right in [("sequence_sha256", "sequence_sha256"), ("pool", "pool"),
                        ("family_status", "family_family_status"), ("qc_pass", "qc_pass")]:
        if not np.array_equal(review[left].to_numpy(), b_index[right].to_numpy()):
            raise ValueError(f"B/C candidate evidence conflict: {left}")
    if set(scores.sequence_id) != set(review.loc[review.pool.eq("main_supported_gvpa"), "sequence_id"]):
        raise ValueError("B score IDs must equal frozen family main pool")
    nearest = table(inputs.read("c", f"{C_BASE}/trusted_similarity.tsv", cm["outputs"]["trusted_similarity.tsv"]))
    if set(nearest.sequence_id) != set(review.sequence_id):
        raise ValueError("Nearest-reference ID set mismatch")
    nearest = nearest.set_index("sequence_id").loc[review.sequence_id]
    if not np.allclose(review.global_distance, nearest.distance, equal_nan=True):
        raise ValueError("C nearest-reference evidence mismatch")
    if not np.array_equal(review.global_distance_status.to_numpy(), nearest.distance_status.to_numpy()):
        raise ValueError("C nearest-reference status mismatch")
    distance = np.load(io.BytesIO(inputs.read("c", f"{C_BASE}/main_distance.npy",
                                             cm["outputs"]["main_distance.npy"])), allow_pickle=False)
    ids = json.loads(inputs.read("c", f"{C_BASE}/main_distance_ids.json",
                                cm["outputs"]["main_distance_ids.json"]))
    if set(ids) != set(scores.sequence_id):
        raise ValueError("Main distance ID set differs from B scored pool")
    scoring_columns = [x for x in scores if x not in review and x not in ("novelty_score", "uniqueness_score")]
    frame = review.merge(scores[["sequence_id", *scoring_columns]], on="sequence_id", how="left", validate="one_to_one")
    if not scores.domain_pass.map(lambda x: isinstance(x, (bool, np.bool_))).all():
        raise ValueError("B domain_pass must be boolean")
    frame["domain_evaluated"] = frame.sequence_id.isin(scores.sequence_id)
    frame["domain_pass"] = frame.domain_pass.eq(True)
    frame["novelty_score"] = frame.global_distance.where(
        frame.qc_pass & frame.family_status.eq("supported_gvpa") & frame.global_distance_status.eq("resolved"))
    frame["uniqueness_score"] = np.nan
    extra = [x for x in metadata if x not in frame]
    frame = frame.merge(metadata[["sequence_id", *extra]], on="sequence_id", validate="one_to_one")
    ranking, ranking_distance, pool = prepare_scored_pool(frame, distance, ids, config["minimum_resolved_fraction"])
    destination.mkdir(parents=True, exist_ok=True)
    input_dir = destination / "inputs"
    input_dir.mkdir(exist_ok=True)
    fasta = input_dir / "frozen_candidates.fasta"
    fasta.write_bytes(b_data["candidates"])
    records_list = list(read_fasta(fasta))
    records = {r.identifier: r for r in records_list}
    if len(records) != len(records_list) or set(records) != set(frame.sequence_id):
        raise ValueError("Candidate FASTA IDs differ")
    for row in frame.itertuples():
        if sequence_sha256(records[row.sequence_id].sequence) != row.sequence_sha256:
            raise ValueError(f"Candidate sequence hash mismatch: {row.sequence_id}")
    def save(name: str, data) -> None:
        value = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
        value.to_csv(destination / f"{name}.tsv", sep="\t", index=False, lineterminator="\n", float_format="%.12g")
    audit = frame.merge(pool[["sequence_id", "score_ready", "score_exclusion_reason"]],
                        on="sequence_id", how="left", validate="one_to_one")
    audit["ranking_eligible"] = audit.sequence_id.isin(ranking.sequence_id)
    audit["exclusion_reason"] = np.select(
        [~audit.qc_pass, ~audit.family_status.eq("supported_gvpa"), ~audit.domain_pass,
         ~audit.ranking_eligible],
        ["qc_failed", "family_not_supported_gvpa", "domain_gate_failed", "insufficient_scoring_evidence"],
        default="",
    )
    save("candidate_audit", audit)
    save("eligible_pool_scores", pool)
    save("ranking_pool_scores", ranking)
    write_fasta((records[x] for x in ranking.sequence_id), destination / "ranking_pool.fasta")
    write_fasta((records[x] for x in frame.loc[frame.pool.eq("ambiguous_exploration"), "sequence_id"]), destination / "ambiguous_exploration.fasta")
    write_fasta((records[x] for x in frame.loc[frame.pool.eq("excluded"), "sequence_id"]), destination / "family_qc_excluded.fasta")
    write_fasta((records[x] for x in audit.loc[~audit.ranking_eligible, "sequence_id"]), destination / "ranking_excluded.fasta")
    np.save(destination / "ranking_distance.npy", ranking_distance, allow_pickle=False)
    write_json(destination / "ranking_distance_ids.json", ranking.sequence_id.tolist())
    position = {x: i for i, x in enumerate(ranking.sequence_id)}
    summaries, overlaps, budgets = [], [], []
    selections = {}
    for k in config["budgets"]:
        budgets.append({"requested_k": k, "available_count": len(ranking), "status": "ok" if k <= len(ranking) else "insufficient"})
        if k > len(ranking):
            continue
        variants = rank_pool(ranking, k)
        for strategy, selected in variants.items():
            selected = selected.copy()
            selected.insert(0, "selection_rank", range(1, len(selected) + 1))
            selections[(strategy, k)] = selected
            save(f"{strategy}_top{k}", selected)
            write_fasta((records[x] for x in selected.sequence_id), destination / f"{strategy}_top{k}.fasta")
            summaries.append({"strategy": strategy, "budget": k, **quality(selected, ranking_distance, position)})
        for left in variants:
            for right in variants:
                overlaps.append({"budget": k, "strategy_a": left, "strategy_b": right,
                                 "jaccard": jaccard(variants[left].sequence_id, variants[right].sequence_id)})
    save("budget_audit", budgets)
    save("strategy_summary", summaries)
    save("strategy_overlap", overlaps)
    k = config["primary_k"]
    primary_ok = k <= len(ranking)
    if primary_ok:
        ablation, ablation_ids = ablation_results(ranking, METRICS, WEIGHTS, k)
        save("ablation_summary", ablation)
        write_json(destination / "ablation_selections.json", ablation_ids)
        weights, frequency = weight_robustness(ranking, METRICS, WEIGHTS, k, config["weight_samples"], 0.2, config["seed"])
        save("weight_robustness", weights)
        save("selection_frequency", [{"sequence_id": x, "frequency": value} for x, value in frequency.items()])
        rng = np.random.default_rng(config["seed"])
        random_rows = [quality(ranking.iloc[rng.choice(len(ranking), k, replace=False)], ranking_distance, position)
                       for _ in range(config["random_samples"])]
        save("random_baseline", random_rows)
        comparisons = []
        for row in summaries:
            if row["budget"] != k:
                continue
            for metric in [*(f"mean_{m}" for m in METRICS), "subset_diversity_observed_mean",
                           "subset_diversity_lower_bound", "subset_diversity_resolved_fraction"]:
                values = np.asarray([r[metric] for r in random_rows], dtype=float)
                values = values[np.isfinite(values)]
                observed = row[metric]
                comparisons.append({"strategy": row["strategy"], "metric": metric, "observed": observed,
                    "valid_random_samples": len(values), "random_mean": float(values.mean()) if len(values) else None,
                    "empirical_p_random_ge_observed": float((np.sum(values >= observed) + 1) / (len(values) + 1))
                    if len(values) and observed is not None else None})
        save("random_comparison", comparisons)
    summary_path = f"{B_SCORE}/summary.json"
    calibration = json.loads(inputs.read("b", summary_path, bm["outputs"][summary_path]))["domain_calibration"]
    baseline = weighted_ranking(ranking, METRICS, WEIGHTS).head(k).sequence_id.tolist() if primary_ok else []
    sensitivity = []
    for score_factor in (0.8, 1.0, 1.2):
        for offset in (-0.1, 0.0, 0.1):
            for fraction in config["coverage_sensitivity"]:
                changed = frame.copy()
                changed["domain_pass"] = (changed.domain_score.ge(calibration["domain_score_cutoff"] * score_factor)
                    & changed.model_coverage.ge(np.clip(calibration["model_coverage_cutoff"] + offset, 0, 1)))
                alternative, _, _ = prepare_scored_pool(changed, distance, ids, fraction)
                enough = len(alternative) >= k
                chosen = weighted_ranking(alternative, METRICS, WEIGHTS).head(k).sequence_id.tolist() if enough else []
                sensitivity.append({"score_factor": score_factor, "coverage_offset": offset,
                    "minimum_resolved_fraction": fraction, "ranking_count": len(alternative),
                    "status": "ok" if enough else "insufficient", "jaccard": jaccard(baseline, chosen) if enough and primary_ok else None})
    save("threshold_coverage_sensitivity", sensitivity)
    old_comparison = []
    for (strategy, budget), selected in selections.items():
        path = f"{B_SCORE}/selections/{strategy}_top{budget}.tsv"
        old = table(inputs.read("b", path))
        old_comparison.append({"strategy": strategy, "budget": budget,
                               "jaccard_with_b_local_alignment": jaccard(selected.sequence_id, old.sequence_id)})
    save("b_baseline_comparison", old_comparison)
    for method in ("pearson", "spearman"):
        ranking[METRICS].corr(method=method).to_csv(destination / f"correlation_{method}.tsv", sep="\t", lineterminator="\n")
    summary = {"status": "completed_on_frozen_provisional_batch", "retrained": False,
        "funnel": eligibility_summary(frame), "ranking_count": len(ranking),
        "score_excluded_count": len(pool) - len(ranking), "primary_budget_status": "ok" if primary_ok else "insufficient",
        "minimum_resolved_fraction": config["minimum_resolved_fraction"], "strategy_summary": summaries,
        "interpretation": "Observed-pair means are conditional on reliable alignment; bounds and coverage must be reported together. Warnings are not automatic exclusions. Reference remains provisional."}
    write_json(destination / "summary.json", summary)
    make_figures(ranking, summaries, destination)
    code_paths = [Path(__file__), root / "src/gv_eval/selection.py", root / "src/gv_eval/analysis.py", root / "src/gv_eval/io.py"]
    manifest = {"source_commits": inputs.revisions, "config": config,
        "input_sha256": inputs.hashes, "config_sha256": sha256_file(config_path),
        "implementation_sha256": {p.relative_to(root).as_posix(): sha256_file(p) for p in code_paths},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
        "output_sha256": {p.relative_to(destination).as_posix(): sha256_file(p)
                          for p in sorted(destination.rglob("*")) if p.is_file()}}
    write_json(destination / "manifest.json", manifest)
    return summary


def make_figures(frame: pd.DataFrame, summaries: list[dict], output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for metric, axis in zip(METRICS, axes.flat):
        axis.hist(frame[metric], bins=20)
        axis.set_title(metric)
    fig.tight_layout()
    fig.savefig(output / "metric_distributions.png", dpi=150)
    plt.close(fig)
    if summaries:
        data = pd.DataFrame(summaries)
        fig, axis = plt.subplots(figsize=(10, 5))
        for strategy, rows in data.groupby("strategy"):
            axis.plot(rows.budget, rows.subset_diversity_observed_mean, marker="o", label=strategy)
        axis.set(xlabel="Top-K", ylabel="Mean distance of resolved pairs (conditional)")
        axis.legend()
        fig.tight_layout()
        fig.savefig(output / "strategy_diversity.png", dpi=150)
        plt.close(fig)


def verify_reproducibility(root: Path, config_path: Path, published: Path) -> dict:
    """Repeat the frozen experiment and compare every artifact, including figures."""
    original = json.loads((published / "manifest.json").read_text())
    for name, digest in original["output_sha256"].items():
        if sha256_file(published / name) != digest:
            raise ValueError(f"Published artifact hash mismatch: {name}")
    with tempfile.TemporaryDirectory(prefix="gv-d-reproduction-") as directory:
        repeated = Path(directory) / "results"
        run_evaluation(root, config_path, repeated)
        manifest = json.loads((repeated / "manifest.json").read_text())
        if manifest != original:
            raise ValueError("Reproduction manifest or artifact hashes differ")
    report = {"passed": True, "compared_artifact_count": len(original["output_sha256"]),
              "source_commits": original["source_commits"],
              "manifest_sha256": sha256_file(published / "manifest.json"),
              "scope": "Same frozen inputs, code, configuration and environment; all output hashes including figures match."}
    write_json(published.parent / f"{published.name}_reproducibility.json", report)
    return report
