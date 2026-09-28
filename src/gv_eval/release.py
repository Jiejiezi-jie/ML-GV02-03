"""Integrated acceptance of the frozen A/B/C/D development batch."""
from __future__ import annotations

import json
import math
import platform
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .d_evaluation import METRICS, run_evaluation, verify_reproducibility
from .io import read_fasta, sha256_file, write_json


def validate_release_config(config: dict) -> None:
    required = {"b_revision", "c_revision", "seed", "budgets", "primary_k", "minimum_resolved_fraction",
                "coverage_sensitivity", "weight_samples", "random_samples", "output_dir"}
    missing = required - config.keys()
    if missing:
        raise ValueError(f"Missing release configuration fields: {sorted(missing)}")
    for name in ("seed", "primary_k", "weight_samples", "random_samples"):
        value = config[name]
        if type(value) is not int or value < (0 if name == "seed" else 1):
            raise ValueError(f"Invalid {name}: expected integer")
    budgets = config["budgets"]
    if (not isinstance(budgets, list) or not budgets or
        any(type(k) is not int or k < 1 for k in budgets) or len(set(budgets)) != len(budgets)):
        raise ValueError("budgets must be unique positive integers")
    if config["primary_k"] not in budgets:
        raise ValueError("primary_k must be present in budgets")
    # B published comparisons exist for these exact budgets; do not silently use other lists.
    if not set(budgets) <= {10, 20, 50}:
        raise ValueError("Frozen B baseline supports only budgets 10, 20 and 50")
    fractions = config["coverage_sensitivity"]
    if not isinstance(fractions, list) or not fractions:
        raise ValueError("coverage_sensitivity must be a nonempty list")
    for value in [config["minimum_resolved_fraction"], *fractions]:
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Resolved fractions must be finite numbers in [0, 1]")
    if config["minimum_resolved_fraction"] not in fractions:
        raise ValueError("coverage_sensitivity must include the baseline fraction")
    for name in ("b_revision", "c_revision", "output_dir"):
        if not isinstance(config[name], str) or not config[name].strip():
            raise ValueError(f"{name} must be a nonempty string")


def audit_artifacts(output: Path, config: dict) -> dict:
    """Independently verify published scores, partitions, matrices and experiments."""
    manifest = json.loads((output / "manifest.json").read_text())
    if manifest["config"] != config:
        raise ValueError("Manifest configuration differs from requested configuration")
    for name, digest in manifest["output_sha256"].items():
        if sha256_file(output / name) != digest:
            raise ValueError(f"Output hash mismatch: {name}")
    def read(name: str) -> pd.DataFrame:
        try:
            return pd.read_csv(output / f"{name}.tsv", sep="\t")
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    audit, pool, rank = (read(name) for name in ("candidate_audit", "eligible_pool_scores", "ranking_pool_scores"))
    for frame in (audit, pool, rank):
        if frame.sequence_id.isna().any() or not frame.sequence_id.is_unique:
            raise ValueError("Artifact IDs must be unique and nonmissing")
    if set(audit.loc[audit.ranking_eligible, "sequence_id"]) != set(rank.sequence_id):
        raise ValueError("Ranking/audit IDs differ")
    if not set(rank.sequence_id) <= set(pool.sequence_id) <= set(audit.sequence_id):
        raise ValueError("Candidate pool containment violated")
    if not (rank.qc_pass & rank.domain_pass & rank.family_status.eq("supported_gvpa")).all():
        raise ValueError("Ranking includes a hard-gate failure")
    values = rank[METRICS].to_numpy(float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("Ranking metrics must be finite and within [0, 1]")
    if not rank.uniqueness_resolved_fraction.ge(config["minimum_resolved_fraction"]).all():
        raise ValueError("Ranking violates resolved-pair requirement")
    ids = json.loads((output / "ranking_distance_ids.json").read_text())
    matrix = np.load(output / "ranking_distance.npy", allow_pickle=False)
    if ids != rank.sequence_id.tolist() or matrix.shape != (len(ids), len(ids)):
        raise ValueError("Ranking matrix order/shape mismatch")
    from .d_evaluation import validate_partial_distance
    validate_partial_distance(matrix, ids)
    selected_count = 0
    for row in read("strategy_summary").itertuples():
        stem = f"{row.strategy}_top{row.budget}"
        selected = read(stem)
        fasta_ids = [record.identifier for record in read_fasta(output / f"{stem}.fasta")]
        if (len(selected) != row.budget or not selected.sequence_id.is_unique or
            not set(selected.sequence_id) <= set(ids) or fasta_ids != selected.sequence_id.tolist()):
            raise ValueError(f"Invalid selection: {stem}")
        selected_count += 1
    expected_selections = 3 * sum(k <= len(rank) for k in config["budgets"])
    if selected_count != expected_selections:
        raise ValueError("Missing selection outputs")
    if config["primary_k"] <= len(rank):
        if len(read("weight_robustness")) != config["weight_samples"]:
            raise ValueError("Weight experiment count mismatch")
        if len(read("random_baseline")) != config["random_samples"]:
            raise ValueError("Random experiment count mismatch")
    if len(read("threshold_coverage_sensitivity")) != 9 * len(config["coverage_sensitivity"]):
        raise ValueError("Threshold/coverage experiment count mismatch")
    return {"passed": True, "hashed_artifacts": len(manifest["output_sha256"]),
            "candidate_count": len(audit), "hard_gate_count": len(pool), "ranking_count": len(rank),
            "selection_count": selected_count}


def run_release(root: Path, config_path: Path, *, output: Path | None = None,
                verify_existing: bool = False) -> dict:
    config = json.loads(config_path.read_text())
    validate_release_config(config)
    destination = output or root / config["output_dir"]
    # Preflight dependency versions and frozen Git objects before creating outputs.
    for name in ("b_revision", "c_revision"):
        subprocess.run(["git", "cat-file", "-e", f"{config[name]}^{{commit}}"], cwd=root, check=True)
    if verify_existing:
        audit = audit_artifacts(destination, config)
        reproduction = verify_reproducibility(root, config_path, destination)
        summary = json.loads((destination / "summary.json").read_text())
    else:
        summary = run_evaluation(root, config_path, destination)
        audit = audit_artifacts(destination, config)
        reproduction = verify_reproducibility(root, config_path, destination)
    paths = [root / "experiments/run_full_experiment.py", Path(__file__)]
    report = {"status": "passed_frozen_batch_engineering_acceptance", "retrained": False,
              "artifact_audit": audit, "reproducibility": reproduction,
              "integration_sha256": {p.relative_to(root).as_posix(): sha256_file(p) for p in paths},
              "python": platform.python_version(),
              "scope": "Frozen upstream artifacts through scoring and selection. Does not replay VAE training or generation, or validate biological function."}
    write_json(destination.parent / f"{destination.name}_acceptance.json", report)
    return {"status": report["status"], "ranking_count": summary["ranking_count"],
            "artifact_audit": audit, "reproducibility": reproduction}
