import json
from pathlib import Path

import pytest

from gv_eval.release import audit_artifacts, run_release, validate_release_config

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("updates", [
    {"budgets": [10, 10, 20]}, {"budgets": [5, 20]}, {"primary_k": 30},
    {"seed": True}, {"weight_samples": 0}, {"random_samples": -1},
    {"minimum_resolved_fraction": float("nan")}, {"coverage_sensitivity": []},
    {"coverage_sensitivity": [0.25, 0.75]}, {"b_revision": ""},
])
def test_release_rejects_invalid_configuration(updates):
    config = json.loads((ROOT / "configs/d_evaluation.json").read_text())
    config.update(updates)
    with pytest.raises(ValueError):
        validate_release_config(config)


def test_release_rejects_missing_configuration():
    with pytest.raises(ValueError, match="Missing release"):
        validate_release_config({})


def test_integrated_frozen_release(tmp_path):
    config_path = ROOT / "configs/d_evaluation.json"
    destination = tmp_path / "batch"
    result = run_release(ROOT, config_path, output=destination)
    assert result["artifact_audit"]["candidate_count"] == 1000
    assert result["artifact_audit"]["hard_gate_count"] == 136
    assert result["ranking_count"] == 114
    assert result["artifact_audit"]["selection_count"] == 9
    assert result["reproducibility"]["passed"]
    acceptance = json.loads((tmp_path / "batch_acceptance.json").read_text())
    assert acceptance["retrained"] is False
    # Tampering must fail before a repeat is performed.
    (destination / "ranking_pool.fasta").write_text("broken")
    with pytest.raises(ValueError, match="hash mismatch"):
        audit_artifacts(destination, json.loads(config_path.read_text()))
