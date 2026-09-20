from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from gv_eval.io import sha256_file


ROOT = Path(__file__).resolve().parents[1]


def test_qc_run_counts_and_filtering_are_auditable():
    processed = ROOT / "data/processed/gv02_03_qc"
    result = ROOT / "results/gv02_03_qc"
    qc = pd.read_csv(processed / "candidate_qc.tsv", sep="\t")
    scores = pd.read_csv(result / "tables/candidate_scores.tsv", sep="\t")
    summary = json.loads((result / "summary.json").read_text(encoding="utf-8"))

    assert len(qc) == 200
    assert qc["sequence_id"].is_unique
    assert int(qc["qc_pass"].sum()) == 170
    assert int((~qc["qc_pass"]).sum()) == 30
    assert int(qc["composition_outlier_warning"].sum()) == 145
    assert int(qc["exact_training_match"].sum()) == 0
    assert int(qc["exact_reference_match"].sum()) == 0
    assert len(scores) == 170
    assert set(scores["sequence_id"]) == set(qc.loc[qc["qc_pass"], "sequence_id"])
    assert summary["quality_control"]["pass_count"] == 170
    assert summary["quality_control"]["failure_count"] == 30


def test_qc_manifest_hashes_match_current_artifacts():
    manifest = json.loads(
        (ROOT / "results/gv02_03_qc/manifest.json").read_text(encoding="utf-8")
    )
    for section in ("inputs", "outputs"):
        for relative, expected in manifest[section].items():
            assert sha256_file(ROOT / relative) == expected


def test_qc_reproducibility_record_is_explicit_about_hmmer_metadata():
    record = json.loads(
        (ROOT / "results/gv02_03_qc/reproducibility_check.json").read_text(
            encoding="utf-8"
        )
    )
    assert record["runs"] == 2
    assert record["derived_outputs_identical"] is True
    assert record["byte_identical_artifacts"] == 55
    assert len(record["non_identical_raw_artifacts"]) == 4
    assert record["manifest_sha256_run_1"] == record["manifest_sha256_run_2"]
    assert record["manifest_sha256_run_2"] == sha256_file(
        ROOT / "results/gv02_03_qc/manifest.json"
    )
