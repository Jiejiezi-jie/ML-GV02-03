"""QC at the real pipeline entry, before external bioinformatics tools run."""
import csv
import json

import pytest
import yaml

from gv_eval.io import read_fasta, sha256_file
from gv_eval.pipeline import run_pipeline


@pytest.mark.parametrize("action", ["warn", "exclude"])
def test_pipeline_audits_matches_and_filters_before_tools(tmp_path, action):
    candidates = tmp_path / "candidates.fasta"
    candidates.write_text(">q\nACDE\n>near\nACDF\n", encoding="utf-8")
    (tmp_path / "train.fasta").write_text(">training_id\nACDE\n", encoding="utf-8")
    (tmp_path / "reference.fasta").write_text(">reference_id\nACDE\n", encoding="utf-8")
    before = sha256_file(candidates)
    config = {
        "seed": 42,
        "inputs": {"candidates": "candidates.fasta", "qc_training_fasta": "train.fasta",
                   "qc_reference_fasta": "reference.fasta"},
        "outputs": {"processed_dir": "processed", "result_dir": "results"},
        "quality": {"minimum_length": 4, "maximum_length": 10,
                    "minimum_unique_residues": 3, "maximum_single_residue_fraction": 0.7,
                    "exact_training_match_action": action},
        # Intentionally above the input size: stop before external tool invocation.
        "selection": {"budgets": [3], "primary_k": 3},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ValueError, match="Insufficient QC-passing"):
        run_pipeline(tmp_path, path)
    with (tmp_path / "processed/candidate_qc.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["exact_training_match"] == "True"
    assert rows[0]["exact_reference_match"] == "True"
    assert json.loads(rows[0]["training_match_ids"]) == ["training_id"]
    assert rows[0]["qc_pass"] == str(action == "warn")
    assert [r.identifier for r in read_fasta(tmp_path / "processed/eligible_candidates.fasta")] == (
        ["q", "near"] if action == "warn" else ["near"])
    assert sha256_file(candidates) == before
