import json
import csv
import subprocess
import sys
from pathlib import Path

import pytest

from gv_eval.io import FastaRecord
from gv_eval.quality import assess_candidates
from gv_eval.io import sha256_file
from gv_eval.pipeline import run_pipeline
import yaml


def assess(sequence, **options):
    return assess_candidates(
        [FastaRecord("q", "q", sequence)], minimum_length=1, maximum_length=1000,
        minimum_unique_residues=1, maximum_single_residue_fraction=1.0,
        **options,
    )["q"]


def test_hydrophobic_runs_have_one_based_inclusive_coordinates_and_warn_only():
    row = assess("deavilfwkavilm", hydrophobic_run_minimum=5)
    assert json.loads(row["hydrophobic_segments"]) == [
        {"start": 3, "end": 8, "length": 6},
        {"start": 10, "end": 14, "length": 5},
    ]
    assert row["hydrophobic_checked"] is True
    assert row["hydrophobic_warning"] is True
    assert row["qc_warnings"] == "hydrophobic_run_warning"
    assert row["qc_reasons"] == row["qc_warnings"]
    assert row["qc_pass"] is True


def test_hydrophobic_threshold_is_inclusive_and_interruptions_split_runs():
    assert assess("AVIL", hydrophobic_run_minimum=5)["hydrophobic_warning"] is False
    assert assess("AVILM", hydrophobic_run_minimum=5)["hydrophobic_warning"] is True
    assert assess("AVDKILM", hydrophobic_run_minimum=5)["hydrophobic_warning"] is False


def test_homopolymers_are_maximal_and_include_terminal_runs():
    row = assess("DAAAAAAGEEEEEE", homopolymer_minimum=6)
    assert json.loads(row["homopolymer_segments"]) == [
        {"start": 2, "end": 7, "length": 6, "residue": "A"},
        {"start": 9, "end": 14, "length": 6, "residue": "E"},
    ]
    assert row["homopolymer_checked"] is True
    assert row["homopolymer_warning"] is True
    assert row["qc_warnings"] == "homopolymer_warning"
    assert row["qc_pass"] is True
    assert not assess("AAAAA", homopolymer_minimum=6)["homopolymer_warning"]


def test_tandem_repeats_report_primitive_motif_and_complete_copies():
    row = assess("KACDACDACDQ", tandem_repeat_minimum_length=9,
                 tandem_repeat_minimum_copies=3, tandem_repeat_maximum_motif_length=6)
    assert json.loads(row["tandem_repeat_segments"]) == [
        {"start": 2, "end": 10, "length": 9, "motif": "ACD", "copies": 3},
    ]
    assert row["tandem_repeat_checked"] is True
    assert row["tandem_repeat_warning"] is True
    assert row["qc_pass"] is True


def test_tandem_repeats_deduplicate_shifted_and_multiple_period_hits():
    row = assess("AC" * 7 + "A", tandem_repeat_minimum_length=12)
    assert json.loads(row["tandem_repeat_segments"]) == [
        {"start": 1, "end": 14, "length": 14, "motif": "AC", "copies": 7},
    ]


@pytest.mark.parametrize("sequence", ["A" * 20, "ACDACD", "ACDACDACD", "ACDACDEACD"])
def test_tandem_repeat_span_and_copy_gates_and_no_homopolymer_duplicates(sequence):
    assert not assess(sequence, tandem_repeat_minimum_length=12)["tandem_repeat_warning"]


def test_tandem_motif_length_limit_is_respected():
    sequence = "ACDE" * 3
    assert not assess(sequence, tandem_repeat_minimum_length=12,
                      tandem_repeat_maximum_motif_length=3)["tandem_repeat_warning"]
    assert assess(sequence, tandem_repeat_minimum_length=12,
                  tandem_repeat_maximum_motif_length=4)["tandem_repeat_warning"]


@pytest.mark.parametrize("sequence", ["", "AXAAAA", "ACDE"])
def test_disabled_checks_are_not_reported_as_passed(sequence):
    row = assess(sequence)
    for prefix in ("hydrophobic", "homopolymer", "tandem_repeat"):
        assert row[f"{prefix}_checked"] is False
        assert row[f"{prefix}_warning"] is False
        assert json.loads(row[f"{prefix}_segments"]) == []


@pytest.mark.parametrize("sequence", ["", "AVILMX", "AAAAAAX"])
def test_invalid_sequences_are_not_pattern_checked(sequence):
    row = assess(sequence, hydrophobic_run_minimum=3, homopolymer_minimum=3,
                 tandem_repeat_minimum_length=6)
    assert row["qc_pass"] is False
    for prefix in ("hydrophobic", "homopolymer", "tandem_repeat"):
        assert row[f"{prefix}_checked"] is False


def test_warning_field_excludes_failures_and_keeps_old_reason_order():
    row = assess("ACDEX", generation_length_cap=5)
    assert row["qc_reasons"] == "nonstandard_amino_acids;generation_cap_warning"
    assert row["qc_warnings"] == "generation_cap_warning"
    assert row["qc_pass"] is False


def test_exact_and_composition_warnings_are_in_separate_field():
    reference = [FastaRecord("a", "a", "AAAA"), FastaRecord("c", "c", "CCCC")]
    row = assess("DEFG", training_records=[FastaRecord("t", "t", "DEFG")],
                 reference_records=reference, composition_outlier_quantile=0.9)
    assert row["qc_warnings"] == "exact_training_match_warning;composition_outlier_warning"
    assert row["qc_pass"] is True
    excluded = assess("DEFG", training_records=[FastaRecord("t", "t", "DEFG")],
                      exact_training_match_action="exclude")
    assert excluded["qc_warnings"] == ""
    assert excluded["qc_reasons"] == "exact_training_match"


@pytest.mark.parametrize("option,value", [
    ("hydrophobic_run_minimum", 0), ("hydrophobic_run_minimum", True),
    ("homopolymer_minimum", 1.5), ("homopolymer_minimum", -1),
    ("tandem_repeat_minimum_length", 0), ("tandem_repeat_minimum_length", "12"),
    ("tandem_repeat_minimum_copies", 1), ("tandem_repeat_minimum_copies", False),
    ("tandem_repeat_maximum_motif_length", 1),
    ("tandem_repeat_maximum_motif_length", float("nan")),
])
def test_invalid_pattern_options_fail_even_with_no_records(option, value):
    with pytest.raises(ValueError, match=option):
        assess_candidates([], minimum_length=1, maximum_length=1000,
                          minimum_unique_residues=1, maximum_single_residue_fraction=1.0,
                          **{option: value})


def test_pipeline_serializes_pattern_evidence_without_filtering_warned_candidate(tmp_path):
    (tmp_path / "input.fa").write_text(">q\nDEAVILM\n", encoding="utf-8")
    config = {
        "seed": 42, "inputs": {"candidates": "input.fa"},
        "outputs": {"processed_dir": "processed", "result_dir": "results"},
        "quality": {"minimum_length": 1, "maximum_length": 100,
                    "minimum_unique_residues": 1, "maximum_single_residue_fraction": 1,
                    "hydrophobic_run_minimum": 5},
        "selection": {"budgets": [2], "primary_k": 2},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ValueError, match="Insufficient QC-passing"):
        run_pipeline(tmp_path, path)
    with (tmp_path / "processed/candidate_qc.tsv").open(encoding="utf-8", newline="") as stream:
        row, = csv.DictReader(stream, delimiter="\t")
    assert row["qc_pass"] == "True"
    assert row["qc_warnings"] == "hydrophobic_run_warning"
    assert json.loads(row["hydrophobic_segments"]) == [{"start": 3, "end": 7, "length": 5}]
    assert "DEAVILM" in (tmp_path / "processed/eligible_candidates.fasta").read_text()


def test_qc_only_command_is_reproducible_and_refuses_overwrite(tmp_path):
    root = Path(__file__).resolve().parents[1]
    (tmp_path / "input.fa").write_text(">q\nDEAVILM\n", encoding="utf-8")
    config = {"inputs": {"candidates": str(tmp_path / "input.fa")},
              "quality": {"minimum_length": 1, "maximum_length": 100,
                          "minimum_unique_residues": 1, "maximum_single_residue_fraction": 1,
                          "hydrophobic_run_minimum": 5}}
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    outputs = [tmp_path / "first", tmp_path / "second"]
    for output in outputs:
        run = subprocess.run([sys.executable, str(root / "experiments/run_quality_audit.py"),
                              "--config", str(path), "--output-dir", str(output)],
                             capture_output=True, text=True)
        assert run.returncode == 0, run.stderr
        summary = json.loads((output / "summary.json").read_text())
        assert summary["input_count"] == summary["pass_count"] == 1
        assert summary["warning_counts"] == {"hydrophobic_run_warning": 1}
        manifest = json.loads((output / "manifest.json").read_text())
        for name, expected in manifest["outputs"].items():
            assert sha256_file(output / name) == expected
    for name in ("candidate_qc.tsv", "eligible_candidates.fasta", "summary.json", "manifest.json"):
        assert (outputs[0] / name).read_bytes() == (outputs[1] / name).read_bytes()
    before = (outputs[1] / "candidate_qc.tsv").read_bytes()
    repeated = subprocess.run([sys.executable, str(root / "experiments/run_quality_audit.py"),
                               "--config", str(path), "--output-dir", str(outputs[1])],
                              capture_output=True, text=True)
    assert repeated.returncode != 0
    assert "not empty" in repeated.stderr
    assert (outputs[1] / "candidate_qc.tsv").read_bytes() == before


def test_committed_pattern_audit_is_consistent_and_does_not_change_eligibility():
    root = Path(__file__).resolve().parents[1]
    output = root / "results/gv02_03_pattern_audit"
    manifest = json.loads((output / "manifest.json").read_text())
    for name, expected in manifest["outputs"].items():
        assert sha256_file(output / name) == expected
    with (output / "candidate_qc.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with (root / "data/processed/gv02_03_qc/candidate_qc.tsv").open(encoding="utf-8", newline="") as stream:
        old_rows = list(csv.DictReader(stream, delimiter="\t"))
    assert {r["sequence_id"]: r["qc_pass"] for r in rows} == {
        r["sequence_id"]: r["qc_pass"] for r in old_rows
    }
    summary = json.loads((output / "summary.json").read_text())
    assert len(rows) == summary["input_count"] == 200
    assert sum(r["qc_pass"] == "True" for r in rows) == summary["pass_count"] == 170
    assert sum(r["homopolymer_warning"] == "True" for r in rows) == 10
    assert summary["warning_counts"]["homopolymer_warning"] == 10
