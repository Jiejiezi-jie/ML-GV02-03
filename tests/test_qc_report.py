from __future__ import annotations

import csv
import importlib
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from gv_eval.io import sha256_file, write_json, write_tsv
from gv_eval.nearest_reference import METRIC_FIELDS


def report_module():
    return importlib.import_module("gv_eval.qc_report")


def qc(identifier, passed=True, reasons="", warnings="", **changes):
    return {"sequence_id": identifier, "length": "20", "qc_pass": str(passed),
            "qc_reasons": reasons, "qc_warnings": warnings, "hydrophobic_checked": "True",
            "custom_note": "原始证据", **changes}


def nearest(identifier, resolved=True, ties=1, **changes):
    row = {"sequence_id": identifier, "query_length": "20",
           "closest_reference_id": "r1" if resolved else "",
           "closest_reference_ids": json.dumps([f"r{i+1}" for i in range(ties)]) if resolved else "[]",
           "closest_reference_tie_count": str(ties if resolved else 0),
           "distance_status": "resolved" if resolved else "unresolved",
           "distance_reason": "" if resolved else "no_reliable_reference",
           "reference_count": "3", "resolved_reference_count": str(ties if resolved else 0),
           **{name: "" for name in METRIC_FIELDS}}
    if resolved:
        row.update(target_length="20", alignment_score="10.0", alignment_length="20",
                   aligned_pairs="20", identical_residues="10", identity="0.5", query_coverage="1.0",
                   target_coverage="1.0", effective_identity="0.5", distance="0.5")
    return row | changes


def example_rows():
    return (
        [qc("clean"), qc("warn", reasons="composition_outlier_warning", warnings="composition_outlier_warning"),
         qc("fail", False, "length_above_maximum;generation_cap_warning", "generation_cap_warning"),
         qc("unknown")],
        [nearest("unknown", False), nearest("fail", False), nearest("warn", ties=2), nearest("clean")],
    )


def test_keyed_join_preserves_all_evidence_and_original_qc_order():
    qc_rows, near_rows = example_rows()
    rows, summary = report_module().merge_qc_rows(qc_rows, near_rows)
    assert [r["sequence_id"] for r in rows] == ["clean", "warn", "fail", "unknown"]
    assert rows[0]["custom_note"] == "原始证据"
    assert rows[2]["qc_reasons"] == "length_above_maximum;generation_cap_warning"
    assert rows[2]["qc_failure_reasons"] == "length_above_maximum"
    assert rows[2]["qc_pass"] == "False"
    assert rows[3]["distance"] == ""
    assert summary["candidate_count"] == 4


def test_review_is_union_not_sum_and_cross_counts_are_explicit():
    rows, summary = report_module().merge_qc_rows(*example_rows())
    assert [r["review_required"] for r in rows] == [False, True, True, True]
    assert rows[1]["review_reasons"] == "qc_warning;reference_tie"
    assert rows[2]["review_reasons"] == "qc_failure;qc_warning;reference_unresolved"
    assert summary["review_count"] == 3
    assert summary["no_review_trigger_count"] == 1
    assert summary["cross_counts"] == {"pass_resolved": 2, "pass_unresolved": 1,
                                       "fail_resolved": 0, "fail_unresolved": 1}
    assert summary["failure_reason_counts"] == {"length_above_maximum": 1}
    assert summary["warning_counts"] == {"composition_outlier_warning": 1, "generation_cap_warning": 1}


def test_unchecked_is_visible_not_claimed_as_passed():
    rows, summary = report_module().merge_qc_rows(
        [qc("q", hydrophobic_checked="False")], [nearest("q")])
    assert rows[0]["qc_checks_not_run"] == "hydrophobic_checked"
    assert summary["not_checked_counts"] == {"hydrophobic_checked": 1}
    assert summary["not_assessed_by_report"] == ["eos_metadata", "transmembrane_topology", "domain_architecture", "family_assignment"]


@pytest.mark.parametrize("which,mutation", [
    ("qc", "duplicate"), ("near", "duplicate"), ("qc", "empty"), ("near", "empty"),
    ("qc", "missing"), ("near", "extra"),
])
def test_bad_id_sets_fail_without_dropping_rows(which, mutation):
    left, right = example_rows()
    target = left if which == "qc" else right
    if mutation == "duplicate":
        target.append(dict(target[0]))
    elif mutation == "empty":
        target[0]["sequence_id"] = ""
    elif mutation == "missing":
        target.pop()
    else:
        target.append(nearest("extra"))
    with pytest.raises(ValueError):
        report_module().merge_qc_rows(left, right)


@pytest.mark.parametrize("qc_changes,near_changes", [
    ({"length": "21"}, {}), ({"qc_pass": "false"}, {}),
    ({"qc_pass": "False"}, {}), ({"qc_reasons": "length_above_maximum"}, {}),
    ({"qc_warnings": "orphan_warning"}, {}), ({"hydrophobic_checked": "maybe"}, {}),
    ({}, {"length": "20"}), ({}, {"distance_status": "maybe"}),
    ({}, {"closest_reference_ids": "not json"}), ({}, {"closest_reference_tie_count": "2"}),
    ({}, {"closest_reference_id": "wrong"}), ({}, {"distance": "nan"}),
    ({}, {"distance": "1.5"}), ({}, {"reference_count": "0"}),
])
def test_invalid_fields_and_conflicting_columns_are_rejected(qc_changes, near_changes):
    with pytest.raises(ValueError):
        report_module().merge_qc_rows([qc("q", **qc_changes)], [nearest("q", **near_changes)])


def test_unresolved_must_not_carry_representative_metrics():
    with pytest.raises(ValueError):
        report_module().merge_qc_rows([qc("q")], [nearest("q", False, distance="1.0")])


@pytest.mark.parametrize("changes", [
    {"identical_residues": "999"}, {"aligned_pairs": "999"},
    {"alignment_length": "21"}, {"distance": "0.0"}, {"identity": "0.7"},
    {"query_coverage": "0.8"}, {"target_coverage": "0.8"}, {"effective_identity": "0.4"},
])
def test_resolved_metric_relationships_are_validated(changes):
    with pytest.raises(ValueError, match="Inconsistent"):
        report_module().merge_qc_rows([qc("q")], [nearest("q", **changes)])


def test_metric_validation_allows_normal_float_rounding_without_rewriting():
    rows, _ = report_module().merge_qc_rows([qc("q")], [nearest("q", identity="0.5000000000000001")])
    assert rows[0]["identity"] == "0.5000000000000001"


def make_runs(tmp_path):
    qc_dir, near_dir = tmp_path / "qc", tmp_path / "near"
    qc_dir.mkdir()
    near_dir.mkdir()
    left, right = example_rows()
    write_tsv(qc_dir / "candidate_qc.tsv", left, list(left[0]))
    write_tsv(near_dir / "nearest_reference.tsv", right, list(right[0]))
    write_json(qc_dir / "manifest.json", {
        "inputs": {"candidates": {"sha256": "a" * 64}},
        "outputs": {"candidate_qc.tsv": sha256_file(qc_dir / "candidate_qc.tsv")},
    })
    write_json(near_dir / "manifest.json", {
        "inputs": {"candidate_fasta": {"sha256": "a" * 64}},
        "outputs": {"nearest_reference.tsv": sha256_file(near_dir / "nearest_reference.tsv")},
    })
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"qc_run_dir": "qc", "nearest_run_dir": "near"}), encoding="utf-8")
    return config, qc_dir, near_dir


def test_run_writes_reproducible_report_and_review_subset_without_overwriting(tmp_path):
    config, _, _ = make_runs(tmp_path)
    module = report_module()
    for output in ("first", "second"):
        summary = module.run_qc_report(tmp_path, config, output)
        assert summary["review_count"] == 3
        with (tmp_path / output / "manual_review.tsv").open(encoding="utf-8", newline="") as stream:
            assert [r["sequence_id"] for r in csv.DictReader(stream, delimiter="\t")] == ["warn", "fail", "unknown"]
        content = (tmp_path / output / "report.md").read_text(encoding="utf-8")
        assert "EOS" in content and "不等于完整 QC 通过" in content
        manifest = json.loads((tmp_path / output / "manifest.json").read_text(encoding="utf-8"))
        for name, expected in manifest["outputs"].items():
            assert sha256_file(tmp_path / output / name) == expected
    for name in ("candidate_qc_summary.tsv", "manual_review.tsv", "summary.json", "report.md", "manifest.json"):
        assert (tmp_path / "first" / name).read_bytes() == (tmp_path / "second" / name).read_bytes()
    before = (tmp_path / "first/report.md").read_bytes()
    with pytest.raises(ValueError, match="not empty"):
        module.run_qc_report(tmp_path, config, "first")
    assert (tmp_path / "first/report.md").read_bytes() == before


@pytest.mark.parametrize("problem", ["table_tamper", "candidate_hash_mismatch", "malformed_tsv", "missing_column", "duplicate_column"])
def test_manifest_or_table_problem_blocks_publication(tmp_path, problem):
    config, qc_dir, near_dir = make_runs(tmp_path)
    if problem == "candidate_hash_mismatch":
        path = near_dir / "manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["inputs"]["candidate_fasta"]["sha256"] = "b" * 64
        write_json(path, manifest)
    else:
        table = qc_dir / "candidate_qc.tsv"
        content = table.read_text(encoding="utf-8")
        if problem == "table_tamper":
            table.write_text(content + "unexpected", encoding="utf-8")
        else:
            if problem == "malformed_tsv":
                content += "incomplete\n"
            elif problem == "missing_column":
                content = content.replace("qc_warnings", "renamed_warnings", 1)
            else:
                content = content.replace("custom_note", "qc_pass", 1)
            table.write_text(content, encoding="utf-8")
            manifest = json.loads((qc_dir / "manifest.json").read_text(encoding="utf-8"))
            manifest["outputs"]["candidate_qc.tsv"] = sha256_file(table)
            write_json(qc_dir / "manifest.json", manifest)
    with pytest.raises(ValueError):
        report_module().run_qc_report(tmp_path, config, "output")
    assert not (tmp_path / "output").exists()


def test_command_runs_and_empty_review_retains_headers(tmp_path):
    config, qc_dir, near_dir = make_runs(tmp_path)
    left, right = [qc("clean")], [nearest("clean")]
    for directory, filename, rows in ((qc_dir, "candidate_qc.tsv", left),
                                       (near_dir, "nearest_reference.tsv", right)):
        write_tsv(directory / filename, rows, list(rows[0]))
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        manifest["outputs"][filename] = sha256_file(directory / filename)
        write_json(directory / "manifest.json", manifest)
    # Absolute input directories permit invoking the repository entry point.
    config.write_text(yaml.safe_dump({"qc_run_dir": str(qc_dir), "nearest_run_dir": str(near_dir)}), encoding="utf-8")
    root = Path(__file__).resolve().parents[1]
    run = subprocess.run([sys.executable, str(root / "experiments/run_qc_report.py"),
                          "--config", str(config), "--output-dir", str(tmp_path / "output")],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    with (tmp_path / "output/manual_review.tsv").open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert "review_reasons" in reader.fieldnames
        assert list(reader) == []


def test_input_changed_during_merge_prevents_output(tmp_path, monkeypatch):
    config, qc_dir, _ = make_runs(tmp_path)
    module = report_module()
    real_merge = module.merge_qc_rows

    def merge_then_change(*args):
        result = real_merge(*args)
        (qc_dir / "manifest.json").write_text("{}", encoding="utf-8")
        return result

    monkeypatch.setattr(module, "merge_qc_rows", merge_then_change)
    with pytest.raises(RuntimeError, match="Input changed"):
        module.run_qc_report(tmp_path, config, "output")
    assert not (tmp_path / "output/candidate_qc_summary.tsv").exists()


def test_output_populated_during_merge_is_preserved(tmp_path, monkeypatch):
    config, _, _ = make_runs(tmp_path)
    module = report_module()
    real_merge = module.merge_qc_rows

    def competing_writer(*args):
        result = real_merge(*args)
        (tmp_path / "output").mkdir()
        (tmp_path / "output/report.md").write_text("other report", encoding="utf-8")
        return result

    monkeypatch.setattr(module, "merge_qc_rows", competing_writer)
    with pytest.raises(ValueError, match="not empty"):
        module.run_qc_report(tmp_path, config, "output")
    assert (tmp_path / "output/report.md").read_text() == "other report"


def test_development_report_preserves_sources_and_audits_all_outputs():
    root = Path(__file__).resolve().parents[1]
    output = root / "results/gv02_03_qc_report"

    def table(path):
        with path.open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    rows = table(output / "candidate_qc_summary.tsv")
    review = table(output / "manual_review.tsv")
    qc_rows = table(root / "results/gv02_03_pattern_audit/candidate_qc.tsv")
    near_rows = table(root / "results/gv02_03_nearest_reference_dev/nearest_reference.tsv")
    by_id = {row["sequence_id"]: row for row in rows}
    assert len(rows) == len(by_id) == 200
    assert [r["sequence_id"] for r in rows] == [r["sequence_id"] for r in qc_rows]
    for original in qc_rows + near_rows:
        assert all(by_id[original["sequence_id"]][key] == value for key, value in original.items())
    assert review == [row for row in rows if row["review_required"] == "True"]
    assert len(review) == 180
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["cross_counts"] == {"pass_resolved": 165, "pass_unresolved": 5,
                                       "fail_resolved": 0, "fail_unresolved": 30}
    assert summary["no_review_trigger_count"] == 20
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["outputs"].items():
        assert sha256_file(output / name) == expected
