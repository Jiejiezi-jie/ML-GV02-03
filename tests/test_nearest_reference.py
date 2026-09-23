from __future__ import annotations

import csv
import importlib
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from gv_eval.io import FastaRecord, read_fasta, sha256_file
from gv_eval.similarity import SimilarityConfig, align_pair


def rec(identifier, sequence):
    return FastaRecord(identifier, identifier, sequence)


def match(queries, references, **settings):
    module = importlib.import_module("gv_eval.nearest_reference")
    return module.nearest_reference_matches(
        queries, references, SimilarityConfig(minimum_aligned_residues=1, **settings)
    )


def test_selects_smallest_reliable_distance_and_preserves_candidate_order():
    rows = match([rec("q2", "ACDEFGHIK"), rec("q1", "ACDEYGHIK")],
                 [rec("far", "ACDEWGHIK"), rec("exact", "acdefghik")])
    assert [row["sequence_id"] for row in rows] == ["q2", "q1"]
    row = rows[0]
    assert row["closest_reference_id"] == "exact"
    assert json.loads(row["closest_reference_ids"]) == ["exact"]
    assert row["distance"] == 0
    assert row["identity"] == row["query_coverage"] == row["target_coverage"] == 1
    assert row["distance_status"] == "resolved"
    assert row["resolved_reference_count"] == row["reference_count"] == 2


def test_no_reliable_match_leaves_id_distance_and_metrics_blank():
    row, = match([rec("q", "AAAAAAAA")], [rec("w", "WWWWWWWW")])
    assert row["closest_reference_id"] is None
    assert json.loads(row["closest_reference_ids"]) == []
    assert row["closest_reference_tie_count"] == 0
    for field in ("distance", "identity", "query_coverage", "target_coverage", "alignment_score"):
        assert row[field] is None
    assert row["distance_status"] == "unresolved"
    assert row["distance_reason"] == "no_reliable_reference"
    assert row["resolved_reference_count"] == 0


def test_unreliable_higher_identity_does_not_displace_reliable_match():
    row, = match([rec("q", "ACDEFGHIK")],
                 [rec("fragment", "ACDEF"), rec("full", "ACDEYGHIK")])
    assert row["closest_reference_id"] == "full"
    assert row["identity"] == pytest.approx(8 / 9)
    assert row["distance"] == pytest.approx(1 / 9)
    assert row["resolved_reference_count"] == 1


def test_all_distance_ties_are_kept_and_primary_id_is_order_independent():
    query = [rec("q", "ACDEFGHIK")]
    references = [rec("z;id", "ACDEYGHIK"), rec("a", "ACDEWGHIK")]
    forward = match(query, references)
    assert forward == match(query, reversed(references))
    row, = forward
    assert row["closest_reference_id"] == "a"
    assert json.loads(row["closest_reference_ids"]) == ["a", "z;id"]
    assert row["closest_reference_tie_count"] == 2
    assert row["identical_residues"] == 8


def test_equal_sequences_with_distinct_reference_ids_are_retained():
    row, = match([rec("q", "ACDE")], [rec("b", "ACDE"), rec("a", "ACDE")])
    assert json.loads(row["closest_reference_ids"]) == ["a", "b"]
    assert row["distance"] == 0


def test_query_and_reference_namespaces_may_share_ids():
    row, = match([rec("same", "ACDE")], [rec("same", "ACDE")])
    assert row["closest_reference_id"] == "same"


def test_ranking_uses_effective_identity_not_raw_identity():
    row, = match([rec("q", "ACDEFGHIK")],
                 [rec("short", "ACDEFGH"), rec("full", "ACDEYGHIK")], minimum_coverage=0.7)
    assert row["resolved_reference_count"] == 2
    assert row["closest_reference_id"] == "full"
    assert row["distance"] == pytest.approx(1 / 9)


def test_ties_with_different_coverages_use_metrics_from_named_representative():
    row, = match([rec("q", "ACDEFGHIK")],
                 [rec("z_full", "ACDEYGHIK"), rec("a_short", "ACDEFGHI")])
    assert json.loads(row["closest_reference_ids"]) == ["a_short", "z_full"]
    assert row["closest_reference_id"] == "a_short"
    assert row["identity"] == row["target_coverage"] == 1
    assert row["query_coverage"] == pytest.approx(8 / 9)
    assert row["target_length"] == 8


@pytest.mark.parametrize("queries,references", [
    ([], [rec("r", "ACDE")]), ([rec("q", "ACDE")], []),
    ([rec("q", "ACDE"), rec("q", "ACDF")], [rec("r", "ACDE")]),
    ([rec("q", "ACDE")], [rec("r", "ACDE"), rec("r", "ACDF")]),
    ([rec("q", "ACDX")], [rec("r", "ACDE")]),
    ([rec("q", "ACDE")], [rec("r", "")]),
])
def test_invalid_or_empty_pools_fail(queries, references):
    with pytest.raises(ValueError):
        match(queries, references)


def test_cross_product_budget_and_length_guard_are_enforced():
    with pytest.raises(ValueError, match="maximum_pairs"):
        match([rec("q", "ACDE")], [rec("a", "ACDE"), rec("b", "ACDF")], maximum_pairs=1)
    with pytest.raises(ValueError, match="maximum length"):
        match([rec("q", "ACDE")], [rec("r", "ACDEF")], maximum_sequence_length=4)


def test_zero_distance_does_not_override_minimum_alignment_length():
    module = importlib.import_module("gv_eval.nearest_reference")
    row, = module.nearest_reference_matches([rec("q", "ACDE")], [rec("r", "ACDE")])
    assert row["distance_status"] == "unresolved"
    assert row["closest_reference_id"] is None


def test_command_outputs_reproducible_evidence_and_refuses_overwrite(tmp_path):
    root = Path(__file__).resolve().parents[1]
    (tmp_path / "q.fa").write_text(">q\nACDEFGHIK\n>unmatched\nWWWWWWWWW\n", encoding="utf-8")
    (tmp_path / "r.fa").write_text(">a\nACDEFGHIK\n>b\nACDEFGHIK\n", encoding="utf-8")
    config = {
        "candidate_fasta": str(tmp_path / "q.fa"), "reference_fasta": str(tmp_path / "r.fa"),
        "pool_label": "test", "reference_label": "test_reference_unverified",
        "alignment": {"minimum_aligned_residues": 1},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    outputs = [tmp_path / "first", tmp_path / "second"]
    command = [sys.executable, str(root / "experiments/run_nearest_reference.py"), "--config", str(path)]
    for output in outputs:
        run = subprocess.run(command + ["--output-dir", str(output)], capture_output=True, text=True)
        assert run.returncode == 0, run.stderr
        with (output / "nearest_reference.tsv").open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        assert rows[0]["closest_reference_id"] == "a"
        assert json.loads(rows[0]["closest_reference_ids"]) == ["a", "b"]
        assert rows[1]["closest_reference_id"] == rows[1]["distance"] == ""
        assert rows[1]["distance_status"] == "unresolved"
        summary = json.loads((output / "summary.json").read_text())
        assert summary["candidate_count"] == summary["reference_count"] == 2
        assert summary["pair_count"] == 4
        assert summary["matched_candidates"] == summary["unresolved_candidates"] == 1
        manifest = json.loads((output / "manifest.json").read_text())
        for name, expected in manifest["outputs"].items():
            assert sha256_file(output / name) == expected
    for name in ("nearest_reference.tsv", "summary.json", "manifest.json"):
        assert (outputs[0] / name).read_bytes() == (outputs[1] / name).read_bytes()
    before = (outputs[0] / "nearest_reference.tsv").read_bytes()
    run = subprocess.run(command + ["--output-dir", str(outputs[0])], capture_output=True, text=True)
    assert run.returncode != 0
    assert "not empty" in run.stderr
    assert (outputs[0] / "nearest_reference.tsv").read_bytes() == before


@pytest.mark.parametrize("changes", [{"unknown": 1}, {"reference_fasta": ""},
                                    {"alignment": []}, {"pool_label": None}])
def test_invalid_run_config_does_not_create_output(tmp_path, changes):
    module = importlib.import_module("gv_eval.nearest_reference")
    config = {"candidate_fasta": "q.fa", "reference_fasta": "r.fa", "pool_label": "test",
              "reference_label": "test", "alignment": {}} | changes
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(ValueError):
        module.run_nearest_reference(tmp_path, path, "output")
    assert not (tmp_path / "output").exists()


def tiny_run_config(tmp_path):
    (tmp_path / "q.fa").write_text(">q\nACDEFGHIK\n", encoding="utf-8")
    (tmp_path / "r.fa").write_text(">r\nACDEFGHIK\n", encoding="utf-8")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({
        "candidate_fasta": "q.fa", "reference_fasta": "r.fa", "pool_label": "test",
        "reference_label": "test", "alignment": {"minimum_aligned_residues": 1},
    }), encoding="utf-8")
    return path


def test_output_populated_during_alignment_is_not_overwritten(tmp_path, monkeypatch):
    module = importlib.import_module("gv_eval.nearest_reference")
    config = tiny_run_config(tmp_path)
    original = module.nearest_reference_matches
    output = tmp_path / "output"

    def competing_writer(*args, **kwargs):
        rows = original(*args, **kwargs)
        output.mkdir()
        (output / "nearest_reference.tsv").write_text("other run", encoding="utf-8")
        return rows

    monkeypatch.setattr(module, "nearest_reference_matches", competing_writer)
    with pytest.raises(ValueError, match="not empty"):
        module.run_nearest_reference(tmp_path, config, output)
    assert (output / "nearest_reference.tsv").read_text() == "other run"
    assert not (output / "manifest.json").exists()


def test_exclusive_output_claim_rejects_simultaneous_writer_and_releases(tmp_path):
    module = importlib.import_module("gv_eval.nearest_reference")
    output = tmp_path / "output"
    with module._output_claim(output):
        with pytest.raises(ValueError, match="claimed"):
            with module._output_claim(output):
                pytest.fail("Two writers entered the same output directory")
    assert list(output.iterdir()) == []
    with pytest.raises(RuntimeError):
        with module._output_claim(output):
            raise RuntimeError("simulated failure")
    assert list(output.iterdir()) == []


def test_source_change_during_alignment_prevents_publication(tmp_path, monkeypatch):
    module = importlib.import_module("gv_eval.nearest_reference")
    config = tiny_run_config(tmp_path)
    real_hash = module.sha256_file
    original = module.nearest_reference_matches
    changed = False

    def align_then_change(*args, **kwargs):
        nonlocal changed
        rows = original(*args, **kwargs)
        changed = True
        return rows

    def hash_changed_source(path):
        if changed and Path(path).name == "nearest_reference.py":
            return "changed_source_hash"
        return real_hash(path)

    monkeypatch.setattr(module, "nearest_reference_matches", align_then_change)
    monkeypatch.setattr(module, "sha256_file", hash_changed_source)
    with pytest.raises(RuntimeError, match="Source changed"):
        module.run_nearest_reference(tmp_path, config, "output")
    assert not (tmp_path / "output/manifest.json").exists()


def test_development_result_hashes_counts_and_selected_pair_evidence():
    root = Path(__file__).resolve().parents[1]
    output = root / "results/gv02_03_nearest_reference_dev"
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["outputs"].items():
        assert sha256_file(output / name) == expected
    with (output / "nearest_reference.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert len(rows) == summary["candidate_count"] == 200
    assert sum(row["distance_status"] == "resolved" for row in rows) == summary["matched_candidates"] == 165
    assert summary["unresolved_candidates"] == 35
    assert sum(int(row["closest_reference_tie_count"]) > 1 for row in rows) == summary["tied_candidates"] == 55
    queries = {r.identifier: r for r in read_fasta(root / manifest["config"]["candidate_fasta"])}
    references = {r.identifier: r for r in read_fasta(root / manifest["config"]["reference_fasta"])}
    settings = SimilarityConfig(**manifest["effective_alignment_parameters"])
    assert list(queries) == [row["sequence_id"] for row in rows]
    for row in rows:
        ids = json.loads(row["closest_reference_ids"])
        assert len(ids) == int(row["closest_reference_tie_count"])
        assert ids == sorted(set(ids))
        if not ids:
            assert row["closest_reference_id"] == row["distance"] == ""
            continue
        assert row["closest_reference_id"] == ids[0]
        for identifier in ids:
            pair = align_pair(queries[row["sequence_id"]], references[identifier], settings)
            assert pair["distance_status"] == "resolved"
            assert pair["distance"] == pytest.approx(float(row["distance"]))
            if identifier == ids[0]:
                for field in ("identity", "query_coverage", "target_coverage", "alignment_score"):
                    assert pair[field] == pytest.approx(float(row[field]))
