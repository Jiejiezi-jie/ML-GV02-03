"""Regression tests for the read-only C result consolidation."""
from __future__ import annotations

import csv
import hashlib
from importlib import import_module
import json
from pathlib import Path

import pytest

from gv_eval.io import sequence_sha256, sha256_file, write_json, write_tsv


IDS = ("c1", "c2", "c3", "c4")
SEQUENCES = ("ACDEFGHIKLMNPQRSTVWY" * 2, "FGHIKLMNPQRSTVWYACDE" * 2,
             "KLMNPQRSTVWYACDEFGHI" * 2, "PQRSTVWYACDEFGHIKLMN" * 2)
POOLS = ("main_supported_gvpa", "excluded", "ambiguous_exploration", "main_supported_gvpa")
BASE_FIELDS = ("sequence_id", "sequence_sha256", "length", "pool", "manual_review",
               "qc_pass", "qc_warnings", "qc_failures", "family_status", "pf00741_status",
               "global_distance_status")
DOMAIN_FIELDS = ("sequence_id", "pool", "sequence_sha256", "length", "pfam_hit_count",
                 "pf00741_hit_segments", "non_target_hit_count", "domain_architecture",
                 "warning_reasons", "review_status")
TM_FIELDS = ("sequence_id", "pool", "sequence_sha256", "length", "tm_segment_count",
             "tm_segments_1based_inclusive", "tm_residues", "max_tm_segment_length",
             "review_status")
FROZEN_HASHES = {"candidates": "a" * 64, "reference_fasta": "b" * 64,
                 "pool_audit": "c" * 64}
CONFIG_HASH = "d" * 64


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _run(base: Path, domains: Path, tm: Path, output: Path):
    try:
        module = import_module("gv_eval.c_final_handoff")
    except ModuleNotFoundError:
        pytest.fail("gv_eval.c_final_handoff consolidation is not implemented")
    return module.run(base, domains, tm, output)


@pytest.fixture
def sources(tmp_path: Path) -> tuple[Path, Path, Path]:
    base, domains, tm = (tmp_path / name for name in ("base", "domains", "tm"))
    for folder in (base, domains, tm):
        folder.mkdir()
    base_rows = []
    domain_rows = []
    tm_rows = []
    for index, (identifier, sequence, pool) in enumerate(zip(IDS, SEQUENCES, POOLS, strict=True)):
        common = dict(sequence_id=identifier, sequence_sha256=sequence_sha256(sequence),
                      length=len(sequence), pool=pool)
        base_rows.append({**common, "manual_review": index == 1, "qc_pass": True,
                          "qc_warnings": "", "qc_failures": "", "family_status": "provisional",
                          "pf00741_status": "pass", "global_distance_status": "resolved"})
        domain_rows.append({**common, "pfam_hit_count": 0 if index == 3 else 2 if index == 1 else 1,
                            "pf00741_hit_segments": 0 if index == 3 else 2 if index == 1 else 1,
                            "non_target_hit_count": 0,
                            "domain_architecture": ("" if index == 3 else
                                                    "PF00741:1-20;PF00741:22-40" if index == 1 else
                                                    "PF00741:1-35"),
                            "warning_reasons": "multiple_pf00741_segments" if index == 1 else "",
                            "review_status": ("manual_review" if index == 1 else
                                              "no_pfam_hit" if index == 3 else
                                              "no_extra_domain_evidence")})
        tm_rows.append({**common, "tm_segment_count": 1 if index == 2 else 0,
                        "tm_segments_1based_inclusive": "10-30" if index == 2 else "",
                        "tm_residues": 21 if index == 2 else 0,
                        "max_tm_segment_length": 21 if index == 2 else 0,
                        "review_status": ("extra_tm_manual_review" if index == 2 else
                                          "no_extra_tm_evidence")})
    write_tsv(base / "candidate_audit.tsv", base_rows, BASE_FIELDS)
    write_tsv(base / "manual_review.tsv", [base_rows[1]], BASE_FIELDS)
    write_tsv(domains / "candidate_domains.tsv", domain_rows, DOMAIN_FIELDS)
    write_tsv(tm / "candidate_tm.tsv", tm_rows, TM_FIELDS)
    write_json(domains / "summary.json", {"reference_max_pf00741_hit_segments": 1})
    write_json(tm / "summary.json", {"reference_max_tm_segments": 0})
    write_json(base / "qc_similarity_manifest.json", {
        "config_sha256": CONFIG_HASH,
        "config": {"expected_counts": {"candidates": 4}},
        "inputs": {key: {"sha256": value} for key, value in FROZEN_HASHES.items()},
        "outputs": {name: sha256_file(base / name)
                    for name in ("candidate_audit.tsv", "manual_review.tsv")},
    })
    for folder, filename in ((domains, "candidate_domains.tsv"), (tm, "candidate_tm.tsv")):
        manifest = {"input_sha256": FROZEN_HASHES,
                    "output_sha256": {filename: sha256_file(folder / filename),
                                      "summary.json": sha256_file(folder / "summary.json")}}
        if folder == domains:
            manifest["config_sha256"] = CONFIG_HASH
        write_json(folder / "manifest.json", manifest)
    return base, domains, tm


def test_merges_review_evidence_without_reassigning_pools(sources, tmp_path: Path):
    base, domains, tm = sources
    output = tmp_path / "final"
    summary = _run(base, domains, tm, output)
    rows = _rows(output / "candidate_c_review.tsv")
    review = _rows(output / "manual_review_c.tsv")

    assert [row["sequence_id"] for row in rows] == list(IDS)
    assert [row["pool"] for row in rows] == list(POOLS)
    assert [row["c_manual_review"] for row in rows] == ["False", "True", "True", "True"]
    assert [row["sequence_id"] for row in review] == ["c2", "c3", "c4"]
    assert rows[0]["manual_review"] == "False"
    assert rows[1]["manual_review"] == "True"
    assert rows[1]["domain_review_status"] == "manual_review"
    assert "domain:multiple_pf00741_segments" in rows[1]["c_review_reasons"].split(";")
    assert "tm:extra_tm_manual_review" in rows[2]["c_review_reasons"].split(";")
    assert "domain:no_pfam_hit_unresolved" in rows[3]["c_review_reasons"].split(";")
    assert summary["candidate_count"] == 4
    assert summary["base_manual_review_count"] == 1
    assert summary["domain_manual_review_count"] == 1
    assert summary["domain_no_pfam_hit_count"] == 1
    assert summary["tm_manual_review_count"] == 1
    assert summary["c_manual_review_count"] == 3
    assert summary["newly_added_review_count"] == 2


def test_manifest_records_source_and_output_bytes(sources, tmp_path: Path):
    base, domains, tm = sources
    output = tmp_path / "final"
    _run(base, domains, tm, output)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["config_sha256"] == CONFIG_HASH
    assert manifest["frozen_input_sha256"] == FROZEN_HASHES
    for name, digest in manifest["output_sha256"].items():
        assert sha256_file(output / name) == digest
    assert manifest["source_manifest_sha256"]["base"] == sha256_file(base / "qc_similarity_manifest.json")
    assert manifest["source_manifest_sha256"]["domains"] == sha256_file(domains / "manifest.json")
    assert manifest["source_manifest_sha256"]["tm"] == sha256_file(tm / "manifest.json")
    implementation = Path(import_module("gv_eval.c_final_handoff").__file__)
    expected_code_hash = hashlib.sha256(implementation.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    assert manifest["implementation_sha256_lf"] == expected_code_hash


def test_tampered_published_table_is_rejected_before_output(sources, tmp_path: Path):
    base, domains, tm = sources
    with (domains / "candidate_domains.tsv").open("a", encoding="utf-8") as stream:
        stream.write("tampered\n")
    output = tmp_path / "final"
    with pytest.raises(ValueError, match="hash mismatch"):
        _run(base, domains, tm, output)
    assert not output.exists()


@pytest.mark.parametrize("field,value", [
    ("sequence_sha256", "f" * 64), ("length", "99"), ("pool", "excluded"),
])
def test_identity_mismatch_is_rejected_after_valid_hash(field, value, sources, tmp_path: Path):
    base, domains, tm = sources
    rows = _rows(domains / "candidate_domains.tsv")
    rows[0][field] = value
    write_tsv(domains / "candidate_domains.tsv", rows, DOMAIN_FIELDS)
    manifest_path = domains / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["output_sha256"]["candidate_domains.tsv"] = sha256_file(domains / "candidate_domains.tsv")
    write_json(manifest_path, manifest)
    output = tmp_path / "final"
    with pytest.raises(ValueError, match="identity mismatch"):
        _run(base, domains, tm, output)
    assert not output.exists()


def test_missing_candidate_is_rejected_after_valid_hash(sources, tmp_path: Path):
    base, domains, tm = sources
    rows = _rows(tm / "candidate_tm.tsv")[:-1]
    write_tsv(tm / "candidate_tm.tsv", rows, TM_FIELDS)
    manifest_path = tm / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["output_sha256"]["candidate_tm.tsv"] = sha256_file(tm / "candidate_tm.tsv")
    write_json(manifest_path, manifest)
    output = tmp_path / "final"
    with pytest.raises(ValueError, match="identifier set mismatch"):
        _run(base, domains, tm, output)
    assert not output.exists()


@pytest.mark.parametrize("field,value", [
    ("warning_reasons", "multiple_pf00741_segments"),
    ("pf00741_hit_segments", "2"),
])
def test_domain_status_cannot_hide_warning_evidence(field, value, sources, tmp_path: Path):
    base, domains, tm = sources
    rows = _rows(domains / "candidate_domains.tsv")
    rows[0][field] = value
    if field == "pf00741_hit_segments":
        rows[0]["pfam_hit_count"] = "2"
    write_tsv(domains / "candidate_domains.tsv", rows, DOMAIN_FIELDS)
    manifest_path = domains / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["output_sha256"]["candidate_domains.tsv"] = sha256_file(domains / "candidate_domains.tsv")
    write_json(manifest_path, manifest)
    output = tmp_path / "final"
    with pytest.raises(ValueError, match="Domain evidence/status mismatch"):
        _run(base, domains, tm, output)
    assert not output.exists()


def test_tm_status_cannot_hide_predicted_segment(sources, tmp_path: Path):
    base, domains, tm = sources
    rows = _rows(tm / "candidate_tm.tsv")
    rows[0]["tm_segment_count"] = "1"
    rows[0]["tm_segments_1based_inclusive"] = "10-30"
    write_tsv(tm / "candidate_tm.tsv", rows, TM_FIELDS)
    manifest_path = tm / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["output_sha256"]["candidate_tm.tsv"] = sha256_file(tm / "candidate_tm.tsv")
    write_json(manifest_path, manifest)
    output = tmp_path / "final"
    with pytest.raises(ValueError, match="TM evidence/status mismatch"):
        _run(base, domains, tm, output)
    assert not output.exists()


def test_tm_segment_count_cannot_hide_listed_segment(sources, tmp_path: Path):
    base, domains, tm = sources
    rows = _rows(tm / "candidate_tm.tsv")
    rows[0]["tm_segments_1based_inclusive"] = "10-30"
    write_tsv(tm / "candidate_tm.tsv", rows, TM_FIELDS)
    manifest_path = tm / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["output_sha256"]["candidate_tm.tsv"] = sha256_file(tm / "candidate_tm.tsv")
    write_json(manifest_path, manifest)
    output = tmp_path / "final"
    with pytest.raises(ValueError, match="TM evidence/status mismatch"):
        _run(base, domains, tm, output)
    assert not output.exists()


def test_existing_output_is_not_overwritten(sources, tmp_path: Path):
    base, domains, tm = sources
    output = tmp_path / "final"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("user data", encoding="utf-8")
    with pytest.raises(FileExistsError, match="not empty"):
        _run(base, domains, tm, output)
    assert marker.read_text(encoding="utf-8") == "user data"
