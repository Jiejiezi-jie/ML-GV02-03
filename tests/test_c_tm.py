from pathlib import Path

import pytest

from gv_eval.c_tm import parse_topology, tm_segments, classify_topology, ensure_empty_output


def test_parse_topology_requires_full_length_binary_prediction(tmp_path: Path):
    top = tmp_path / "one.top"
    top.write_text(">one\nACDEFGHI\n00111100\n", encoding="utf-8")
    assert parse_topology(top, "ACDEFGHI") == "00111100"
    top.write_text(">one\nACDEFGHI\n\n", encoding="utf-8")
    with pytest.raises(ValueError, match="length"):
        parse_topology(top, "ACDEFGHI")
    top.write_text(">one\nACDEFGHI\n0011x100\n", encoding="utf-8")
    with pytest.raises(ValueError, match="binary"):
        parse_topology(top, "ACDEFGHI")


def test_parse_topology_rejects_sequence_mismatch(tmp_path: Path):
    top = tmp_path / "one.top"
    top.write_text(">one\nACDE\n0000\n", encoding="utf-8")
    with pytest.raises(ValueError, match="sequence"):
        parse_topology(top, "ACDF")


def test_tm_segments_are_one_based_inclusive():
    assert tm_segments("0011100110111") == [(3, 5), (8, 9), (11, 13)]
    assert tm_segments("00000") == []
    assert tm_segments("111") == [(1, 3)]


def test_classification_only_warns_when_exceeding_reference_maximum():
    assert classify_topology([(4, 20)], reference_max_segments=1) == "no_extra_tm_evidence"
    assert classify_topology([(4, 20), (30, 48)], reference_max_segments=1) == "extra_tm_manual_review"
    assert classify_topology([], reference_max_segments=0) == "no_extra_tm_evidence"


def test_existing_results_are_not_overwritten(tmp_path: Path):
    destination = tmp_path / "result"
    destination.mkdir()
    (destination / "existing.tsv").write_text("keep", encoding="utf-8")
    with pytest.raises(FileExistsError, match="not empty"):
        ensure_empty_output(destination)
    assert (destination / "existing.tsv").read_text(encoding="utf-8") == "keep"
