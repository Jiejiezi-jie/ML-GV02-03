from __future__ import annotations

from gv_eval.io import FastaRecord
import pytest
from gv_eval.quality import assess_candidates, passing_candidates


def _records(*rows: tuple[str, str]) -> list[FastaRecord]:
    return [FastaRecord(identifier, identifier, sequence) for identifier, sequence in rows]


def _qc(records: list[FastaRecord]) -> dict[str, dict[str, object]]:
    return assess_candidates(records, minimum_length=4, maximum_length=10,
        minimum_unique_residues=3, maximum_single_residue_fraction=0.70,
        generation_length_cap=10)


def test_qc_passes_a_standard_candidate_and_warns_at_generation_cap():
    result = _qc(_records(("normal", "ACDE"), ("at_cap", "ACDEFGHIKL")))
    assert result["normal"]["qc_pass"] is True
    assert result["normal"]["qc_reasons"] == ""
    assert result["at_cap"]["qc_pass"] is True
    assert result["at_cap"]["generation_cap_warning"] is True
    assert result["at_cap"]["qc_reasons"] == "generation_cap_warning"


def test_qc_rejects_nonstandard_and_out_of_range_sequences():
    result = _qc(_records(("bad_residue", "ACDZ"), ("too_short", "ACD"), ("too_long", "ACDEFGHIKLM")))
    assert result["bad_residue"]["qc_pass"] is False
    assert result["bad_residue"]["qc_reasons"] == "nonstandard_amino_acids"
    assert result["too_short"]["qc_reasons"] == "length_below_minimum"
    assert result["too_long"]["qc_reasons"] == "length_above_maximum;generation_cap_warning"


def test_qc_rejects_duplicate_and_low_complexity_candidates():
    result = _qc(_records(("duplicate_a", "ACDE"), ("duplicate_b", "ACDE"), ("low_complexity", "AAAAAC")))
    assert result["duplicate_a"]["duplicate_candidate"] is True
    assert result["duplicate_b"]["duplicate_candidate"] is True
    assert result["duplicate_a"]["qc_reasons"] == "duplicate_candidate"
    assert result["low_complexity"]["low_complexity_warning"] is True
    assert result["low_complexity"]["qc_pass"] is False
    assert result["low_complexity"]["qc_reasons"] == "low_complexity"


def test_passing_candidates_excludes_failed_records_without_reordering():
    records = _records(("pass", "ACDE"), ("fail", "AAAA"), ("later", "ACDF"))
    results = _qc(records)

    assert [record.identifier for record in passing_candidates(records, results)] == ["pass", "later"]


def test_duplicate_detection_is_case_insensitive():
    result = _qc(_records(("one", "acde"), ("two", "ACDE")))
    assert not result["one"]["qc_pass"]
    assert not result["two"]["qc_pass"]


def test_duplicate_identifiers_cannot_overwrite_qc_evidence():
    with pytest.raises(ValueError, match="duplicate identifiers"):
        _qc(_records(("one", "ACDE"), ("one", "ACDF")))


def test_empty_sequence_is_audited_as_failure():
    result = _qc(_records(("empty", "")))["empty"]
    assert not result["qc_pass"]
    assert "empty_sequence" in result["qc_reasons"]


@pytest.mark.parametrize("override", [
    {"minimum_length": 0}, {"maximum_length": 2},
    {"minimum_unique_residues": 0}, {"maximum_single_residue_fraction": float("nan")},
    {"maximum_single_residue_fraction": 1.1}, {"generation_length_cap": 0},
])
def test_invalid_qc_thresholds_fail(override):
    settings = dict(minimum_length=4, maximum_length=10, minimum_unique_residues=3,
                    maximum_single_residue_fraction=0.7)
    with pytest.raises(ValueError):
        assess_candidates(_records(("one", "ACDE")), **(settings | override))
