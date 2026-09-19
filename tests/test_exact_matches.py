import csv
import json

import pytest

from gv_eval.io import FastaRecord, read_fasta, sha256_file, write_tsv
from gv_eval import quality


def records(*pairs):
    return [FastaRecord(i, i, s) for i, s in pairs]


def assess(**kwargs):
    return quality.assess_candidates(
        records(("q", "ACDE"), ("near", "ACDF")), minimum_length=4,
        maximum_length=10, minimum_unique_residues=3,
        maximum_single_residue_fraction=0.7, **kwargs)


def test_exact_matches_report_all_sorted_ids_and_do_not_match_near_sequences():
    result = assess(training_records=records(("z", "acde"), ("a", "ACDE")),
                    reference_records=records(("ref", "ACDE")))
    assert result["q"]["exact_training_match"] is True
    assert result["q"]["exact_reference_match"] is True
    assert json.loads(result["q"]["training_match_ids"]) == ["a", "z"]
    assert json.loads(result["q"]["reference_match_ids"]) == ["ref"]
    assert result["q"]["qc_pass"] is True
    assert "exact_training_match_warning" in result["q"]["qc_reasons"]
    assert "exact_reference_match_warning" in result["q"]["qc_reasons"]
    assert result["near"]["exact_training_match"] is False
    assert json.loads(result["near"]["training_match_ids"]) == []


@pytest.mark.parametrize("source", ["training", "reference"])
def test_exclusion_is_independent_for_each_source(source):
    result = assess(**{f"{source}_records": records(("match", "ACDE")),
                       f"exact_{source}_match_action": "exclude"})
    assert not result["q"]["qc_pass"]
    assert f"exact_{source}_match" in result["q"]["qc_reasons"].split(";")
    assert result["near"]["qc_pass"]


def test_absent_source_is_not_reported_as_checked():
    row = assess()["q"]
    assert row["training_match_checked"] is False
    assert row["reference_match_checked"] is False


@pytest.mark.parametrize("options", [
    {"exact_training_match_action": "invalid"},
    {"exact_reference_match_action": True},
    {"exact_training_match_action": "exclude"},
    {"training_records": []},
    {"reference_records": records(("r", ""))},
    {"training_records": records(("r", "ACDZ"))},
    {"training_records": records(("r", "ACDE"), ("r", "ACDF"))},
])
def test_invalid_sources_and_policies_fail_explicitly(options):
    with pytest.raises(ValueError):
        assess(**options)


def test_file_inputs_and_tsv_preserve_evidence_and_original_fasta(tmp_path):
    candidate = tmp_path / "candidate.fasta"
    training = tmp_path / "training.fasta"
    reference = tmp_path / "reference.fasta"
    candidate.write_text(">q\nACDE\n", encoding="utf-8")
    training.write_text(">train;1\nacde\n", encoding="utf-8")
    reference.write_text(">ref\nACDF\n", encoding="utf-8")
    before = sha256_file(candidate)
    result = quality.assess_candidate_files(list(read_fasta(candidate)),
        training_fasta=training, reference_fasta=reference,
        minimum_length=4, maximum_length=10, minimum_unique_residues=3,
        maximum_single_residue_fraction=0.7)
    row = {"sequence_id": "q", **result["q"]}
    output = tmp_path / "qc.tsv"
    write_tsv(output, [row], list(row))
    with output.open(encoding="utf-8", newline="") as stream:
        loaded = next(csv.DictReader(stream, delimiter="\t"))
    assert json.loads(loaded["training_match_ids"]) == ["train;1"]
    assert loaded["reference_match_checked"] == "True"
    assert loaded["exact_reference_match"] == "False"
    assert sha256_file(candidate) == before
    with pytest.raises(FileNotFoundError):
        quality.assess_candidate_files([], training_fasta=tmp_path / "missing.fasta")
