from __future__ import annotations

from gv_eval.io import FastaRecord, read_fasta, write_fasta


def test_fasta_round_trip(tmp_path):
    records = [
        FastaRecord("seq_a", "original source A", "ACDEFGHIK"),
        FastaRecord("seq_b", "seq_b", "LMNPQRSTV"),
    ]
    target = tmp_path / "sequences.fasta"
    write_fasta(records, target, width=4)
    loaded = list(read_fasta(target))
    assert [row.identifier for row in loaded] == ["seq_a", "seq_b"]
    assert [row.sequence for row in loaded] == ["ACDEFGHIK", "LMNPQRSTV"]
    assert loaded[0].description == "seq_a original source A"


def test_sequence_before_header_is_rejected(tmp_path):
    target = tmp_path / "bad.fasta"
    target.write_text("ACDE\n", encoding="utf-8")
    try:
        list(read_fasta(target))
    except ValueError as error:
        assert "sequence before header" in str(error)
    else:
        raise AssertionError("Malformed FASTA was accepted")

