from __future__ import annotations

from gv_eval.reference import build_clean_reference


def test_reference_cleaning_is_label_aware_and_traceable(tmp_path):
    source = tmp_path / "source.fasta"
    source.write_text(
        ">a protein GvpA\n" + "A" * 60 + "\n"
        ">duplicate protein GvpA\n" + "A" * 60 + "\n"
        ">j protein GvpJ\n" + "C" * 60 + "\n"
        ">partial protein GvpA, partial\n" + "D" * 60 + "\n"
        ">unknown protein\n" + "E" * 60 + "\n",
        encoding="utf-8",
    )
    records, mapping, summary = build_clean_reference(
        [source], required_label="gvpa", excluded_labels=["gvpj"],
        excluded_header_terms=["partial", "fragment"], minimum_length=50, maximum_length=180,
    )
    assert len(records) == 1
    assert records[0].identifier == "ref_000001"
    assert len(mapping) == 5
    assert summary["clean_unique_sequences"] == 1
    assert summary["reason_counts"]["exact_duplicate"] == 1
    assert summary["reason_counts"]["excluded_family_label"] == 1

