import pytest
from gv_eval.adjudication import parse_genbank_proteins, adjudicate_conflict, apply_adjudication
from gv_eval.io import sequence_sha256


def evidence(name="GvpA", source=""):
    return dict(accession="WP_1.1", definition="gas vesicle protein " + name,
        sequence_sha256=sequence_sha256("ACDE"), gene_names=name.lower(),
        name_evidence_category="HMM",name_evidence_accession="example",
        source_identifier=source,region_names="example")


def row():
    return dict(sequence_id="x",sequence_sha256=sequence_sha256("ACDE"),
                source_accessions="WP_1.1",expected_family="gvpj",family_status="supported_gvpa")


def test_official_correction_requires_identical_sequence_and_accession():
    assert adjudicate_conflict(row(),evidence())["biological_identity_resolved"]
    for update in [dict(sequence_sha256="wrong"),dict(accession="WP_2.1")]:
        with pytest.raises(ValueError):
            adjudicate_conflict(row(),dict(evidence(),**update))


def test_pfam_named_gvpj_is_ambiguous_not_relabelled_gvpa():
    d=adjudicate_conflict(row(),evidence("GvpJ","PF00741.24"))
    assert d["operational_family_status"] == "gvpa_gvpj_ambiguous"
    assert not d["biological_identity_resolved"]
    result=apply_adjudication(row(),{d["sequence_sha256"]:d})
    assert result["family_status"] == result["blind_family_status"] == "supported_gvpa"
    assert result["operational_family_status"] == "gvpa_gvpj_ambiguous"
    novel=dict(row(),sequence_sha256=sequence_sha256("ACDF"))
    assert apply_adjudication(novel,{d["sequence_sha256"]:d})["adjudication_disposition"] == "not_reviewed"


def test_genbank_parser_rejects_incomplete_record():
    with pytest.raises(ValueError):
        parse_genbank_proteins("VERSION WP_1.1\nORIGIN\n1 acde\n//")
