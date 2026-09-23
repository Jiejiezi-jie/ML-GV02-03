from dataclasses import replace

import pytest

from gv_eval.classification import FamilyThresholds, classify_family
from gv_eval.external import DomainHit
from gv_eval.io import FastaRecord
from gv_eval.reference import audit_family_sources
from gv_eval.classification import apply_similarity_evidence


RECORD = FastaRecord("query", "untrusted GvpA label", "ACDEFGHIKLMNPQRSTVWY" * 4)


def hit(score=60, model="gvpa", length=80, start=1, end=80):
    return DomainHit("query", 80, model, length, 1e-12, score, 1e-12, score,
                     start, end, start, end)


def classify(hits, pfam=True, **kwargs):
    return classify_family(RECORD, hit(model="Gas_vesicle") if pfam else None,
                           hits, available_profiles={"gvpa", "gvpj", "gvpn"}, **kwargs)


def test_pfam_alone_never_proves_gvpa():
    assert classify({})["family_status"] == "unsupported"


def test_gvpa_positive_and_gvpj_near_negative():
    assert classify({"gvpa": hit(80), "gvpj": hit(30)})["family_status"] == "supported_gvpa"
    row = classify({"gvpa": hit(30), "gvpj": hit(80)})
    assert row["family_status"] == "other_gvp"
    assert row["best_family"] == "gvpj"


@pytest.mark.parametrize("scores", [(50, 50), (50, 49), (49, 50)])
def test_close_competition_is_ambiguous(scores):
    row = classify({"gvpa": hit(scores[0]), "gvpj": hit(scores[1])})
    assert row["family_status"] == "gvpa_gvpj_ambiguous"


def test_other_family_competes_without_pfam_requirement():
    row = classify({"gvpa": hit(30), "gvpn": hit(80)}, pfam=False)
    assert row["family_status"] == "other_gvp"
    assert row["best_family"] == "gvpn"


def test_gvpa_requires_domain_and_coverage():
    assert classify({"gvpa": hit()}, pfam=False)["family_status"] == "unsupported"
    assert classify({"gvpa": hit(90, end=15)})["family_status"] == "unsupported"
    pf = replace(hit(), sequence_score=24)
    assert classify_family(RECORD, pf, {"gvpa":hit()}, available_profiles={"gvpa", "gvpj"})["family_status"] == "unsupported"


def test_missing_profile_is_not_a_negative_search():
    row = classify_family(RECORD, hit(), {"gvpa": hit()}, available_profiles={"gvpa"})
    assert row["classification_reason"] == "missing_competitive_profile"
    assert row["family_status"] == "unsupported"


def test_no_runner_up_has_no_fabricated_margin():
    row = classify({"gvpa": hit()})
    assert row["family_status"] == "supported_gvpa"
    assert row["family_margin"] is None
    assert "no_reportable_runner_up" in row["classification_reason"]


def test_short_competitor_prevents_overconfident_gvpa():
    row = classify({"gvpa":hit(80), "gvpj":hit(90, end=20)})
    assert row["family_status"] != "supported_gvpa"


@pytest.mark.parametrize("seq", ["", "AXX", "ACDEFGHIKLMNPQRSTVWY" * 20])
def test_invalid_or_overlong_cannot_enter_gvpa(seq):
    row = classify_family(FastaRecord("x", "", seq), hit(), {"gvpa":hit()}, available_profiles={"gvpa", "gvpj"})
    assert row["family_status"] != "supported_gvpa"


@pytest.mark.parametrize("kwargs", [{"pfam_coverage":1.1}, {"family_score":float("nan")},
                                    {"margin_bits":0}, {"minimum_gvpa_length":200}])
def test_invalid_thresholds_fail(kwargs):
    with pytest.raises(ValueError):
        FamilyThresholds(**kwargs)


def test_inventory_quarantines_conflicts_and_never_trusts_filename(tmp_path):
    source = tmp_path / "GvpA.fasta"
    source.write_text(
        ">WP_1.1 GvpA\n" + "A"*60 + "\n"
        ">WP_2.1 GvpJ\n" + "A"*60 + "\n"
        ">WP_3.1 unknown protein\n" + "C"*60 + "\n"
        ">WP_4.1 GvpA, partial\n" + "D"*60 + "\n"
        ">WP_5.1 GvpN\n" + "E"*200 + "\n"
        ">WP_6.1 GvpA\n" + "F"*60 + "\n"
        ">WP_6.2 GvpJ\n" + "G"*60 + "\n", encoding="utf-8")
    records, rows = audit_family_sources([source])
    assert len(records) == 6 and len(rows) == 7
    assert [r["seed_eligible"] for r in rows] == [False,False,False,False,True,False,False]
    assert "conflicting_family_annotation" in rows[0]["annotation_reasons"]
    assert rows[4]["annotated_family"] == "gvpn"


def similarity(family, identity, coverage=0.9):
    return dict(family=family, target_id=family + "_ref", effective_identity=identity,
                query_coverage=coverage, target_coverage=coverage)


def test_sequence_competition_vetoes_wrong_hmm_family():
    row = classify({"gvpa":hit(100), "gvpj":hit(50)})
    result = apply_similarity_evidence(row, [similarity("gvpj",0.8), similarity("gvpa",0.5)], FamilyThresholds())
    assert result["hmm_family_status"] == "supported_gvpa"
    assert result["family_status"] == "gvpa_gvpj_ambiguous"
    assert "profile_sequence_family_conflict" in result["classification_reason"]


def test_sequence_agreement_requires_coverage_and_strength():
    row = classify({"gvpa":hit(100), "gvpj":hit(50)})
    for hits in [[], [similarity("gvpa",0.3)], [similarity("gvpa",0.9,0.2)]]:
        assert apply_similarity_evidence(row,hits,FamilyThresholds())["family_status"] == "unsupported"
    assert apply_similarity_evidence(row,[similarity("gvpa",0.8),similarity("gvpj",0.4)],FamilyThresholds())["family_status"] == "supported_gvpa"
    assert apply_similarity_evidence(row,[similarity("gvpa",0.8),similarity("gvpj",0.78)],FamilyThresholds())["family_status"] == "gvpa_gvpj_ambiguous"


def test_similarity_alone_never_promotes_failed_hmm():
    row = classify({})
    assert apply_similarity_evidence(row,[similarity("gvpa",1.0)],FamilyThresholds())["family_status"] == "unsupported"
