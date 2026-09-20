from __future__ import annotations

import pytest

from gv_eval.io import FastaRecord
from gv_eval.quality import assess_candidates


def records(*pairs: tuple[str, str]) -> list[FastaRecord]:
    return [FastaRecord(identifier, identifier, sequence) for identifier, sequence in pairs]


def assess(candidate_sequence: str, **overrides) -> dict[str, object]:
    options = {
        "minimum_length": 4,
        "maximum_length": 100,
        "minimum_unique_residues": 1,
        "maximum_single_residue_fraction": 1.0,
        "reference_records": records(
            ("balanced_a", "ACDEFGHIKLMNPQRSTVWY"),
            ("balanced_b", "YWVTSRQPNMLKIHGFEDCA"),
            ("mixed", "AACDEFGHIKLMNPQRSTVW"),
        ),
        "composition_outlier_quantile": 0.95,
    }
    options.update(overrides)
    return assess_candidates(records(("candidate", candidate_sequence)), **options)["candidate"]


def test_composition_outlier_is_calibrated_from_reference_distribution():
    ordinary = assess("ACDEFGHIKLMNPQRSTVWY")
    outlier = assess("AAAAAAAAAAAAAAAAAAAA")

    assert ordinary["composition_outlier_checked"] is True
    assert ordinary["composition_outlier_warning"] is False
    assert outlier["composition_outlier_warning"] is True
    assert outlier["composition_js_divergence"] > outlier["composition_js_threshold"]
    assert outlier["qc_pass"] is True
    assert "composition_outlier_warning" in outlier["qc_reasons"]


def test_composition_outlier_can_be_a_hard_failure():
    outlier = assess(
        "AAAAAAAAAAAAAAAAAAAA",
        composition_outlier_action="exclude",
    )
    assert outlier["composition_outlier_warning"] is True
    assert outlier["qc_pass"] is False
    assert "composition_outlier" in outlier["qc_reasons"].split(";")


def test_composition_check_is_disabled_by_default():
    row = assess_candidates(
        records(("candidate", "ACDE")),
        minimum_length=4,
        maximum_length=10,
        minimum_unique_residues=1,
        maximum_single_residue_fraction=1.0,
    )["candidate"]
    assert row["composition_outlier_checked"] is False
    assert row["composition_js_divergence"] is None
    assert row["composition_outlier_warning"] is False


@pytest.mark.parametrize(
    "overrides",
    [
        {"composition_outlier_quantile": 0.0},
        {"composition_outlier_quantile": 1.0},
        {"composition_outlier_action": "invalid"},
        {"reference_records": None},
        {"reference_records": records(("only_one", "ACDE"))},
    ],
)
def test_invalid_composition_configuration_fails(overrides):
    with pytest.raises(ValueError):
        assess("ACDE", **overrides)
