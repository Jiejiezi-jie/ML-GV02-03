from __future__ import annotations

import itertools
import json
import math

import numpy as np
import pandas as pd
import pytest

from gv_eval.pareto_features import (
    AA_COLUMNS, FEATURES, METRICS, analyze_pareto_features, cliffs_delta,
    length_match_sensitivity, match_by_length, read_fasta_strict, sequence_features,
)


def _inputs():
    sequences = {
        "a": "AAADEKRG", "b": "ACDEFGHIKL", "c": "CCCCDEKR",
        "d": "ACDEFGHIK", "e": "ACDEFGHIKLM", "f": "ACDEFGHIKLMN",
    }
    scores = pd.DataFrame({
        "sequence_id": list(sequences),
        "constraint_score": [1, 0.3, 0.2, 0.1, 0.2, 0.05],
        "conservation_score": [0.3, 1, 0.2, 0.1, 0.1, 0.05],
        "novelty_score": [0.5, 0.5, 0.5, 0.4, 0.3, 0.2],
        "uniqueness_score": [0.5, 0.5, 0.5, 0.4, 0.3, 0.2],
        "length": [len(value) for value in sequences.values()],
        "manual_review": [True, False, False, False, False, False],
        "c_manual_review": [True, True, False, False, False, False],
        "uniqueness_resolved_fraction": [0.5, 0.8, 0.7, 0.9, 0.8, 1.0],
    })
    return scores, sequences


def test_known_front_layers_and_tied_percentiles():
    scores, sequences = _inputs()
    result = analyze_pareto_features(scores, sequences)
    candidates = result["candidate_features"].set_index("sequence_id")
    assert candidates["pareto_rank"].to_dict() == {"a": 0, "b": 0, "c": 1, "d": 2, "e": 2, "f": 3}
    assert result["front_details"]["sequence_id"].tolist() == ["a", "b"]
    assert candidates.loc["a", "novelty_score_percentile"] == pytest.approx(100 * 5 / 6)
    assert candidates.loc["c", "novelty_score_percentile"] == pytest.approx(100 * 5 / 6)
    assert candidates.loc["a", "strongest_metrics"] == "constraint_score"
    assert candidates.loc["b", "strongest_metrics"] == "conservation_score"
    assert candidates.loc["c", "strongest_metrics"] == "novelty_score;uniqueness_score"
    assert candidates.loc["a", "weakest_metrics"] == "conservation_score;novelty_score;uniqueness_score"
    summary = result["summary"]
    assert summary["candidate_count"] == 6
    assert summary["pareto_layer_count"] == 4
    assert summary["evidence"]["front"]["manual_review_counts"] == {"manual_review": 1, "c_manual_review": 2}
    json.dumps(summary, allow_nan=False)


def test_sequence_composition_entropy_and_runs_have_known_values():
    uniform = sequence_features("ACDEFGHIKLMNPQRSTVWY")
    assert uniform["shannon_entropy_bits"] == pytest.approx(math.log2(20))
    assert uniform["max_residue_fraction"] == 0.05
    assert uniform["charged_fraction"] == 0.2
    assert uniform["gly_pro_fraction"] == 0.1
    assert uniform["longest_homopolymer"] == 1
    assert sum(uniform[column] for column in AA_COLUMNS) == pytest.approx(1)
    binary = sequence_features("AAAACCCC")
    assert binary["shannon_entropy_bits"] == 1
    assert binary["longest_homopolymer"] == 4
    assert binary["max_residue_fraction"] == 0.5
    assert binary["aa_W"] == 0
    assert sequence_features("AAAGAAAAG")["longest_homopolymer"] == 4
    assert sequence_features("A")["shannon_entropy_bits"] == 0
    assert sequence_features("H")["charged_fraction"] == 0  # DEKR definition excludes H.


def test_cliffs_delta_direction_ties_and_symmetry():
    assert cliffs_delta([1, 2], [2, 3]) == -0.75
    assert cliffs_delta([2, 3], [1, 2]) == 0.75
    assert cliffs_delta([2, 2], [2, 2]) == 0
    assert cliffs_delta([9], [1, 2, 3]) == 1
    with pytest.raises(ValueError, match="nonempty finite"):
        cliffs_delta([], [1])
    with pytest.raises(ValueError, match="nonempty finite"):
        cliffs_delta([np.nan], [1])


def test_length_assignment_matches_exhaustive_optimum_and_avoids_greedy_trap():
    front = pd.DataFrame({"sequence_id": ["a", "b"], "length": [8, 10]})
    controls = pd.DataFrame({"sequence_id": ["x", "y", "z"], "length": [6, 9, 100]})
    result = match_by_length(front, controls)
    oracle = min(sum(abs(left - right) for left, right in zip([8, 10], perm))
                 for perm in itertools.permutations([6, 9, 100], 2))
    assert oracle == 3
    assert result["length_gap"].sum() == oracle
    assert result["control_id"].tolist() == ["x", "y"]
    assert result["control_id"].is_unique
    pd.testing.assert_frame_equal(result, match_by_length(front.iloc[::-1], controls.iloc[::-1]))


def test_tied_length_assignment_is_invariant_to_input_order():
    front = pd.DataFrame({"sequence_id": ["a", "b"], "length": [10, 10]})
    controls = pd.DataFrame({"sequence_id": ["x", "y", "z"], "length": [10, 10, 10]})
    baseline = match_by_length(front, controls)
    for permutation in itertools.permutations(range(3)):
        shuffled = controls.iloc[list(permutation)]
        pd.testing.assert_frame_equal(baseline, match_by_length(front.iloc[::-1], shuffled))
    assert baseline["length_gap"].sum() == 0


def test_descriptive_statistics_pair_differences_and_complete_feature_reporting():
    scores, sequences = _inputs()
    result = analyze_pareto_features(scores, sequences)
    comparisons = result["group_comparison"].set_index(["comparison", "feature"])
    length = comparisons.loc[("front_vs_rest", "length")]
    assert length["front_median"] == 9
    assert length["front_q1"] == 8.5
    assert length["front_q3"] == 9.5
    assert length["control_median"] == 10
    assert length["mean_difference"] == -1
    assert length["cliffs_delta"] == -0.375
    assert len(comparisons) == 2 * (4 + 6 + 20)
    assert len(result["layer_summary"]) == 4 * (4 + 6 + 20)
    assert set(comparisons.index.get_level_values("feature")) == set((*METRICS, *FEATURES, *AA_COLUMNS))
    pairs = result["matched_pairs"]
    assert pairs["length_gap"].sum() == 1
    matched = result["matched_feature_differences"]
    assert len(matched) == 2 * 30
    length_differences = matched.loc[matched["feature"].eq("length"), "difference"]
    assert length_differences.tolist() == [0, 1]
    paired_summary = result["matched_difference_summary"].set_index("feature")
    assert paired_summary.loc["length", "mean"] == 0.5
    assert paired_summary.loc["length", "median"] == 0.5
    assert comparisons.loc[("front_vs_length_matched", "length"), "mean_difference"] == 0.5
    assert result["summary"]["length_matching"]["exact_match_count"] == 1


def test_analysis_is_order_invariant_and_does_not_mutate_input():
    scores, sequences = _inputs()
    before = scores.copy(deep=True)
    original = analyze_pareto_features(scores, sequences)
    shuffled = analyze_pareto_features(scores.iloc[::-1], dict(reversed(list(sequences.items()))))
    for key, value in original.items():
        if isinstance(value, pd.DataFrame):
            pd.testing.assert_frame_equal(value, shuffled[key])
        else:
            assert value == shuffled[key]
    pd.testing.assert_frame_equal(scores, before)


@pytest.mark.parametrize("text, message", [
    (">a\nACD\n>a\nEFG\n", "Duplicate FASTA identifier"),
    (">a\nACD\n>b\nACD\n", "Duplicate FASTA sequence"),
    (">a\nACX\n", "Nonstandard"),
    (">a\nacD\n", "Nonstandard"),
    (">a\nAC D\n", "Nonstandard"),
    (">a\n>b\nDEF\n", "Empty or non-string"),
    (">\nACD\n", "Empty FASTA header"),
    ("ACD\n>a\nDEF\n", "before first header"),
    ("\n", "at least one"),
])
def test_fasta_rejects_ambiguous_records(tmp_path, text, message):
    path = tmp_path / "input.fasta"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_fasta_strict(path)


def test_fasta_accepts_wrapping_and_header_descriptions(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">a description\nACD\nEFG\n\n>b second\nGHIK\n", encoding="utf-8")
    assert read_fasta_strict(path) == {"a": "ACDEFG", "b": "GHIK"}


@pytest.mark.parametrize("column,value,message", [
    ("constraint_score", np.nan, "finite values"),
    ("conservation_score", np.inf, "finite values"),
    ("novelty_score", -0.1, "finite values"),
    ("uniqueness_score", 1.01, "finite values"),
    ("uniqueness_resolved_fraction", 1.01, "finite values"),
    ("length", 99, "length mismatch"),
    ("manual_review", "False", "actual booleans"),
    ("c_manual_review", None, "actual booleans"),
])
def test_analysis_rejects_invalid_metadata(column, value, message):
    scores, sequences = _inputs()
    if isinstance(value, str) or value is None:
        scores[column] = scores[column].astype(object)
    scores.loc[0, column] = value
    with pytest.raises(ValueError, match=message):
        analyze_pareto_features(scores, sequences)


def test_analysis_rejects_identity_and_hash_mismatches():
    scores, sequences = _inputs()
    missing = dict(sequences)
    missing.pop("a")
    with pytest.raises(ValueError, match="ID mismatch"):
        analyze_pareto_features(scores, missing)
    with pytest.raises(ValueError, match="unique sequence IDs"):
        analyze_pareto_features(pd.concat([scores, scores.iloc[:1]]), sequences)
    with pytest.raises(ValueError, match="Duplicate FASTA identifier"):
        analyze_pareto_features(scores, list(sequences.items()) + [("a", "ACD")])
    scores["sequence_sha256"] = "incorrect"
    with pytest.raises(ValueError, match="sequence_sha256 mismatch"):
        analyze_pareto_features(scores, sequences)


def test_all_front_and_insufficient_matching_controls_fail_clearly():
    scores, sequences = _inputs()
    with pytest.raises(ValueError, match="nonempty front and nonfront"):
        analyze_pareto_features(scores.iloc[:2], {key: sequences[key] for key in ["a", "b"]})
    with pytest.raises(ValueError, match="Insufficient nonfront"):
        analyze_pareto_features(scores.iloc[:3], {key: sequences[key] for key in ["a", "b", "c"]})


def test_matching_rejects_overlap_and_invalid_lengths():
    front = pd.DataFrame({"sequence_id": ["a"], "length": [8]})
    with pytest.raises(ValueError, match="disjoint"):
        match_by_length(front, front)
    with pytest.raises(ValueError, match="positive integer"):
        match_by_length(front, pd.DataFrame({"sequence_id": ["b"], "length": [8.5]}))


def test_tie_sensitivity_keeps_exhaustive_primary_optimum_and_unique_controls():
    front = pd.DataFrame({"sequence_id": ["a", "b"], "length": [8, 10], "feature": [0.4, 0.6]})
    control = pd.DataFrame({
        "sequence_id": ["w", "x", "y", "z"], "length": [6, 9, 9, 100],
        "feature": [0.2, 0.3, 0.5, 0.9],
    })
    oracle = min(sum(abs(left - right) for left, right in zip([8, 10], perm))
                 for perm in itertools.permutations([6, 9, 9, 100], 2))
    assert oracle == 2
    result = length_match_sensitivity(front, control, ["length", "feature"], runs=100, seed=42)
    pairs = result["tie_match_pairs"]
    for _, group in pairs.groupby("run"):
        assert set(group["front_id"]) == {"a", "b"}
        assert group["control_id"].is_unique
        assert set(group["control_id"]) == {"x", "y"}
        assert group["length_gap"].sum() == oracle
    assert result["summary"]["baseline_total_gap"] == oracle
    assert result["summary"]["unique_assignments"] == 2
    assert result["summary"]["unique_control_sets"] == 1
    differences = result["tie_match_differences"]
    assert differences.loc[differences["feature"].eq("length"), "mean_difference"].eq(0).all()
    assert np.allclose(differences.loc[differences["feature"].eq("feature"), "mean_difference"], 0.1)


def test_tie_sensitivity_seed_and_input_order_reproducibility():
    front = pd.DataFrame({"sequence_id": ["a"], "length": [10], "feature": [0.5]})
    controls = pd.DataFrame({
        "sequence_id": ["x", "y", "z"], "length": [10, 10, 10], "feature": [0.0, 0.5, 1.0],
    })
    one = length_match_sensitivity(front, controls, ["feature"], seed=42)
    again = length_match_sensitivity(front, controls.iloc[::-1], ["feature"], seed=42)
    different = length_match_sensitivity(front, controls, ["feature"], seed=43)
    for key in ("tie_match_pairs", "tie_match_differences", "tie_match_summary"):
        pd.testing.assert_frame_equal(one[key], again[key])
    assert one["summary"] == again["summary"]
    assert not one["tie_match_pairs"].equals(different["tie_match_pairs"])
    assert one["summary"]["unique_control_sets"] == 3
    summary = one["tie_match_summary"].iloc[0]
    assert summary["min_mean_difference"] == -0.5
    assert summary["max_mean_difference"] == 0.5
    assert summary["negative_count"] + summary["zero_count"] + summary["positive_count"] == 100
    assert min(summary["negative_count"], summary["zero_count"], summary["positive_count"]) > 0
    json.dumps(one["summary"], allow_nan=False)


def test_default_analysis_retains_all_tie_sensitivity_runs_and_features():
    scores, sequences = _inputs()
    result = analyze_pareto_features(scores, sequences)
    assert len(result["tie_match_pairs"]) == 100 * 2
    assert len(result["tie_match_differences"]) == 100 * 30
    assert len(result["tie_match_summary"]) == 30
    assert result["summary"]["length_match_sensitivity"]["runs"] == 100
    assert result["summary"]["length_match_sensitivity"]["seed"] == 42
    assert result["tie_match_pairs"].groupby("run")["length_gap"].sum().eq(1).all()
