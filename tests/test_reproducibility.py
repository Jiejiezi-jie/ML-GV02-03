import json

import numpy as np
import pytest

from gv_eval.io import sha256_file
from gv_eval.reproducibility import compare_artifacts, verify_hashes


def test_portable_comparison_preserves_ids_missing_distances_and_tamper_detection(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir(); b.mkdir()
    (a / "scores.tsv").write_bytes(b"sequence_id\tscore\n001\t0.123456789\n")
    (b / "scores.tsv").write_bytes(b"sequence_id\tscore\r\n001\t0.12345678900001\r\n")
    matrix = np.array([[0., np.nan], [np.nan, 0.]])
    for directory in (a, b):
        np.save(directory / "distance.npy", matrix)
        (directory / "ids.json").write_text(json.dumps(["001", "002"]))
    def hashes(directory):
        return {p.name: sha256_file(p) for p in directory.iterdir()}
    result = compare_artifacts(a, b, hashes(a), hashes(b))
    assert result["passed"]
    assert result["equivalent_data"] == ["scores.tsv"]
    (b / "scores.tsv").write_bytes(b"sequence_id\tscore\n1\t0.123456789\n")
    assert not compare_artifacts(a, b, hashes(a), hashes(b))["passed"]
    (b / "scores.tsv").write_bytes((a / "scores.tsv").read_bytes())
    np.save(b / "distance.npy", np.nan_to_num(matrix, nan=1.))
    assert compare_artifacts(a, b, hashes(a), hashes(b))["failed"] == ["distance.npy"]
    original = hashes(a)
    (a / "scores.tsv").write_text("tampered")
    with pytest.raises(ValueError, match="hash mismatch"):
        compare_artifacts(a, b, original, hashes(b))


def test_missing_artifact_and_path_escape_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="Missing or unsafe"):
        verify_hashes(tmp_path, {"../outside": "irrelevant"})
