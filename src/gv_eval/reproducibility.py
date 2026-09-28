"""Compare frozen scientific outputs across environments without hiding drift."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from .io import sha256_file


def verify_hashes(directory: Path, hashes: dict[str, str]) -> None:
    base = directory.resolve()
    for name, digest in hashes.items():
        path = (base / name).resolve()
        if not path.is_relative_to(base) or not path.is_file():
            raise ValueError(f"Missing or unsafe artifact: {name}")
        if sha256_file(path) != digest:
            raise ValueError(f"Artifact hash mismatch: {name}")


def _json_equal(left, right, rtol: float, atol: float) -> bool:
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(
            _json_equal(left[key], right[key], rtol, atol) for key in left
        )
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(
            _json_equal(a, b, rtol, atol) for a, b in zip(left, right)
        )
    if type(left) in (int, float) and type(right) in (int, float):
        if type(left) is type(right) is int:
            return left == right
        return math.isclose(left, right, rel_tol=rtol, abs_tol=atol)
    return type(left) is type(right) and left == right


def _table_equal(left: Path, right: Path, rtol: float, atol: float) -> bool:
    def read(path):
        with path.open(encoding="utf-8-sig", newline="") as stream:
            return list(csv.reader(stream, delimiter="\t"))
    a, b = read(left), read(right)
    if len(a) != len(b) or not a or a[0] != b[0]:
        return a == b
    for row_a, row_b in zip(a[1:], b[1:]):
        if len(row_a) != len(row_b) or len(row_a) != len(a[0]):
            return False
        for column, x, y in zip(a[0], row_a, row_b):
            if x == y:
                continue
            # IDs, labels, hashes and integer counts must never be rounded away.
            if (column.endswith(("id", "ids")) or "sha256" in column
                    or not any(c in x + y for c in ".eE")):
                return False
            try:
                if not math.isclose(float(x), float(y), rel_tol=rtol, abs_tol=atol):
                    return False
            except ValueError:
                return False
    return True


def compare_artifacts(reference: Path, repeated: Path,
                      reference_hashes: dict[str, str], repeated_hashes: dict[str, str],
                      *, rtol: float = 1e-9, atol: float = 1e-12) -> dict:
    """Verify bytes first, then compare data; report image byte changes separately.

    Cross-environment tolerance applies only to floating point data. Table row
    order, IDs, integer counts, JSON structure and missing matrix entries remain
    fixed. PNGs are never used to establish numeric reproducibility.
    """
    if not math.isfinite(rtol) or not math.isfinite(atol) or min(rtol, atol) < 0:
        raise ValueError("Comparison tolerances must be finite and nonnegative")
    verify_hashes(reference, reference_hashes)
    verify_hashes(repeated, repeated_hashes)
    if reference_hashes.keys() != repeated_hashes.keys():
        raise ValueError("Artifact sets differ")
    identical, equivalent, rendering, failed = [], [], [], []
    for name, digest in reference_hashes.items():
        if digest == repeated_hashes[name]:
            identical.append(name)
            continue
        left, right = reference / name, repeated / name
        suffix = left.suffix.lower()
        if suffix == ".png":
            rendering.append(name)
            continue
        if suffix == ".npy":
            a, b = np.load(left, allow_pickle=False), np.load(right, allow_pickle=False)
            equal = a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
        elif suffix == ".tsv":
            equal = _table_equal(left, right, rtol, atol)
        elif suffix == ".json":
            equal = _json_equal(json.loads(left.read_text(encoding="utf-8")),
                                json.loads(right.read_text(encoding="utf-8")), rtol, atol)
        elif suffix in (".fasta", ".faa", ".md"):
            equal = left.read_bytes().replace(b"\r\n", b"\n") == right.read_bytes().replace(b"\r\n", b"\n")
        else:
            equal = False
        (equivalent if equal else failed).append(name)
    return {
        "passed": not failed, "compared_artifact_count": len(reference_hashes),
        "byte_identical_count": len(identical), "byte_identical": identical,
        "equivalent_data": equivalent, "rendering_byte_differences": rendering,
        "failed": failed, "rtol": rtol, "atol": atol,
        "scope": "Scientific data and ordered selections; PNG differences are recorded, not treated as numerical evidence.",
    }
