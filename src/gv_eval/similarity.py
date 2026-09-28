"""Auditable global protein alignments with explicit unresolved distances.

This module expects one preselected candidate pool. It does not establish family
membership, and does not turn missing distances into diversity rewards.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from numbers import Real
from pathlib import Path
import platform

import Bio
from Bio.Align import PairwiseAligner, substitution_matrices
import numpy as np
import yaml

from .io import FastaRecord, STANDARD_AA, read_fasta, sha256_file, write_json, write_tsv


@dataclass(frozen=True)
class SimilarityConfig:
    gap_open_score: float = -10.0
    gap_extend_score: float = -0.5
    minimum_coverage: float = 0.80
    minimum_identity: float = 0.20
    minimum_score: float = 0.0
    minimum_aligned_residues: int = 20
    maximum_sequence_length: int = 2000
    maximum_pairs: int = 500000

    def __post_init__(self):
        for name in ("gap_open_score", "gap_extend_score", "minimum_coverage",
                     "minimum_identity", "minimum_score"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
                raise ValueError(f"{name} must be a finite number")
        if not self.gap_open_score <= self.gap_extend_score <= 0:
            raise ValueError("Gap scores must satisfy open <= extend <= 0")
        if not 0 <= self.minimum_coverage <= 1 or not 0 <= self.minimum_identity <= 1:
            raise ValueError("Coverage and identity thresholds must be in [0, 1]")
        for name in ("minimum_aligned_residues", "maximum_sequence_length", "maximum_pairs"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")


def _aligner(config: SimilarityConfig) -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # Penalize internal AND terminal gaps with the same affine scores.
    aligner.open_gap_score = config.gap_open_score
    aligner.extend_gap_score = config.gap_extend_score
    aligner.epsilon = 1e-6
    return aligner


def _validated(records, config):
    rows = list(records)
    if not rows:
        raise ValueError("Candidate pool is empty")
    seen = set()
    normalized = []
    for record in rows:
        if not record.identifier or record.identifier in seen:
            raise ValueError("Candidate identifiers must be nonempty and unique")
        seen.add(record.identifier)
        sequence = record.sequence.upper()
        if not sequence or set(sequence) - STANDARD_AA:
            raise ValueError(f"Invalid protein sequence: {record.identifier}")
        if len(sequence) > config.maximum_sequence_length:
            raise ValueError(f"Sequence exceeds maximum length: {record.identifier}")
        normalized.append(FastaRecord(record.identifier, record.description, sequence))
    return normalized


def _pair(query, target, config, aligner):
    # A fixed sequence orientation makes tied alignments independent of input order.
    left, right = sorted((query.sequence, target.sequence))
    alignment = aligner.align(left, right)[0]
    a, b = alignment[0], alignment[1]
    paired = sum(x != "-" and y != "-" for x, y in zip(a, b))
    identical = sum(x == y and x != "-" for x, y in zip(a, b))
    identity = identical / paired if paired else 0.0
    query_coverage = paired / len(query.sequence)
    target_coverage = paired / len(target.sequence)
    effective_identity = identical / max(len(query.sequence), len(target.sequence))
    reasons = []
    if paired < config.minimum_aligned_residues:
        reasons.append("too_few_aligned_residues")
    if query_coverage < config.minimum_coverage:
        reasons.append("low_query_coverage")
    if target_coverage < config.minimum_coverage:
        reasons.append("low_target_coverage")
    if identity < config.minimum_identity:
        reasons.append("low_identity")
    if alignment.score <= config.minimum_score:
        reasons.append("insufficient_alignment_score")
    return {
        "query_id": query.identifier, "target_id": target.identifier,
        "query_length": len(query.sequence), "target_length": len(target.sequence),
        "alignment_score": float(alignment.score), "alignment_length": len(a),
        "aligned_pairs": paired, "identical_residues": identical,
        "identity": identity, "query_coverage": query_coverage,
        "target_coverage": target_coverage, "effective_identity": effective_identity,
        "distance": None if reasons else 1.0 - effective_identity,
        "distance_status": "unresolved" if reasons else "resolved",
        "distance_reason": ";".join(reasons),
    }


def align_pair(query: FastaRecord, target: FastaRecord, config: SimilarityConfig | None = None) -> dict:
    """Compare full proteins; return evidence even if distance is unresolved."""
    config = config or SimilarityConfig()
    query = _validated([query], config)[0]
    target = _validated([target], config)[0]
    return _pair(query, target, config, _aligner(config))


def build_distance_matrix(records, config: SimilarityConfig | None = None):
    """Preserve input ID order, zero diagonal, NaN unresolved off-diagonal pairs."""
    config = config or SimilarityConfig()
    rows = _validated(records, config)
    count = len(rows)
    if count * (count - 1) // 2 > config.maximum_pairs:
        raise ValueError("Number of pairs exceeds maximum_pairs; split pools or raise the explicit budget")
    matrix = np.full((count, count), np.nan)
    np.fill_diagonal(matrix, 0.0)
    pairs = []
    aligner = _aligner(config)
    for left in range(count):
        for right in range(left + 1, count):
            result = _pair(rows[left], rows[right], config, aligner)
            if result["distance"] is not None:
                matrix[left, right] = matrix[right, left] = result["distance"]
            pairs.append(result)
    return [row.identifier for row in rows], matrix, pairs


PAIR_FIELDS = [
    "query_id", "target_id", "query_length", "target_length", "alignment_score",
    "alignment_length", "aligned_pairs", "identical_residues", "identity",
    "query_coverage", "target_coverage", "effective_identity", "distance",
    "distance_status", "distance_reason",
]


def run_similarity(root: str | Path, config_path: str | Path) -> dict:
    """Run a standalone same-pool analysis; refuse to overwrite an existing run."""
    root = Path(root).resolve()
    config_path = (root / config_path).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    required = {"input_fasta", "output_dir", "pool_label", "alignment"}
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError(f"Configuration must contain exactly {sorted(required)}")
    for name in ("input_fasta", "output_dir", "pool_label"):
        if not isinstance(config[name], str) or not config[name].strip():
            raise ValueError(f"{name} must be a nonempty string")
    if not isinstance(config["alignment"], dict):
        raise ValueError("alignment must be a mapping")
    settings = SimilarityConfig(**config["alignment"])
    fasta = (root / config["input_fasta"]).resolve()
    output = (root / config["output_dir"]).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output directory is not empty; choose a new run directory: {output}")
    inputs = {
        "fasta": {"path": config["input_fasta"], "sha256": sha256_file(fasta)},
        "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
    }
    ids, matrix, pairs = build_distance_matrix(read_fasta(fasta), settings)
    if sha256_file(fasta) != inputs["fasta"]["sha256"] or sha256_file(config_path) != inputs["config"]["sha256"]:
        raise RuntimeError("Input changed during alignment; refusing to write inconsistent results")
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "pairwise_similarity.tsv", pairs, PAIR_FIELDS)
    np.save(output / "candidate_distance.npy", matrix, allow_pickle=False)
    write_json(output / "candidate_distance_ids.json", ids)
    resolved = sum(row["distance_status"] == "resolved" for row in pairs)
    summary = {
        "pool_label": config["pool_label"], "candidate_count": len(ids),
        "pair_count": len(pairs), "resolved_pairs": resolved,
        "unresolved_pairs": len(pairs) - resolved,
        "resolved_fraction": resolved / len(pairs) if pairs else None,
        "missing_distance_policy": "NaN in NPY; blank in TSV; never impute 1 or silently drop",
        "family_verification": "not_performed_by_this_module",
    }
    write_json(output / "summary.json", summary)
    artifacts = ["pairwise_similarity.tsv", "candidate_distance.npy", "candidate_distance_ids.json", "summary.json"]
    write_json(output / "manifest.json", {
        "module_version": "1.0.0", "inputs": inputs,
        "outputs": {name: sha256_file(output / name) for name in artifacts},
        "config": config, "effective_alignment_parameters": asdict(settings),
        "backend": {"mode": "global", "matrix": "BLOSUM62", "epsilon": 1e-6,
                    "algorithm": _aligner(settings).algorithm,
                    "tie_policy": "lexicographic sequence orientation; first optimal alignment"},
        "versions": {"python": platform.python_version(), "biopython": Bio.__version__, "numpy": np.__version__},
        "source_sha256": {name: sha256_file(Path(__file__).with_name(name)) for name in ("similarity.py", "io.py")},
    })
    return summary
