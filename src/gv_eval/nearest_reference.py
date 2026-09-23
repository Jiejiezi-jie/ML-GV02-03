"""Nearest reliable reference using the existing global-alignment distance."""
from __future__ import annotations

from dataclasses import asdict
from contextlib import contextmanager
from fractions import Fraction
import json
from pathlib import Path
import platform

import Bio
import yaml

from .io import read_fasta, sha256_file, write_json, write_tsv
from .similarity import PAIR_FIELDS, SimilarityConfig, _aligner, _pair, _validated


METRIC_FIELDS = [field for field in PAIR_FIELDS
                 if field not in ("query_id", "target_id", "query_length", "distance_status", "distance_reason")]
MATCH_FIELDS = ["sequence_id", "query_length", "closest_reference_id", "closest_reference_ids",
                "closest_reference_tie_count", *METRIC_FIELDS, "distance_status", "distance_reason",
                "reference_count", "resolved_reference_count"]


@contextmanager
def _output_claim(output: Path):
    """Serialize this command's publishers, then recheck for prior artifacts."""
    output.mkdir(parents=True, exist_ok=True)
    lock = output / ".nearest_reference.lock"
    try:
        claim = lock.open("x", encoding="utf-8")
    except FileExistsError as exc:
        raise ValueError(f"Output directory is claimed by another run: {output}") from exc
    try:
        with claim:
            if any(path != lock for path in output.iterdir()):
                raise ValueError(f"Output directory is not empty: {output}")
            yield
    finally:
        lock.unlink()


def nearest_reference_matches(candidates, references, config: SimilarityConfig | None = None) -> list[dict]:
    """Exhaustive search; all exact distance ties retained, representative by ID."""
    config = config or SimilarityConfig()
    queries = _validated(candidates, config)
    targets = sorted(_validated(references, config), key=lambda record: record.identifier)
    if len(queries) * len(targets) > config.maximum_pairs:
        raise ValueError("Number of candidate-reference pairs exceeds maximum_pairs")
    aligner = _aligner(config)
    rows = []
    for query in queries:
        best = None
        best_identity = None
        tied_ids = []
        resolved = 0
        for target in targets:
            pair = _pair(query, target, config, aligner)
            if pair["distance_status"] != "resolved":
                continue
            resolved += 1
            # Compare exact rational M/max(Lq,Lr); no floating-point tie tolerance.
            effective_identity = Fraction(pair["identical_residues"], max(pair["query_length"], pair["target_length"]))
            if best is None or effective_identity > best_identity:
                best, best_identity, tied_ids = pair, effective_identity, [target.identifier]
            elif effective_identity == best_identity:
                tied_ids.append(target.identifier)
        rows.append({
            "sequence_id": query.identifier, "query_length": len(query.sequence),
            "closest_reference_id": tied_ids[0] if tied_ids else None,
            "closest_reference_ids": json.dumps(tied_ids, ensure_ascii=False),
            "closest_reference_tie_count": len(tied_ids),
            **{field: best[field] if best is not None else None for field in METRIC_FIELDS},
            "distance_status": "resolved" if best is not None else "unresolved",
            "distance_reason": "" if best is not None else "no_reliable_reference",
            "reference_count": len(targets), "resolved_reference_count": resolved,
        })
    return rows


def run_nearest_reference(root: str | Path, config_path: str | Path, output_dir: str | Path) -> dict:
    """Write one auditable nearest-reference run to a new/empty directory."""
    root = Path(root).resolve()
    config_path = (root / config_path).resolve()
    config_hash = sha256_file(config_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    required = {"candidate_fasta", "reference_fasta", "pool_label", "reference_label", "alignment"}
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError(f"Configuration must contain exactly {sorted(required)}")
    for name in required - {"alignment"}:
        if not isinstance(config[name], str) or not config[name].strip():
            raise ValueError(f"{name} must be a nonempty string")
    if not isinstance(config["alignment"], dict):
        raise ValueError("alignment must be a mapping")
    settings = SimilarityConfig(**config["alignment"])
    output = (root / output_dir).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output directory is not empty: {output}")
    sources = {name: (root / config[name]).resolve() for name in ("candidate_fasta", "reference_fasta")}
    inputs = {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in sources.items()}
    inputs["config"] = {"path": str(config_path), "sha256": config_hash}
    source_paths = {name: Path(__file__).with_name(name)
                    for name in ("nearest_reference.py", "similarity.py", "io.py")}
    source_hashes = {name: sha256_file(path) for name, path in source_paths.items()}
    rows = nearest_reference_matches(read_fasta(sources["candidate_fasta"]),
                                     read_fasta(sources["reference_fasta"]), settings)
    if any(sha256_file(item["path"]) != item["sha256"] for item in inputs.values()):
        raise RuntimeError("Input changed during alignment; refusing inconsistent results")
    if any(sha256_file(path) != source_hashes[name] for name, path in source_paths.items()):
        raise RuntimeError("Source changed during alignment; refusing inconsistent audit")
    matched = sum(row["distance_status"] == "resolved" for row in rows)
    summary = {
        "pool_label": config["pool_label"], "reference_label": config["reference_label"],
        "candidate_count": len(rows), "reference_count": rows[0]["reference_count"],
        "pair_count": len(rows) * rows[0]["reference_count"],
        "matched_candidates": matched, "unresolved_candidates": len(rows) - matched,
        "tied_candidates": sum(row["closest_reference_tie_count"] > 1 for row in rows),
        "exact_match_candidates": sum(row["distance"] == 0 for row in rows),
        "family_verification": "not_performed_by_this_module",
        "selection_policy": "minimum distance among reliable pairs; all exact ties; lexicographic ID representative",
        "missing_distance_policy": "blank ID, metrics and distance when no reliable reference; never impute 1",
    }
    with _output_claim(output):
        write_tsv(output / "nearest_reference.tsv", rows, MATCH_FIELDS)
        write_json(output / "summary.json", summary)
        write_json(output / "manifest.json", {
            "module_version": "1.0.0", "inputs": inputs,
            "outputs": {name: sha256_file(output / name) for name in ("nearest_reference.tsv", "summary.json")},
            "config": config, "effective_alignment_parameters": asdict(settings),
            "backend": {"mode": "global", "matrix": "BLOSUM62", "epsilon": 1e-6,
                        "algorithm": _aligner(settings).algorithm,
                        "alignment_tie_policy": "lexicographic sequence orientation; first optimal alignment",
                        "reference_tie_policy": "exact rational effective identity; representative by lexicographic ID"},
            "versions": {"python": platform.python_version(), "biopython": Bio.__version__, "pyyaml": yaml.__version__},
            "source_sha256": source_hashes,
        })
    return summary
