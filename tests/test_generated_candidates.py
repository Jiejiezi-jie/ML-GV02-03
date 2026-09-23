from __future__ import annotations

import csv
import json
from pathlib import Path

import yaml

from gv_eval.generated_candidates import build_generated_candidate_pools, prepare_generated_qc
from gv_eval.io import read_fasta, sequence_sha256


def _write_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_generated_qc_and_family_pools_are_traceable(tmp_path):
    sequences = {"a": "ACDEFGHIK", "b": "LMNPQRSTV", "c": "ACDFGHIKL"}
    (tmp_path / "candidates.fasta").write_text(
        "".join(f">{key}\n{sequence}\n" for key, sequence in sequences.items()), encoding="utf-8"
    )
    (tmp_path / "train.fasta").write_text(">t\nACDEFGHIL\n", encoding="utf-8")
    (tmp_path / "reference.fasta").write_text(">r1\nACDEFGHIM\n>r2\nLMNPQRSTW\n", encoding="utf-8")
    metadata = [
        {
            "sequence_id": key,
            "sequence_length": len(sequence),
            "sequence_sha256": sequence_sha256(sequence),
            "terminated_by_eos": str(key != "c").lower(),
            "hit_generation_cap": str(key == "c").lower(),
            "stop_reason": "length_cap" if key == "c" else "eos",
        }
        for key, sequence in sequences.items()
    ]
    _write_table(tmp_path / "metadata.tsv", metadata)
    config = {
        "version": 1,
        "seed": 42,
        "inputs": {
            "candidates": "candidates.fasta",
            "generation_metadata": "metadata.tsv",
            "training_fasta": "train.fasta",
            "reference_fasta": "reference.fasta",
            "family_classification": "processed/family/candidate_family_classification.tsv",
        },
        "quality": {
            "minimum_length": 8,
            "maximum_length": 10,
            "minimum_unique_residues": 4,
            "maximum_single_residue_fraction": 0.5,
            "generation_length_cap": None,
            "exact_training_match_action": "exclude",
            "exact_reference_match_action": "exclude",
            "composition_outlier_quantile": 0.99,
            "composition_outlier_action": "warn",
        },
        "outputs": {"processed_dir": "processed"},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    qc = prepare_generated_qc(tmp_path, config_path)
    assert qc["pass_count"] == 3
    assert qc["generation_cap_count"] == 1

    family_rows = []
    statuses = {"a": ("supported_gvpa", "pass"), "b": ("gvpa_gvpj_ambiguous", "pass"), "c": ("unsupported", "fail")}
    for key, sequence in sequences.items():
        status, pfam = statuses[key]
        family_rows.append({"sequence_id": key, "sequence_sha256": sequence_sha256(sequence), "family_status": status, "pf00741_status": pfam})
    _write_table(tmp_path / "processed/family/candidate_family_classification.tsv", family_rows)
    pools = build_generated_candidate_pools(tmp_path, config_path)
    assert pools["pool_counts"] == {"main_supported_gvpa": 1, "ambiguous_exploration": 1, "excluded": 1}
    assert [record.identifier for record in read_fasta(tmp_path / "processed/main_supported_gvpa.fasta")] == ["a"]
    manifest = json.loads((tmp_path / "processed/candidate_pool_manifest.json").read_text(encoding="utf-8"))
    assert manifest["outputs"]
