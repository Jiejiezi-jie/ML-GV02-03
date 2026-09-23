from __future__ import annotations

import csv
import json
from pathlib import Path

from gv_eval.io import read_fasta, sha256_file


ROOT = Path(__file__).resolve().parents[1]


def _table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_member_a_v1_vae_outputs_are_complete_gated_and_hash_pinned():
    generated = ROOT / "data/generated/sequence_vae/member_a_v1_seed42"
    generation = json.loads((generated / "generation_manifest.json").read_text(encoding="utf-8"))
    assert generation["candidate_count"] == 1000
    assert generation["generation_diagnostics"]["unique_sequence_count"] == 1000
    for name, digest in generation["output_sha256"].items():
        assert sha256_file(generated / name) == digest

    processed = ROOT / "data/processed/gv02_03_v2/vae_member_a_v1"
    qc = json.loads((processed / "qc_summary.json").read_text(encoding="utf-8"))
    pools = json.loads((processed / "candidate_pool_summary.json").read_text(encoding="utf-8"))
    assert qc["input_count"] == qc["pass_count"] == 1000
    assert pools["pool_counts"] == {
        "main_supported_gvpa": 147,
        "ambiguous_exploration": 47,
        "excluded": 806,
    }
    pool_ids = {
        name: {record.identifier for record in read_fasta(processed / f"{name}.fasta")}
        for name in ("main_supported_gvpa", "ambiguous_exploration", "excluded")
    }
    assert len(set.union(*pool_ids.values())) == 1000
    assert not any(pool_ids[left] & pool_ids[right] for left in pool_ids for right in pool_ids if left < right)

    for manifest_name in ("qc_manifest.json", "candidate_pool_manifest.json"):
        manifest = json.loads((processed / manifest_name).read_text(encoding="utf-8"))
        for relative, digest in manifest["inputs"].items():
            assert sha256_file(ROOT / relative) == digest
        for relative, digest in manifest["outputs"].items():
            assert sha256_file(ROOT / relative) == digest

    scores_dir = ROOT / "results/gv02_03_v2/vae_member_a_v1/scoring"
    summary = json.loads((scores_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["candidate_count"] == 147
    assert summary["ranking_candidate_count"] == 136
    scores = _table(scores_dir / "tables/candidate_scores.tsv")
    assert max(float(row["closest_reference_query_coverage"]) for row in scores) <= 1.0
    assert max(float(row["closest_reference_target_coverage"]) for row in scores) <= 1.0
    eligible = {row["sequence_id"] for row in scores if row["domain_pass"] == "True"}
    for strategy in ("weighted_sum", "pareto", "dimension_round_robin"):
        for budget in (10, 20, 50):
            selected = _table(scores_dir / f"selections/{strategy}_top{budget}.tsv")
            assert len(selected) == budget
            assert {row["sequence_id"] for row in selected} <= eligible

    manifest = json.loads((scores_dir / "manifest.json").read_text(encoding="utf-8"))
    for relative, digest in manifest["inputs"].items():
        assert sha256_file(ROOT / relative) == digest
    for relative, digest in manifest["outputs"].items():
        assert sha256_file(ROOT / relative) == digest
