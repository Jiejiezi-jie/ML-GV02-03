#!/usr/bin/env python3
"""Replay frozen-batch QC, global alignments, review consolidation and evaluation."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.c_final_handoff import run as consolidate
from gv_eval.c_handoff import run_handoff
from gv_eval.d_evaluation import FrozenInputs, C_BASE, C_FINAL
from gv_eval.io import sha256_file, write_json
from gv_eval.release import run_release, validate_release_config
from gv_eval.reproducibility import compare_artifacts, verify_hashes


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def reproduce(output: Path) -> dict:
    output = output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be a new or empty directory")
    config_path = ROOT / "configs/d_evaluation.json"
    config = read_json(config_path)
    validate_release_config(config)
    inputs = FrozenInputs(ROOT, {"b": config["b_revision"], "c": config["c_revision"]})
    handoff = ROOT / "data/processed/gv02_03_v2/member_b_handoff"
    release = read_json(handoff / "release_manifest.json")
    split = read_json(handoff / "development_split/split_manifest.json")
    for manifest in (read_json(ROOT / "results/gv02_03_v2/family/manifest.json"), release):
        verify_hashes(ROOT, manifest["inputs"])
        verify_hashes(ROOT, manifest["outputs"])
    verify_hashes(ROOT, split["outputs"])
    if split["release_sha256"] != sha256_file(handoff / "release_manifest.json"):
        raise ValueError("B split is not linked to the A release")
    output.mkdir(parents=True, exist_ok=True)
    c_config_path = ROOT / "configs/c_handoff_member_a_v1.json"
    c_config = read_json(c_config_path)
    snapshot = output / "frozen_b_inputs"
    for spec in c_config["inputs"].values():
        path = (snapshot / spec["path"]).resolve()
        if not path.is_relative_to(snapshot.resolve()):
            raise ValueError("Unsafe frozen input path")
        data = inputs.read("b", spec["path"], spec["sha256"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    comparisons = {}
    first, second = output / "c_run_1", output / "c_run_2"
    for number, destination in enumerate((first, second), 1):
        print(f"C QC/alignment run {number}/2", flush=True)
        run_handoff(snapshot, c_config_path, destination, progress=lambda s: print(s, flush=True))
    c1, c2 = (read_json(p / "qc_similarity_manifest.json") for p in (first, second))
    if c1 != c2:
        raise ValueError("C same-environment runs differ")
    comparisons["c_same_environment"] = compare_artifacts(first, second, c1["outputs"], c2["outputs"])
    published_c = read_json(ROOT / C_BASE / "qc_similarity_manifest.json")
    comparisons["c_published"] = compare_artifacts(ROOT / C_BASE, first, published_c["outputs"], c1["outputs"])
    final = output / "c_review"
    consolidate(first, ROOT / "results/gv02_03_c_domains_member_a_v1",
                ROOT / "results/gv02_03_c_tm_member_a_v1", final)
    fm = read_json(final / "manifest.json")
    old_fm = read_json(ROOT / C_FINAL / "manifest.json")
    comparisons["c_review_published"] = compare_artifacts(ROOT / C_FINAL, final, old_fm["output_sha256"], fm["output_sha256"])
    print("D evaluation and independent repeat", flush=True)
    evaluation = output / "evaluation"
    result = run_release(ROOT, config_path, output=evaluation)
    dm = read_json(evaluation / "manifest.json")
    published_d = ROOT / config["output_dir"]
    old_dm = read_json(published_d / "manifest.json")
    for key in ("source_commits", "config", "config_sha256", "input_sha256"):
        if dm[key] != old_dm[key]:
            raise ValueError(f"Published D input provenance differs: {key}")
    comparisons["d_published"] = compare_artifacts(published_d, evaluation, old_dm["output_sha256"], dm["output_sha256"])
    passed = all(comparison["passed"] for comparison in comparisons.values())
    report = {
        "passed": passed, "status": "passed_frozen_batch_reproduction" if passed else "failed",
        "source_commits": inputs.revisions,
        "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {p.relative_to(ROOT).as_posix(): sha256_file(p) for p in sorted((ROOT / "src/gv_eval").glob("*.py"))},
        "entry_sha256": sha256_file(Path(__file__)),
        "versions": {"python": platform.python_version(), "platform": platform.platform(),
                     **{name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy", "matplotlib", "biopython", "PyYAML")}},
        "reference_release_status": release["status"],
        "split_sizes": {key: value["sequences"] for key, value in split["splits"].items()},
        "evaluation": result, "comparisons": comparisons,
        "replayed": ["B frozen-input/metadata audit", "C baseline/extended QC twice",
                     "C global reference and within-pool alignments twice", "C evidence consolidation",
                     "D evaluation, selections and robustness twice"],
        "verified_existing_only": ["A reference/family evidence", "B generated candidates", "C full-Pfam and PureseqTM results"],
        "not_replayed": ["A external-tool family-model construction", "B original training and sampling", "C full-Pfam scan and PureseqTM predictions"],
        "checkpoint": {"path": "models/generator/sequence_vae/member_a_v1_seed42/best.pt",
                       "available": (ROOT / "models/generator/sequence_vae/member_a_v1_seed42/best.pt").is_file()},
        "scope": "Frozen provisional batch engineering reproduction; existing family/domain/TM evidence is hash-verified, not recomputed. No biological validation.",
    }
    write_json(output / "reproduction_report.json", report)
    write_json(output / "manifest.json", {"output_sha256": {
        p.relative_to(output).as_posix(): sha256_file(p) for p in sorted(output.rglob("*")) if p.is_file()}})
    if not passed:
        raise ValueError("Scientific outputs differ; see reproduction_report.json")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = reproduce(args.output)
    print(json.dumps({"status": report["status"], "ranking_count": report["evaluation"]["ranking_count"],
                      "report": str(args.output / "reproduction_report.json")}, indent=2))


if __name__ == "__main__":
    main()
