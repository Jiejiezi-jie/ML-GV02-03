#!/usr/bin/env python3
"""Run only candidate QC, without external tools or overwriting prior results."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import platform
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.io import read_fasta, sha256_file, write_fasta, write_json, write_tsv
from gv_eval.quality import assess_candidate_files, passing_candidates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gv02_03_pattern_qc.yaml")
    parser.add_argument("--output-dir", required=True, help="New or empty directory; relative to repo root")
    args = parser.parse_args()
    config_path = (ROOT / args.config).resolve()
    output = (ROOT / args.output_dir).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output directory is not empty: {output}")
    config_hash = sha256_file(config_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    sources = {"candidates": (ROOT / config["inputs"]["candidates"]).resolve()}
    match_paths = {}
    for source in ("training", "reference"):
        value = config["inputs"].get(f"qc_{source}_fasta")
        if value is not None:
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"qc_{source}_fasta must be a nonempty path or null")
            sources[source] = (ROOT / value).resolve()
            match_paths[f"{source}_fasta"] = sources[source]
    input_hashes = {name: sha256_file(path) for name, path in sources.items()}
    records = list(read_fasta(sources["candidates"]))
    if not records:
        raise ValueError("Candidate FASTA is empty")
    qc = assess_candidate_files(records, **match_paths, **config["quality"])
    rows = [{"sequence_id": record.identifier, **qc[record.identifier]} for record in records]
    passing = passing_candidates(records, qc)
    if sha256_file(config_path) != config_hash or any(
        sha256_file(path) != input_hashes[name] for name, path in sources.items()
    ):
        raise ValueError("QC input changed during the run")
    summary = {
        "analysis_scope": "QC only; no family, topology, domain or functional validation",
        "input_count": len(records), "pass_count": len(passing),
        "failure_count": len(records) - len(passing),
        "warning_counts": dict(sorted(Counter(
            warning for row in rows for warning in row["qc_warnings"].split(";") if warning
        ).items())),
        "pattern_checked_counts": {
            name: sum(row[f"{name}_checked"] for row in rows)
            for name in ("hydrophobic", "homopolymer", "tandem_repeat")
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "candidate_qc.tsv", rows, list(rows[0]))
    write_fasta(passing, output / "eligible_candidates.fasta")
    write_json(output / "summary.json", summary)
    manifest = {
        "config": {"path": str(config_path), "sha256": config_hash},
        "quality_options": config["quality"],
        "inputs": {name: {"path": str(path), "sha256": input_hashes[name]}
                   for name, path in sources.items()},
        "outputs": {name: sha256_file(output / name) for name in
                    ("candidate_qc.tsv", "eligible_candidates.fasta", "summary.json")},
        "source_sha256": {name: sha256_file(ROOT / name) for name in
                          ("src/gv_eval/quality.py", "src/gv_eval/sequence_patterns.py",
                           "src/gv_eval/io.py", "experiments/run_quality_audit.py")},
        "python_version": platform.python_version(), "pyyaml_version": yaml.__version__,
    }
    write_json(output / "manifest.json", manifest)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
