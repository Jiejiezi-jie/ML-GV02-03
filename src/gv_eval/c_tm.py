"""Local, warning-only transmembrane topology review for the frozen C handoff.

Runs PureseqTM fast mode in Linux/WSL. No web service or email is used.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
import subprocess

from .io import read_fasta, sequence_sha256, sha256_file, write_json, write_tsv


TOOL_URL = "https://github.com/PureseqTM/PureseqTM_Package"
TOOL_REVISION = "4bb463f5b973f6b7b48bdab55b5cfd150f1741e3"
FIELDS = ("sequence_id", "pool", "sequence_sha256", "length", "tm_segment_count",
          "tm_segments_1based_inclusive", "tm_residues", "max_tm_segment_length", "review_status")


def parse_topology(path: Path, expected_sequence: str) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 3 or not lines[0].startswith(">"):
        raise ValueError(f"Malformed topology file: {path}")
    if lines[1] != expected_sequence:
        raise ValueError(f"Predicted sequence mismatch: {path}")
    if len(lines[2]) != len(expected_sequence):
        raise ValueError(f"Prediction length mismatch: {path}")
    if set(lines[2]) - {"0", "1"}:
        raise ValueError(f"Prediction is not binary: {path}")
    return lines[2]


def tm_segments(topology: str) -> list[tuple[int, int]]:
    """Return contiguous predicted TM spans in 1-based inclusive coordinates."""
    spans: list[tuple[int, int]] = []
    start = None
    for position, label in enumerate(topology + "0", start=1):
        if label == "1" and start is None:
            start = position
        elif label == "0" and start is not None:
            spans.append((start, position - 1))
            start = None
    return spans


def classify_topology(spans: list[tuple[int, int]], reference_max_segments: int) -> str:
    return ("extra_tm_manual_review" if len(spans) > reference_max_segments
            else "no_extra_tm_evidence")


def ensure_empty_output(output: Path) -> None:
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError(f"Output is not empty; existing results will not be overwritten: {output}")
    output.mkdir(parents=True, exist_ok=True)


def _load_records(handoff_root: Path, config: dict):
    records = {}
    for role, key in (("candidate", "candidates"), ("reference", "reference_fasta")):
        entry = config["inputs"][key]
        relative = Path(entry["path"])
        source = (handoff_root / relative).resolve()
        if relative.is_absolute() or not source.is_relative_to(handoff_root.resolve()):
            raise ValueError(f"Input escapes handoff root: {key}")
        if sha256_file(source) != entry["sha256"]:
            raise ValueError(f"Frozen input hash mismatch: {key}")
        values = list(read_fasta(source))
        if len(values) != config["expected_counts"]["candidates" if role == "candidate" else "references"]:
            raise ValueError(f"Frozen input count mismatch: {key}")
        if len({item.identifier for item in values}) != len(values):
            raise ValueError(f"Duplicate FASTA identifier: {key}")
        records[role] = values
    return records


def _load_pools(handoff_root: Path, config: dict) -> dict[str, str]:
    entry = config["inputs"]["pool_audit"]
    path = (handoff_root / entry["path"]).resolve()
    if not path.is_relative_to(handoff_root.resolve()) or sha256_file(path) != entry["sha256"]:
        raise ValueError("Frozen pool audit hash mismatch")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    pools = {row["sequence_id"]: row["pool"] for row in rows}
    if len(pools) != len(rows) or len(rows) != config["expected_counts"]["candidates"]:
        raise ValueError("Invalid pool audit identifiers or count")
    return pools


def _predict_one(tool: Path, cache: Path, role: str, index: int, record) -> tuple[str, str, str]:
    stem = ("c" if role == "candidate" else "r") + f"{index:06d}"
    fasta = cache / "inputs" / f"{stem}.fasta"
    result_dir = cache / "predictions" / stem
    top = result_dir / f"{stem}.top"
    if not top.exists():
        fasta.write_text(f">{stem}\n{record.sequence}\n", encoding="ascii")
        result_dir.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        process = subprocess.run(
            ["bash", str(tool / "PureseqTM.sh"), "-i", str(fasta), "-o", str(result_dir),
             "-m", "0", "-K", "1", "-H", str(tool)],
            cwd=tool, env=env, capture_output=True, text=True,
        )
        if process.returncode:
            raise RuntimeError(f"PureseqTM failed for {role} {record.identifier}: "
                               f"{process.stderr[-1000:]}")
    return role, record.identifier, parse_topology(top, record.sequence)


def _row(record, topology: str, pool: str, reference_max: int) -> dict:
    spans = tm_segments(topology)
    return {
        "sequence_id": record.identifier,
        "pool": pool,
        "sequence_sha256": sequence_sha256(record.sequence),
        "length": len(record.sequence),
        "tm_segment_count": len(spans),
        "tm_segments_1based_inclusive": ";".join(f"{a}-{b}" for a, b in spans),
        "tm_residues": topology.count("1"),
        "max_tm_segment_length": max((b - a + 1 for a, b in spans), default=0),
        "review_status": classify_topology(spans, reference_max),
    }


def run(config_path: Path, handoff_root: Path, tool: Path, cache: Path,
        output: Path, workers: int) -> dict:
    if workers < 1 or workers > 32:
        raise ValueError("workers must be between 1 and 32")
    ensure_empty_output(output)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    records = _load_records(handoff_root, config)
    pools = _load_pools(handoff_root, config)
    if set(pools) != {row.identifier for row in records["candidate"]}:
        raise ValueError("Candidate/pool identifiers disagree")
    for required in (tool / "PureseqTM.sh", tool / "bin" / "DeepCNF_Pred",
                     tool / "param" / "detect.model_fast", tool / "param" / "puretm.model_fast"):
        if not required.is_file():
            raise FileNotFoundError(required)
    cache.joinpath("inputs").mkdir(parents=True, exist_ok=True)

    predictions = {}
    jobs = [(role, index, record) for role in ("reference", "candidate")
            for index, record in enumerate(records[role], start=1)]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_predict_one, tool, cache, role, index, record)
                   for role, index, record in jobs]
        for completed, future in enumerate(as_completed(futures), start=1):
            role, identifier, topology = future.result()
            predictions[(role, identifier)] = topology
            if completed % 100 == 0 or completed == len(jobs):
                print(f"Validated {completed}/{len(jobs)} local topology predictions", flush=True)

    ref_max = max(len(tm_segments(predictions[("reference", row.identifier)]))
                  for row in records["reference"])
    reference_rows = [_row(row, predictions[("reference", row.identifier)], "reference", ref_max)
                      for row in records["reference"]]
    candidate_rows = [_row(row, predictions[("candidate", row.identifier)], pools[row.identifier], ref_max)
                      for row in records["candidate"]]
    write_tsv(output / "reference_tm.tsv", reference_rows, FIELDS)
    write_tsv(output / "candidate_tm.tsv", candidate_rows, FIELDS)
    counts = Counter(row["review_status"] for row in candidate_rows)
    by_pool = {pool: dict(Counter(row["review_status"] for row in candidate_rows if row["pool"] == pool))
               for pool in sorted(set(pools.values()))}
    summary = {
        "method": "PureseqTM fast mode (-m 0), local WSL/Ubuntu, warning-only",
        "tool_url": TOOL_URL, "tool_revision": TOOL_REVISION,
        "reference_tm_segment_histogram": dict(sorted(Counter(row["tm_segment_count"] for row in reference_rows).items())),
        "reference_max_tm_segments": ref_max,
        "candidate_review_status_counts": dict(counts),
        "candidate_review_status_by_pool": by_pool,
        "decision_rule": "Manual review only when predicted TM segment count exceeds maximum of 346 frozen references; no exclusion or pool changes.",
        "limitations": "Predicted 2-state topology is not experimental evidence. Reference set is provisional; absence of a warning does not prove absence of abnormal TM regions.",
    }
    write_json(output / "summary.json", summary)
    manifest = {
        "input_sha256": {key: config["inputs"][key]["sha256"] for key in
                         ("candidates", "reference_fasta", "pool_audit")},
        "tool_url": TOOL_URL, "tool_revision": TOOL_REVISION, "tool_mode": 0,
        "tool_file_sha256": {name: sha256_file(tool / name) for name in
                             ("PureseqTM.sh", "param/detect.model_fast", "param/puretm.model_fast", "param/model_fast.2")},
        "output_sha256": {name: sha256_file(output / name) for name in
                          ("candidate_tm.tsv", "reference_tm.tsv", "summary.json")},
    }
    write_json(output / "manifest.json", manifest)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("config", "handoff-root", "tool", "cache", "output"):
        parser.add_argument(f"--{arg}", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run(args.config.resolve(), args.handoff_root.resolve(), args.tool.resolve(),
        args.cache.resolve(), args.output.resolve(), args.workers)


if __name__ == "__main__":
    main()
