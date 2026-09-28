"""Join frozen C QC, Pfam and TM evidence without changing B's candidate pools.

The output is a handoff/review index, not a new family or function verdict.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from .c_tm import ensure_empty_output
from .io import sha256_file, write_json, write_tsv


BASE_REQUIRED = {"sequence_id", "sequence_sha256", "length", "pool", "manual_review"}
DOMAIN_REQUIRED = {"sequence_id", "sequence_sha256", "length", "pool", "review_status",
                   "warning_reasons", "pfam_hit_count", "pf00741_hit_segments", "non_target_hit_count",
                   "domain_architecture"}
TM_REQUIRED = {"sequence_id", "sequence_sha256", "length", "pool", "review_status",
               "tm_segment_count", "tm_segments_1based_inclusive"}
ADDED_FIELDS = ("domain_review_status", "domain_warning_reasons", "pf00741_hit_segments",
                "non_target_hit_count", "domain_architecture", "tm_review_status",
                "tm_segment_count", "tm_segments_1based_inclusive", "c_review_reasons",
                "c_manual_review")
DOMAIN_STATUSES = {"manual_review", "no_pfam_hit", "no_extra_domain_evidence"}
TM_STATUSES = {"extra_tm_manual_review", "no_extra_tm_evidence"}


def _manifest(root: Path, filename: str, output_key: str,
              required_output: str) -> tuple[dict, dict[Path, str]]:
    path = root / filename
    manifest = json.loads(path.read_text(encoding="utf-8"))
    expected = manifest[output_key]
    if required_output not in expected:
        raise ValueError(f"Missing required output in {path}: {required_output}")
    verified = {path: sha256_file(path)}
    for name, digest in expected.items():
        relative = Path(name)
        source = (root / relative).resolve()
        if relative.is_absolute() or not source.is_relative_to(root.resolve()):
            raise ValueError(f"Manifest output escapes source directory: {name}")
        if sha256_file(source) != digest:
            raise ValueError(f"Published output hash mismatch: {source}")
        verified[source] = digest
    return manifest, verified


def _read_table(path: Path, required: set[str]) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not required <= set(fields):
            raise ValueError(f"Missing or duplicate TSV fields: {path}")
        rows = list(reader)
    if any(None in row or not row["sequence_id"] for row in rows):
        raise ValueError(f"Malformed TSV row: {path}")
    identifiers = [row["sequence_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate candidate identifier: {path}")
    return rows, fields


def _index(rows: list[dict[str, str]], expected_ids: set[str], label: str) -> dict[str, dict[str, str]]:
    indexed = {row["sequence_id"]: row for row in rows}
    if set(indexed) != expected_ids:
        raise ValueError(f"Candidate identifier set mismatch: {label}")
    return indexed


def _reference_max(root: Path, verified: dict[Path, str], key: str) -> int:
    summary = root / "summary.json"
    if summary not in verified:
        raise ValueError(f"Published summary is not hash-verified: {summary}")
    value = json.loads(summary.read_text(encoding="utf-8"))[key]
    if type(value) is not int or value < 0:
        raise ValueError(f"Invalid reference maximum in {summary}: {key}")
    return value


def _validate_evidence(domain: dict[str, str], tm: dict[str, str], identifier: str,
                       max_target: int, max_tm: int) -> None:
    try:
        hits, target, non_target = (int(domain[key]) for key in
                                    ("pfam_hit_count", "pf00741_hit_segments", "non_target_hit_count"))
    except ValueError as exc:
        raise ValueError(f"Domain evidence/status mismatch: {identifier}") from exc
    warnings = set(filter(None, domain["warning_reasons"].split(";")))
    domain_warning = target > max_target or non_target > 0
    expected_domain_status = ("manual_review" if domain_warning else
                              "no_pfam_hit" if hits == 0 else "no_extra_domain_evidence")
    target_warning = "multiple_pf00741_segments" in warnings
    other_warning = bool(warnings & {"overlapping_alternative_family", "non_target_extra_domain"})
    if (min(hits, target, non_target) < 0 or hits != target + non_target or
            domain["review_status"] != expected_domain_status or
            bool(warnings) != domain_warning or target_warning != (target > max_target) or
            other_warning != (non_target > 0) or
            warnings - {"multiple_pf00741_segments", "overlapping_alternative_family",
                        "non_target_extra_domain"}):
        raise ValueError(f"Domain evidence/status mismatch: {identifier}")

    try:
        count = int(tm["tm_segment_count"])
        spans = [] if not tm["tm_segments_1based_inclusive"] else [
            tuple(map(int, span.split("-")))
            for span in tm["tm_segments_1based_inclusive"].split(";")
        ]
        valid_spans = all(len(span) == 2 and 1 <= span[0] <= span[1] <= int(tm["length"])
                          for span in spans)
    except ValueError as exc:
        raise ValueError(f"TM evidence/status mismatch: {identifier}") from exc
    expected_tm_status = ("extra_tm_manual_review" if count > max_tm else
                          "no_extra_tm_evidence")
    if count < 0 or count != len(spans) or not valid_spans or tm["review_status"] != expected_tm_status:
        raise ValueError(f"TM evidence/status mismatch: {identifier}")


def _reasons(base: dict[str, str], domain: dict[str, str], tm: dict[str, str]) -> list[str]:
    reasons = []
    if base["manual_review"] == "True":
        reasons.append("base:manual_review")
    if domain["review_status"] == "manual_review":
        warnings = [reason for reason in domain["warning_reasons"].split(";") if reason]
        reasons.extend(f"domain:{reason}" for reason in warnings or ["manual_review"])
    elif domain["review_status"] == "no_pfam_hit":
        reasons.append("domain:no_pfam_hit_unresolved")
    if tm["review_status"] == "extra_tm_manual_review":
        reasons.append("tm:extra_tm_manual_review")
    return reasons


def run(base: Path, domains: Path, tm: Path, output: Path) -> dict:
    """Validate three frozen C outputs, then write a warning-only joined handoff."""
    base, domains, tm, output = (Path(path).resolve() for path in (base, domains, tm, output))
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError(f"Output is not empty; existing results will not be overwritten: {output}")
    if any(output == source or output.is_relative_to(source) for source in (base, domains, tm)):
        raise ValueError("Output must be outside the source result directories")

    base_manifest, base_files = _manifest(base, "qc_similarity_manifest.json", "outputs",
                                          "candidate_audit.tsv")
    domain_manifest, domain_files = _manifest(domains, "manifest.json", "output_sha256",
                                              "candidate_domains.tsv")
    tm_manifest, tm_files = _manifest(tm, "manifest.json", "output_sha256", "candidate_tm.tsv")
    verified = base_files | domain_files | tm_files
    max_target = _reference_max(domains, domain_files, "reference_max_pf00741_hit_segments")
    max_tm = _reference_max(tm, tm_files, "reference_max_tm_segments")
    frozen_inputs = {key: base_manifest["inputs"][key]["sha256"]
                     for key in ("candidates", "reference_fasta", "pool_audit")}
    for label, manifest in (("domains", domain_manifest), ("tm", tm_manifest)):
        if manifest["input_sha256"] != frozen_inputs:
            raise ValueError(f"Frozen input hash mismatch: {label}")
        if (label == "domains" or "config_sha256" in manifest) and \
                manifest.get("config_sha256") != base_manifest["config_sha256"]:
            raise ValueError(f"Frozen config hash mismatch: {label}")

    base_rows, base_fields = _read_table(base / "candidate_audit.tsv", BASE_REQUIRED)
    old_review, _ = _read_table(base / "manual_review.tsv", BASE_REQUIRED)
    domain_rows, _ = _read_table(domains / "candidate_domains.tsv", DOMAIN_REQUIRED)
    tm_rows, _ = _read_table(tm / "candidate_tm.tsv", TM_REQUIRED)
    expected_count = base_manifest["config"]["expected_counts"]["candidates"]
    if len(base_rows) != expected_count:
        raise ValueError("Base candidate count mismatch")
    if any(row["manual_review"] not in {"True", "False"} for row in base_rows):
        raise ValueError("Invalid base manual_review value")
    if old_review != [row for row in base_rows if row["manual_review"] == "True"]:
        raise ValueError("Base manual review subset mismatch")
    ids = {row["sequence_id"] for row in base_rows}
    domain_by_id = _index(domain_rows, ids, "domains")
    tm_by_id = _index(tm_rows, ids, "tm")

    merged = []
    for row in base_rows:
        identifier = row["sequence_id"]
        domain = domain_by_id[identifier]
        topology = tm_by_id[identifier]
        for label, other in (("domains", domain), ("tm", topology)):
            if any(row[key] != other[key] for key in
                   ("sequence_sha256", "length", "pool")):
                raise ValueError(f"Candidate identity mismatch: {label}/{identifier}")
        if domain["review_status"] not in DOMAIN_STATUSES:
            raise ValueError(f"Unknown domain review status: {identifier}")
        if topology["review_status"] not in TM_STATUSES:
            raise ValueError(f"Unknown TM review status: {identifier}")
        _validate_evidence(domain, topology, identifier, max_target, max_tm)
        reasons = _reasons(row, domain, topology)
        merged.append({**row,
                       "domain_review_status": domain["review_status"],
                       "domain_warning_reasons": domain["warning_reasons"],
                       "pf00741_hit_segments": domain["pf00741_hit_segments"],
                       "non_target_hit_count": domain["non_target_hit_count"],
                       "domain_architecture": domain["domain_architecture"],
                       "tm_review_status": topology["review_status"],
                       "tm_segment_count": topology["tm_segment_count"],
                       "tm_segments_1based_inclusive": topology["tm_segments_1based_inclusive"],
                       "c_review_reasons": ";".join(reasons),
                       "c_manual_review": bool(reasons)})

    if any(sha256_file(path) != digest for path, digest in verified.items()):
        raise RuntimeError("Source changed during consolidation")
    review = [row for row in merged if row["c_manual_review"]]
    domain_warned = {row["sequence_id"] for row in domain_rows
                     if row["review_status"] == "manual_review"}
    tm_warned = {row["sequence_id"] for row in tm_rows
                 if row["review_status"] == "extra_tm_manual_review"}
    pool_counts = Counter(row["pool"] for row in merged)
    summary = {
        "status": "warning_only_provisional",
        "candidate_count": len(merged),
        "pool_counts": dict(sorted(pool_counts.items())),
        "base_manual_review_count": len(old_review),
        "domain_manual_review_count": len(domain_warned),
        "domain_no_pfam_hit_count": sum(row["review_status"] == "no_pfam_hit"
                                         for row in domain_rows),
        "tm_manual_review_count": len(tm_warned),
        "domain_tm_warning_overlap_count": len(domain_warned & tm_warned),
        "new_warning_id_count": len(domain_warned | tm_warned),
        "c_manual_review_count": len(review),
        "newly_added_review_count": sum(row["manual_review"] == "False"
                                        for row in review),
        "c_manual_review_by_pool": dict(sorted(Counter(row["pool"] for row in review).items())),
        "interpretation": "Review flags are a union of existing C review, Pfam warnings/no-hit uncertainty, and TM warnings; no B pool, QC, family or functional verdict changes.",
    }
    ensure_empty_output(output)
    fields = [*base_fields, *ADDED_FIELDS]
    outputs = {
        "candidate_c_review.tsv": (merged, fields),
        "manual_review_c.tsv": (review, fields),
    }
    for name, (rows, columns) in outputs.items():
        write_tsv(output / name, rows, columns)
    write_json(output / "summary.json", summary)
    write_json(output / "manifest.json", {
        "implementation_sha256_lf": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest(),
        "config_sha256": base_manifest["config_sha256"],
        "frozen_input_sha256": frozen_inputs,
        "source_manifest_sha256": {
            "base": base_files[base / "qc_similarity_manifest.json"],
            "domains": domain_files[domains / "manifest.json"],
            "tm": tm_files[tm / "manifest.json"],
        },
        "output_sha256": {name: sha256_file(output / name)
                          for name in (*outputs, "summary.json")},
    })
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "domains", "tm", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    run(args.base, args.domains, args.tm, args.output)


if __name__ == "__main__":
    main()
