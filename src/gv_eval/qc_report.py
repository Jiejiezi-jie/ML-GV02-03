"""Join frozen QC and nearest-reference evidence without changing eligibility."""
from __future__ import annotations

from collections import Counter
import csv
import json
import math
from pathlib import Path
import re

import yaml

from .io import sha256_file, write_json, write_tsv
from .nearest_reference import MATCH_FIELDS, METRIC_FIELDS, _output_claim


QC_FIELDS = {"sequence_id", "length", "qc_pass", "qc_reasons", "qc_warnings"}
DERIVED_FIELDS = {"qc_failure_reasons", "qc_checks_not_run", "review_required", "review_reasons"}
NOT_ASSESSED = ["eos_metadata", "transmembrane_topology", "domain_architecture", "family_assignment"]


def _boolean(value, field):
    if value not in ("True", "False"):
        raise ValueError(f"{field} must be True or False, got {value!r}")
    return value == "True"


def _integer(value, field, minimum=0):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value) or int(value) < minimum:
        raise ValueError(f"Invalid integer in {field}: {value!r}")
    return int(value)


def _tokens(value, field):
    tokens = value.split(";") if value else []
    if any(not token.strip() or token != token.strip() for token in tokens) or len(set(tokens)) != len(tokens):
        raise ValueError(f"Invalid reason list in {field}")
    return tokens


def _index(rows, required, source):
    if not rows:
        raise ValueError(f"Empty {source} table")
    fields = set(rows[0])
    if not required <= fields or fields & DERIVED_FIELDS:
        raise ValueError(f"Missing required or conflicting derived columns in {source}")
    result = {}
    for row in rows:
        if set(row) != fields or any(not isinstance(value, str) for value in row.values()):
            raise ValueError(f"Inconsistent columns/values in {source}")
        identifier = row["sequence_id"]
        if not identifier.strip() or identifier in result:
            raise ValueError(f"Empty or duplicate sequence_id in {source}: {identifier!r}")
        result[identifier] = row
    return result


def _check_nearest(row):
    status = row["distance_status"]
    if status not in ("resolved", "unresolved"):
        raise ValueError("Unknown distance_status")
    references = _integer(row["reference_count"], "reference_count", 1)
    resolved = _integer(row["resolved_reference_count"], "resolved_reference_count")
    ties = _integer(row["closest_reference_tie_count"], "closest_reference_tie_count")
    if not 0 <= ties <= resolved <= references:
        raise ValueError("Inconsistent reference counts")
    ids = json.loads(row["closest_reference_ids"])
    if not isinstance(ids, list) or any(not isinstance(item, str) or not item.strip() for item in ids):
        raise ValueError("Reference IDs must be a JSON string array")
    if len(ids) != ties or ids != sorted(set(ids)):
        raise ValueError("Inconsistent or unsorted tied reference IDs")
    if status == "unresolved":
        if ids or resolved or row["closest_reference_id"] or not row["distance_reason"]:
            raise ValueError("Unresolved row carries a selected reference or lacks a reason")
        if any(row[field] != "" for field in METRIC_FIELDS):
            raise ValueError("Unresolved row must have blank representative metrics")
    else:
        if not ids or row["closest_reference_id"] != ids[0] or row["distance_reason"]:
            raise ValueError("Resolved representative does not match reference IDs")
        for field in ("target_length", "alignment_length", "aligned_pairs", "identical_residues"):
            _integer(row[field], field, 0 if field == "identical_residues" else 1)
        for field in ("alignment_score", "identity", "query_coverage", "target_coverage", "effective_identity", "distance"):
            value = float(row[field])
            if not math.isfinite(value) or (field != "alignment_score" and not 0 <= value <= 1):
                raise ValueError(f"Invalid finite metric: {field}")
        query_length, target_length = int(row["query_length"]), int(row["target_length"])
        paired, identical = int(row["aligned_pairs"]), int(row["identical_residues"])
        if not 0 <= identical <= paired <= min(query_length, target_length):
            raise ValueError("Inconsistent residue counts or sequence lengths")
        if int(row["alignment_length"]) != query_length + target_length - paired:
            raise ValueError("Inconsistent global alignment length")
        effective_identity = identical / max(query_length, target_length)
        expected = {"identity": identical / paired, "query_coverage": paired / query_length,
                    "target_coverage": paired / target_length, "effective_identity": effective_identity,
                    "distance": 1.0 - effective_identity}
        if any(not math.isclose(float(row[field]), value, rel_tol=1e-9, abs_tol=1e-12)
               for field, value in expected.items()):
            raise ValueError("Inconsistent alignment metric formulas")
    return status, ties


def merge_qc_rows(qc_rows: list[dict], nearest_rows: list[dict]) -> tuple[list[dict], dict]:
    """Merge parsed TSV string rows by ID, retaining QC order and all evidence."""
    left = _index(qc_rows, QC_FIELDS, "QC")
    right = _index(nearest_rows, set(MATCH_FIELDS), "nearest reference")
    if left.keys() != right.keys():
        raise ValueError("Candidate ID sets differ; refusing missing or extra rows")
    if (set(qc_rows[0]) & set(nearest_rows[0])) != {"sequence_id"}:
        raise ValueError("Conflicting source columns; refusing to overwrite evidence")
    merged = []
    failure_counts, warning_counts, unchecked_counts, review_counts = (Counter() for _ in range(4))
    cross = {"pass_resolved": 0, "pass_unresolved": 0, "fail_resolved": 0, "fail_unresolved": 0}
    for identifier, qc in left.items():
        near = right[identifier]
        if _integer(qc["length"], "length", 1) != _integer(near["query_length"], "query_length", 1):
            raise ValueError(f"Candidate length mismatch: {identifier}")
        passed = _boolean(qc["qc_pass"], "qc_pass")
        reasons = _tokens(qc["qc_reasons"], "qc_reasons")
        warnings = _tokens(qc["qc_warnings"], "qc_warnings")
        if not set(warnings) <= set(reasons):
            raise ValueError(f"Warnings missing from legacy qc_reasons: {identifier}")
        failures = [reason for reason in reasons if reason not in warnings]
        if passed != (not failures):
            raise ValueError(f"qc_pass conflicts with hard failure reasons: {identifier}")
        unchecked = sorted(field for field, value in qc.items()
                           if field.endswith("_checked") and not _boolean(value, field))
        status, ties = _check_nearest(near)
        review = []
        for trigger, label in ((not passed, "qc_failure"), (bool(warnings), "qc_warning"),
                               (status == "unresolved", "reference_unresolved"), (ties > 1, "reference_tie")):
            if trigger:
                review.append(label)
        merged.append({**qc, **{key: value for key, value in near.items() if key != "sequence_id"},
                       "qc_failure_reasons": ";".join(failures), "qc_checks_not_run": ";".join(unchecked),
                       "review_required": bool(review), "review_reasons": ";".join(review)})
        failure_counts.update(failures)
        warning_counts.update(warnings)
        unchecked_counts.update(unchecked)
        review_counts.update(review)
        cross[f"{'pass' if passed else 'fail'}_{status}"] += 1
    review_count = sum(row["review_required"] for row in merged)
    summary = {
        "candidate_count": len(merged),
        "qc_pass_count": cross["pass_resolved"] + cross["pass_unresolved"],
        "qc_fail_count": cross["fail_resolved"] + cross["fail_unresolved"],
        "resolved_count": cross["pass_resolved"] + cross["fail_resolved"],
        "unresolved_count": cross["pass_unresolved"] + cross["fail_unresolved"],
        "warning_candidate_count": sum(bool(row["qc_warnings"]) for row in merged),
        "tie_candidate_count": review_counts["reference_tie"],
        "review_count": review_count, "no_review_trigger_count": len(merged) - review_count,
        "cross_counts": cross, "failure_reason_counts": dict(sorted(failure_counts.items())),
        "warning_counts": dict(sorted(warning_counts.items())),
        "not_checked_counts": dict(sorted(unchecked_counts.items())),
        "review_reason_counts": dict(sorted(review_counts.items())),
        "not_assessed_by_report": NOT_ASSESSED,
        "eligibility_policy": "preserve upstream qc_pass; review flags are not exclusions",
    }
    return merged, summary


def _read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields) or any(not name or name.strip() != name for name in fields):
            raise ValueError(f"Invalid or duplicate TSV columns: {path}")
        try:
            rows = list(reader)
        except csv.Error as exc:
            raise ValueError(f"Invalid TSV: {path}") from exc
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"TSV row has missing/extra values: {path}")
    return rows


def _candidate_hash(manifest, key):
    try:
        value = manifest["inputs"][key]["sha256"]
    except (KeyError, TypeError) as exc:
        raise ValueError("Manifest is missing candidate SHA-256") from exc
    if not isinstance(value, str) or not re.fullmatch("[0-9a-fA-F]{64}", value):
        raise ValueError("Invalid candidate SHA-256 in manifest")
    return value.lower()


def _markdown(value):
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("|", "&#124;").replace("\n", " ").replace("\r", " ")


def _report(summary):
    cross = summary["cross_counts"]
    text = ["# C 阶段 QC 与参考匹配汇总", "",
            "本报告汇总已有结果，不重跑比对、不新增淘汰规则，也不推断家族或功能。", "",
            f"共 {summary['candidate_count']} 条候选：QC 通过 {summary['qc_pass_count']} 条、失败 {summary['qc_fail_count']} 条；"
            f"有可靠参考 {summary['resolved_count']} 条、未解析 {summary['unresolved_count']} 条。", "",
            f"人工复核清单共 {summary['review_count']} 条；未触发复核规则 {summary['no_review_trigger_count']} 条。"
            "未触发复核不等于完整 QC 通过。", "",
            "## QC 与参考匹配交叉统计", "", "| QC 状态 | 可靠参考 | 未解析 |", "| --- | ---: | ---: |",
            f"| 通过 | {cross['pass_resolved']} | {cross['pass_unresolved']} |",
            f"| 失败 | {cross['fail_resolved']} | {cross['fail_unresolved']} |", "",
            "## 人工复核口径", "",
            "QC 失败、有任意 QC 警告、参考匹配未解析、最近参考并列：四类取并集，每条候选只列一次。"
            "原因可以重叠，各原因计数不能直接相加。复核标记不会改变 qc_pass。", "",
            f"有警告的候选 {summary['warning_candidate_count']} 条；最近参考并列 {summary['tie_candidate_count']} 条。", ""]
    for title, key in (("硬失败原因", "failure_reason_counts"), ("警告原因", "warning_counts"),
                       ("逐条未检查项", "not_checked_counts")):
        text.extend([f"## {title}", "", "| 项目 | 候选数 |", "| --- | ---: |"])
        text.extend(f"| {_markdown(name)} | {count} |" for name, count in summary[key].items())
        if not summary[key]:
            text.append("| 当前输入无此类记录 | 0 |")
        text.append("")
    text.extend(["## 检查边界与交付文件", "",
                 "EOS 生成元数据、实际跨膜拓扑、非目标/重复结构域、正式家族归属：本汇总均未执行检查。"
                 "疏水/重复模式告警不能替代跨膜或结构域分析。参考未解析不等于高新颖性，可靠参考匹配不等于功能确认。", "",
                 "原始 qc_reasons 包含警告；本报告仅用 qc_reasons 扣除 qc_warnings 后的项目统计硬失败。"
                 "qc_checks_not_run 保留源表 checked=false 的项目；这些项目不会被当作检查通过。", "",
                 "- [逐条总表](candidate_qc_summary.tsv)：保留两份输入全部原始列，新增硬失败原因及复核标记。",
                 "- [人工复核清单](manual_review.tsv)：总表中触发复核规则的完整行，保持原始 QC 顺序。",
                 "- [结构化汇总](summary.json)与[审计记录](manifest.json)。", "",
                 "两份表均已核对对应 manifest 的输出哈希、候选来源哈希、ID 集合和长度；"
                 "哈希核验用于一致性检查，不是来源真实性或生物学正确性的证明。", "",
                 f"候选池标签：{_markdown(summary['source_labels']['pool'])}；"
                 f"参考库标签：{_markdown(summary['source_labels']['reference'])}。", ""])
    return "\n".join(text)


def run_qc_report(root: str | Path, config_path: str | Path, output_dir: str | Path) -> dict:
    root = Path(root).resolve()
    config_path = (root / config_path).resolve()
    config_hash = sha256_file(config_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or set(config) != {"qc_run_dir", "nearest_run_dir"}:
        raise ValueError("Configuration must contain exactly qc_run_dir and nearest_run_dir")
    if any(not isinstance(value, str) or not value.strip() for value in config.values()):
        raise ValueError("Input run directories must be nonempty strings")
    output = (root / output_dir).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output directory is not empty: {output}")
    inputs = {"config": {"path": str(config_path), "sha256": config_hash}}
    manifests, tables = {}, {}
    for name, directory_key, filename in (("qc", "qc_run_dir", "candidate_qc.tsv"),
                                           ("nearest", "nearest_run_dir", "nearest_reference.tsv")):
        directory = (root / config[directory_key]).resolve()
        for label, path in ((f"{name}_manifest", directory / "manifest.json"), (f"{name}_table", directory / filename)):
            inputs[label] = {"path": str(path), "sha256": sha256_file(path)}
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        if not isinstance(manifest, dict) or manifest.get("outputs", {}).get(filename) != inputs[f"{name}_table"]["sha256"]:
            raise ValueError(f"{name} table SHA-256 does not match manifest")
        manifests[name], tables[name] = manifest, _read_tsv(directory / filename)
    candidate_hash = _candidate_hash(manifests["qc"], "candidates")
    if candidate_hash != _candidate_hash(manifests["nearest"], "candidate_fasta"):
        raise ValueError("Candidate source hashes differ between the two runs")
    source_paths = {name: Path(__file__).with_name(name) for name in ("qc_report.py", "nearest_reference.py", "io.py")}
    source_hashes = {name: sha256_file(path) for name, path in source_paths.items()}
    rows, summary = merge_qc_rows(tables["qc"], tables["nearest"])
    labels = manifests["nearest"].get("config", {})
    summary["source_labels"] = {"pool": labels.get("pool_label", "未提供"), "reference": labels.get("reference_label", "未提供")}
    summary["candidate_source_sha256"] = candidate_hash
    report = _report(summary)
    with _output_claim(output):
        if any(sha256_file(item["path"]) != item["sha256"] for item in inputs.values()):
            raise RuntimeError("Input changed while building report")
        if any(sha256_file(path) != source_hashes[name] for name, path in source_paths.items()):
            raise RuntimeError("Source changed while building report")
        fields = list(rows[0])
        write_tsv(output / "candidate_qc_summary.tsv", rows, fields)
        write_tsv(output / "manual_review.tsv", [row for row in rows if row["review_required"]], fields)
        write_json(output / "summary.json", summary)
        (output / "report.md").write_text(report, encoding="utf-8", newline="\n")
        write_json(output / "manifest.json", {
            "module_version": "1.0.0", "inputs": inputs, "source_sha256": source_hashes,
            "candidate_source_sha256": candidate_hash,
            "outputs": {name: sha256_file(output / name) for name in
                        ("candidate_qc_summary.tsv", "manual_review.tsv", "summary.json", "report.md")},
        })
    return summary
