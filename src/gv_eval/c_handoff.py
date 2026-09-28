"""Read-only B handoff audit and standalone C QC/global-similarity delivery."""
from __future__ import annotations

from collections import Counter
import csv
from dataclasses import asdict
import json
from pathlib import Path
import platform

import Bio
import numpy as np
import yaml

from .io import STANDARD_AA, read_fasta, sequence_sha256, sha256_file, write_fasta, write_json, write_tsv
from .nearest_reference import MATCH_FIELDS, _output_claim, nearest_reference_matches
from .quality import assess_candidates
from .similarity import PAIR_FIELDS, SimilarityConfig, build_distance_matrix

POOLS = ("main_supported_gvpa", "ambiguous_exploration", "excluded")
PATTERNS = {"hydrophobic_run_minimum", "homopolymer_minimum", "tandem_repeat_minimum_length",
            "tandem_repeat_minimum_copies", "tandem_repeat_maximum_motif_length"}
REQUIRED_INPUTS = {"candidates", "generation_metadata", "generation_manifest", "training_fasta",
                   "reference_fasta", "b_config", "candidate_qc", "family_classification", "pool_audit", *POOLS}


def _index(rows, key="sequence_id"):
    result = {}
    for row in rows:
        identifier = row[key]
        if not identifier or identifier in result:
            raise ValueError(f"Empty or duplicate {key}: {identifier}")
        result[identifier] = row
    return result


def _table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Invalid TSV headers: {path}")
        rows = list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError(f"Malformed TSV row: {path}")
        return rows


def _bool(value):
    if value not in ("true", "false", "True", "False"):
        raise ValueError(f"Invalid boolean: {value!r}")
    return value in ("true", "True")


def _int(value):
    if isinstance(value, bool) or not str(value).isdigit():
        raise ValueError(f"Invalid nonnegative integer: {value!r}")
    return int(value)


def _tokens(value):
    return set(filter(None, value.split(";")))


def audit_generation(records, metadata, manifest):
    """Validate exported metadata consistency, not checkpoint replay/raw tokens."""
    ids = [r.identifier for r in records]
    if not ids or len(set(ids)) != len(ids) or any(not i for i in ids):
        raise ValueError("Empty/duplicate candidate identifiers")
    rows = _index(metadata)
    if set(rows) != set(ids) or _int(manifest["candidate_count"]) != len(ids):
        raise ValueError("Generation candidate/metadata ID or count mismatch")
    parameters = manifest["generation_parameters"]
    minimum, maximum = (_int(parameters[k]) for k in ("minimum_sequence_length", "maximum_sequence_length"))
    if minimum < 1 or maximum < minimum:
        raise ValueError("Invalid generation length limits")
    result = {}
    for record in records:
        row, length = rows[record.identifier], len(record.sequence)
        if not record.sequence or set(record.sequence) - STANDARD_AA:
            raise ValueError(f"Invalid candidate sequence: {record.identifier}")
        if row["sequence_sha256"] != sequence_sha256(record.sequence) or _int(row["sequence_length"]) != length:
            raise ValueError(f"Generation sequence hash/length mismatch: {record.identifier}")
        eos, cap = _bool(row["terminated_by_eos"]), _bool(row["hit_generation_cap"])
        if eos == cap or row["stop_reason"] != ("eos" if eos else "length_cap"):
            raise ValueError(f"Inconsistent EOS/cap flags: {record.identifier}")
        if not minimum <= length <= maximum or (cap and length != maximum):
            raise ValueError(f"Generation length/cap mismatch: {record.identifier}")
        if _int(row["raw_token_length"]) != length + int(eos):
            raise ValueError(f"Raw token length mismatch: {record.identifier}")
        if "data_provenance" in manifest:
            provenance = manifest["data_provenance"]
            for field, expected in (("reference_release_id", provenance["release_id"]),
                                    ("split_manifest_sha256", provenance["split_manifest_sha256"]),
                                    ("generation_seed", str(manifest["seed"]))):
                if row[field] != expected:
                    raise ValueError(f"Generation provenance mismatch: {field}")
            for field, suffix in (("checkpoint_sha256", ".pt"), ("vocabulary_sha256", "vocabulary.json")):
                hashes = {v for k, v in manifest["input_sha256"].items() if k.endswith(suffix)}
                if len(hashes) != 1 or row[field] not in hashes:
                    raise ValueError(f"Generation provenance mismatch: {field}")
        result[record.identifier] = {"sequence_sha256": row["sequence_sha256"], "terminated_by_eos": eos,
                                     "hit_generation_cap": cap, "stop_reason": row["stop_reason"],
                                     "raw_token_length": _int(row["raw_token_length"])}
    for key in ("terminated_by_eos", "hit_generation_cap"):
        if sum(row[key] for row in result.values()) != _int(manifest[key]):
            raise ValueError(f"Generation manifest {key} count mismatch")
    return result


def _load_inputs(root, config):
    if set(config) != {"inputs", "expected_counts", "patterns", "alignment"}:
        raise ValueError("Invalid C handoff configuration keys")
    if not REQUIRED_INPUTS <= config["inputs"].keys():
        raise ValueError("Missing required handoff inputs")
    if set(config["patterns"]) - PATTERNS:
        raise ValueError("Only warning-pattern options may override B QC")
    paths = {}
    for key, entry in config["inputs"].items():
        relative = Path(entry["path"])
        path = (root / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(root):
            raise ValueError(f"Input path escapes handoff root: {key}")
        if sha256_file(path) != entry["sha256"]:
            raise ValueError(f"Frozen input hash mismatch: {key}")
        paths[key] = path
    return paths


def _quality_and_pools(records, paths, metadata, quality, patterns):
    baseline = assess_candidates(records, training_records=read_fasta(paths["training_fasta"]),
                                  reference_records=read_fasta(paths["reference_fasta"]), **quality)
    extended = assess_candidates(records, training_records=read_fasta(paths["training_fasta"]),
                                  reference_records=read_fasta(paths["reference_fasta"]), **quality, **patterns)
    frozen, family, pools = (_index(_table(paths[k])) for k in ("candidate_qc", "family_classification", "pool_audit"))
    ids = set(baseline)
    if any(set(table) != ids for table in (frozen, family, pools)):
        raise ValueError("QC/family/pool table ID mismatch")
    rows = []
    for record in records:
        identifier, digest = record.identifier, sequence_sha256(record.sequence)
        base, extra, old, fam, meta = baseline[identifier], extended[identifier], frozen[identifier], family[identifier], metadata[identifier]
        reasons = _tokens(base["qc_reasons"])
        if meta["hit_generation_cap"]: reasons.add("generation_cap_warning")
        if old["sequence_sha256"] != digest or _bool(old["qc_pass"]) != base["qc_pass"] or _tokens(old["qc_reasons"]) != reasons:
            raise ValueError(f"B baseline QC mismatch: {identifier}")
        if fam["sequence_sha256"] != digest or _int(fam["sequence_length"]) != len(record.sequence):
            raise ValueError(f"Family sequence hash/length mismatch: {identifier}")
        expected_pool = "excluded"
        if base["qc_pass"]:
            if fam["family_status"] == "supported_gvpa" and fam["pf00741_status"] == "pass":
                expected_pool = "main_supported_gvpa"
            elif fam["family_status"] == "gvpa_gvpj_ambiguous":
                expected_pool = "ambiguous_exploration"
        if pools[identifier]["pool"] != expected_pool:
            raise ValueError(f"B pool classification mismatch: {identifier}")
        warnings = _tokens(extra["qc_warnings"])
        failures = _tokens(extra["qc_reasons"]) - warnings
        if meta["hit_generation_cap"]: warnings.add("generation_cap_warning")
        extra.update(generation_cap_warning=meta["hit_generation_cap"], qc_warnings=";".join(sorted(warnings)),
                     qc_failures=";".join(sorted(failures)), qc_reasons=";".join(sorted(failures | warnings)))
        rows.append(dict(sequence_id=identifier, **meta, **extra, pool=expected_pool,
                         family_status=fam["family_status"], pf00741_status=fam["pf00741_status"]))
    pool_records = {}
    for pool in POOLS:
        members = list(read_fasta(paths[pool]))
        actual = {r.identifier: r.sequence for r in members}
        expected = {r.identifier: r.sequence for r in records if pools[r.identifier]["pool"] == pool}
        if len(actual) != len(members) or actual != expected:
            raise ValueError(f"Pool FASTA membership/sequence mismatch: {pool}")
        pool_records[pool] = members
    return rows, pool_records


def run_handoff(handoff_root, config_path, output_dir, progress=None):
    """Run into a new directory; never mutate handoff files or reassign pools."""
    root, config_path, output = (Path(p).resolve() for p in (handoff_root, config_path, output_dir))
    if output == root or output.is_relative_to(root):
        raise ValueError("Output must be outside the read-only handoff root")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError(f"Output directory is not empty: {output}")
    log = progress or (lambda message: None)
    config_hash = sha256_file(config_path)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    paths = _load_inputs(root, config)
    settings = SimilarityConfig(**config["alignment"])
    source_paths = {name: Path(__file__).with_name(name) for name in
                    ("c_handoff.py", "quality.py", "sequence_patterns.py", "nearest_reference.py", "similarity.py", "io.py")}
    source_hashes = {name: sha256_file(path) for name, path in source_paths.items()}
    generation = json.loads(paths["generation_manifest"].read_text(encoding="utf-8"))
    for key in ("candidates", "generation_metadata"):
        if generation["output_sha256"].get(paths[key].name) != config["inputs"][key]["sha256"]:
            raise ValueError(f"Generation manifest output mismatch: {key}")
    records = list(read_fasta(paths["candidates"]))
    references = list(read_fasta(paths["reference_fasta"]))
    metadata = audit_generation(records, _table(paths["generation_metadata"]), generation)
    quality = yaml.safe_load(paths["b_config"].read_text(encoding="utf-8"))["quality"]
    if set(quality) & PATTERNS or quality.get("generation_length_cap") is not None:
        raise ValueError("B QC must use explicit metadata cap and no C pattern overrides")
    rows, pools = _quality_and_pools(records, paths, metadata, quality, config["patterns"])
    counts = dict(candidates=len(records), references=len(references), **{k: len(v) for k, v in pools.items()})
    if counts != config["expected_counts"]:
        raise ValueError(f"Handoff count mismatch: {counts}")
    log(f"Input audit passed: {counts}; baseline QC reproduced")
    matches = []
    # Keep the total budget guard even though progress is emitted per chunk.
    if len(records) * len(references) > settings.maximum_pairs:
        raise ValueError("Candidate-reference pairs exceed maximum_pairs")
    for start in range(0, len(records), 50):
        matches.extend(nearest_reference_matches(records[start:start + 50], references, settings))
        log(f"Global nearest reference: {len(matches)}/{len(records)}")
    matrices, pair_summaries = {}, {}
    for pool, short in zip(POOLS[:2], ("main", "ambiguous")):
        if pools[pool]:
            ids, matrix, pairs = build_distance_matrix(pools[pool], settings)
        else:
            ids, matrix, pairs = [], np.empty((0, 0)), []
        matrices[short] = (ids, matrix, pairs)
        resolved = sum(p["distance_status"] == "resolved" for p in pairs)
        pair_summaries[pool] = dict(candidate_count=len(ids), pair_count=len(pairs), resolved_pairs=resolved,
                                    unresolved_pairs=len(pairs) - resolved)
        log(f"Pool matrix complete: {pool}, {len(pairs)} pairs")
    merged = [{**row, **{f"global_{k}": v for k, v in match.items() if k != "sequence_id"}}
              for row, match in zip(rows, matches, strict=True)]
    for row in merged:
        row["manual_review"] = (not row["qc_pass"] or bool(row["qc_warnings"]) or row["pool"] != POOLS[0]
                                or row["global_distance_status"] != "resolved" or row["global_closest_reference_tie_count"] > 1)
    warnings = Counter(w for row in rows for w in _tokens(row["qc_warnings"]))
    summary = dict(status="development_provisional", candidate_count=len(rows), reference_count=len(references),
                   qc_pass=sum(r["qc_pass"] for r in rows), qc_fail=sum(not r["qc_pass"] for r in rows),
                   warning_candidates=sum(bool(r["qc_warnings"]) for r in rows), warning_counts=dict(sorted(warnings.items())),
                   terminated_by_eos=sum(r["terminated_by_eos"] for r in rows), hit_generation_cap=sum(r["hit_generation_cap"] for r in rows),
                   eos_audit="exported metadata/hash/length consistency only; checkpoint/raw token replay not performed",
                   baseline_qc_reproduced=True, pool_counts={k: len(v) for k, v in pools.items()},
                   pool_assignment="B frozen labels preserved; not reclassified by C",
                   nearest_reference=dict(pair_count=len(records) * len(references), resolved=sum(m["distance_status"] == "resolved" for m in matches),
                                          unresolved=sum(m["distance_status"] != "resolved" for m in matches),
                                          tied=sum(m["closest_reference_tie_count"] > 1 for m in matches),
                                          exact=sum(m["distance"] == 0 for m in matches)),
                   matrices=pair_summaries, manual_review_count=sum(r["manual_review"] for r in merged),
                   not_assessed=["additional_domains", "transmembrane_topology"],
                   missing_distance_policy="NaN in matrices; blank in TSV; never impute 1",
                   interpretation="QC pass is not family/function validation; hydrophobic runs are not TM predictions; no final selection")
    if sha256_file(config_path) != config_hash or any(sha256_file(paths[k]) != v["sha256"] for k, v in config["inputs"].items()):
        raise RuntimeError("Input changed during run")
    if any(sha256_file(path) != source_hashes[name] for name, path in source_paths.items()):
        raise RuntimeError("Source changed during run")
    with _output_claim(output):
        write_tsv(output / "candidate_qc_extended.tsv", rows, list(rows[0]))
        write_fasta([r for r, qc in zip(records, rows, strict=True) if qc["qc_pass"]], output / "qc_pass_extended.fasta")
        write_tsv(output / "trusted_similarity.tsv", matches, MATCH_FIELDS)
        write_tsv(output / "candidate_audit.tsv", merged, list(merged[0]))
        write_tsv(output / "manual_review.tsv", [r for r in merged if r["manual_review"]], list(merged[0]))
        for short, (ids, matrix, pairs) in matrices.items():
            np.save(output / f"{short}_distance.npy", matrix, allow_pickle=False)
            write_json(output / f"{short}_distance_ids.json", ids)
            write_tsv(output / f"{short}_pairwise.tsv", pairs, PAIR_FIELDS)
        write_json(output / "qc_similarity_summary.json", summary)
        report = ("# C：真实 VAE 候选交接检查\n\n"
                  f"开发性结果，非最终科学结论。候选 {len(rows)} 条，参考 {len(references)} 条。\n\n"
                  f"基础 QC：{summary['qc_pass']} 通过，{summary['qc_fail']} 失败；{summary['warning_candidates']} 条有警告。\n"
                  f"EOS：{summary['terminated_by_eos']}；生成上限：{summary['hit_generation_cap']}。仅核验导出元数据、哈希和长度一致性，未重放模型。\n\n"
                  "## 结果\n\n"
                  f"候选池（沿用 B）：{json.dumps(summary['pool_counts'], ensure_ascii=False)}。\n\n"
                  f"全局最近参考：{summary['nearest_reference']['resolved']} 可解析，{summary['nearest_reference']['unresolved']} 未解析，"
                  f"{summary['nearest_reference']['tied']} 条存在并列。\n\n"
                  f"警告计数（可重叠）：{json.dumps(summary['warning_counts'], ensure_ascii=False)}。\n\n"
                  f"人工复核清单 {summary['manual_review_count']} 条：QC 失败、警告、非主池、未解析或并列之一。\n\n"
                  "## 使用边界\n\n"
                  "全局 BLOSUM62 比对；distance = 1 - M/max(Lq,Lr)，只对可靠比对赋值。"
                  "不是 B 的局部比对结果；最近参考不能替代家族鉴定。\n\n"
                  "main / ambiguous 分别提供矩阵、顺序 ID 和逐对证据；未解析为 NaN（TSV 留空），不可补成 1。"
                  "qc_pass_extended.fasta 含所有基础 QC 通过者，不等于主候选池。\n\n"
                  "**尚未检查：额外结构域、跨膜拓扑。** 连续疏水片段只作警告，不是跨膜区预测；"
                  "家族标签沿用 B，没有重新验证或做最终候选推荐。\n")
        (output / "report.md").write_text(report, encoding="utf-8", newline="\n")
        outputs = {p.name: sha256_file(p) for p in sorted(output.iterdir()) if p.is_file() and not p.name.startswith(".")}
        write_json(output / "qc_similarity_manifest.json", dict(module_version="1.0.0", config=config,
                   config_sha256=config_hash, inputs=config["inputs"], outputs=outputs, source_sha256=source_hashes,
                   effective_quality=quality, effective_alignment=asdict(settings),
                   backend=dict(mode="global", matrix="BLOSUM62", epsilon=1e-6,
                                alignment_tie_policy="lexicographic sequence orientation; first optimal alignment",
                                reference_tie_policy="exact rational effective identity; lexicographic ID representative"),
                   versions=dict(python=platform.python_version(), biopython=Bio.__version__, numpy=np.__version__, pyyaml=yaml.__version__)))
    return summary
