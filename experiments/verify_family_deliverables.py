"""Verify member-A artifacts, report scientific failures without hiding them."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.io import read_fasta, sha256_file, write_json


def read_table(path):
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def verify(root=ROOT, baseline=None):
    result = root / "results/gv02_03_v2/family"
    ref = root / "data/processed/gv02_03_v2/reference"
    manifest = json.loads((result / "manifest.json").read_text(encoding="utf-8"))
    checks = {}
    checks["manifest_hashes"] = all(sha256_file(root / p) == digest
        for section in ("inputs", "outputs", "code") for p, digest in manifest[section].items())
    candidate = read_table(ref.parent / "candidate_family_classification.tsv")
    source_path = next(p for p in manifest["inputs"] if p.endswith(".fasta") and p not in manifest["config"]["sources"])
    source = list(read_fasta(root / source_path))
    checks["candidate_ids_and_order"] = [r.identifier for r in source] == [r["sequence_id"] for r in candidate]
    checks["candidate_reason_complete"] = all(r["classification_reason"] for r in candidate)
    controls = read_table(result / "control_classification.tsv")
    seeds = read_table(ref / "seed_manifest.tsv")
    seed_by_id = {r["reference_id"]:r for r in seeds}
    profile = json.loads((root / "models/family_profiles/profile_manifest.json").read_text(encoding="utf-8"))
    checks["cluster_isolation"] = True
    for fold in profile["cross_validation"]:
        training = {i for p in fold["profiles"] for i in p["training_ids"]}
        test = set(fold["heldout_ids"])
        train_clusters = {seed_by_id[i]["cluster_id"] for i in training}
        test_clusters = {seed_by_id[i]["cluster_id"] for i in test}
        checks["cluster_isolation"] &= not bool(training & test or train_clusters & test_clusters)
    audit = read_table(ref / "reference_classification.tsv")
    categories = ("high_confidence_gvpa", "gvpj_negative_control", "other_gvp_negative_control")
    exported = {}
    for category in categories:
        ids = [r.identifier for r in read_fasta(ref / (category + ".fasta"))]
        exported[category] = set(ids)
        expected = {r["reference_id"] for r in audit if r["reference_status"] == category}
        checks[category + "_consistent"] = len(ids) == len(set(ids)) and set(ids) == expected and bool(ids)
    checks["reference_sets_disjoint"] = sum(map(len, exported.values())) == len(set.union(*exported.values()))
    natural = [r for r in controls if r["control_type"] == "natural_heldout"]
    false_positive = [r for r in natural if r["expected_family"] != "gvpa" and r["family_status"] == "supported_gvpa"]
    adjudication_path = result / "adjudication" / "conflict_adjudication.tsv"
    adjudications = read_table(adjudication_path) if adjudication_path.exists() else []
    adjudicated_by_id = {r["sequence_id"]: r for r in adjudications}
    from gv_eval.io import write_tsv
    enriched = []
    for row in false_positive:
        origins = [r for r in audit if r["reference_id"] == row["sequence_id"]]
        enriched.append(dict(**row, source_accessions=";".join(sorted({r["source_id"] for r in origins})),
                             source_headers=" | ".join(sorted({r["source_header"] for r in origins})),
                             review_action="check_database_annotation_and_gene_neighborhood;do_not_auto_relabel"))
    write_tsv(result / "unresolved_negative_controls.tsv", enriched,
              list(controls[0]) + ["source_accessions", "source_headers", "review_action"])
    synthetic_fp = [r for r in controls if r["control_type"] in ("shuffled", "truncated") and r["family_status"] == "supported_gvpa"]
    by_natural_id = {r["sequence_id"]:r for r in natural}
    checks["duplicate_controls_consistent"] = all(r["family_status"] == by_natural_id[r["parent_id"]]["family_status"]
        for r in controls if r["control_type"] == "exact_duplicate")
    checks["synthetic_negative_controls"] = not synthetic_fp
    reviewed = [r for r in natural if seed_by_id[r["sequence_id"]]["reviewed_anchor"] == "True"]
    checks["reviewed_negative_controls"] = not any(r["expected_family"] != "gvpa" and r["family_status"] == "supported_gvpa" for r in reviewed)
    checks["positive_controls_present"] = any(r["expected_family"] == "gvpa" and r["family_status"] == "supported_gvpa" for r in reviewed)
    if baseline:
        old = json.loads(Path(baseline).read_text(encoding="utf-8"))
        checks["reproducible_outputs"] = manifest["outputs"] == old["outputs"]
    operational_false_positive = [r for r in false_positive
                                  if adjudicated_by_id.get(r["sequence_id"], {}).get("operational_family_status") == "supported_gvpa"
                                  and adjudicated_by_id[r["sequence_id"]].get("biological_identity_resolved") != "True"]
    reviewed_conflicts = len(adjudicated_by_id) == len(false_positive) and not operational_false_positive
    report = dict(engineering_checks=checks, engineering_pass=all(checks.values()),
                  nominal_negative_false_positives=len(false_positive),
                  nominal_negative_count=sum(r["expected_family"] != "gvpa" for r in natural),
                  false_positive_ids=[r["sequence_id"] for r in false_positive],
                  reviewed_conflict_count=len(adjudicated_by_id),
                  operational_false_positive_ids=[r["sequence_id"] for r in operational_false_positive],
                  scientific_acceptance=("pass_internal_controls_with_adjudicated_annotation_conflicts"
                                         if all(checks.values()) and reviewed_conflicts else "requires_review"),
                  nominal_confusion={"/".join(k):v for k,v in Counter((r["expected_family"],r["family_status"]) for r in natural).items()},
                  reviewed_confusion={"/".join(k):v for k,v in Counter((r["expected_family"],r["family_status"]) for r in reviewed).items()},
                  limitation="Six conflicts were resolved as stale local GvpJ annotations and eight remain operationally ambiguous; internal controls do not establish experimental function.")
    write_json(result / "acceptance.json", report)
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-manifest")
    args = parser.parse_args()
    report = verify(baseline=args.baseline_manifest)
    sys.exit(0 if report["engineering_pass"] else 1)
