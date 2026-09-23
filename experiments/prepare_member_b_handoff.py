"""Freeze a conservative, explicitly provisional A-to-B reference release.

No relabeling and no performance claims: quarantine supported references in
mixed-annotation clusters and clusters containing the best A match of any
unresolved negative control. Preserve the full original classifier artifacts.
"""
from __future__ import annotations

import csv
from collections import defaultdict, Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.io import read_fasta, write_fasta, write_tsv, write_json, sha256_file


def read_table(path):
    with path.open(encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def quarantine_reasons(seeds, conflicts):
    by_id = {r["reference_id"]: r for r in seeds}
    labels = defaultdict(set)
    reasons = defaultdict(set)
    for row in seeds:
        labels[row["cluster_id"]].add(row["annotated_family"])
    for cluster, families in labels.items():
        if len(families) > 1:
            reasons[cluster].add("mixed_annotation_cluster")
    for row in conflicts:
        for field in ("sequence_id", "similarity_best_reference"):
            identifier = row[field]
            if identifier not in by_id:
                raise ValueError(f"Conflict reference missing from seed audit: {identifier}")
            reason = "unresolved_negative_cluster" if field == "sequence_id" else "nearest_gvpa_to_unresolved_negative_cluster"
            reasons[by_id[identifier]["cluster_id"]].add(reason)
    return by_id, reasons


def prepare(root=ROOT):
    ref = root / "data/processed/gv02_03_v2/reference"
    result = root / "results/gv02_03_v2/family"
    output = root / "data/processed/gv02_03_v2/member_b_handoff"
    seeds = read_table(ref / "seed_manifest.tsv")
    conflicts = read_table(result / "unresolved_negative_controls.tsv")
    audit = read_table(ref / "reference_classification.tsv")
    by_id, quarantine = quarantine_reasons(seeds, conflicts)
    references = list(read_fasta(ref / "high_confidence_gvpa.fasta"))
    rows, kept, rejected = [], [], []
    for record in references:
        cluster = by_id[record.identifier]["cluster_id"]
        reasons = sorted(quarantine.get(cluster, set()))
        origins = [r for r in audit if r["reference_id"] == record.identifier and r["reference_status"] == "high_confidence_gvpa"]
        if not origins:
            raise ValueError(f"No accepted provenance: {record.identifier}")
        (rejected if reasons else kept).append(record)
        rows.append(dict(sequence_id=record.identifier, sequence_sha256=origins[0]["sequence_sha256"],
                         cluster_id=cluster, length=len(record.sequence),
                         release_status="quarantined" if reasons else "provisional_training_eligible",
                         reasons=";".join(reasons) if reasons else "annotation_and_out_of_cluster_evidence_agree;no_known_cluster_conflict",
                         source_accessions=";".join(sorted({r["source_id"] for r in origins}))))
    if not kept:
        raise ValueError("No references remain after conservative quarantine")
    write_fasta(kept, output / "gvpa_training_eligible.fasta")
    write_fasta(rejected, output / "quarantined_gvpa.fasta")
    write_tsv(output / "release_audit.tsv", rows, list(rows[0]))
    # Carry the entire conflict table, including strong A-like matches, not just
    # the negative examples that the classifier already handles successfully.
    write_tsv(output / "unresolved_negative_controls.tsv", conflicts, list(conflicts[0]))
    files = [output / name for name in ("gvpa_training_eligible.fasta", "quarantined_gvpa.fasta",
                                       "release_audit.tsv", "unresolved_negative_controls.tsv")]
    inputs = [ref / "high_confidence_gvpa.fasta", ref / "seed_manifest.tsv",
              ref / "reference_classification.tsv", result / "unresolved_negative_controls.tsv",
              result / "acceptance.json"]
    manifest = dict(release_id="member_a_v1_conservative_provisional", status="provisional_not_scientifically_final",
        permitted_use="development_and_explicitly_provisional_training",
        final_training_gate="requires resolution or documented expert adjudication of family conflicts",
        total_supported_gvpa=len(references), training_eligible=len(kept), quarantined=len(rejected),
        unresolved_negative_controls=len(conflicts),
        quarantine_reason_counts=dict(Counter(reason for row in rows if row["release_status"] == "quarantined"
                                               for reason in row["reasons"].split(";"))),
        inputs={p.relative_to(root).as_posix():sha256_file(p) for p in inputs},
        outputs={p.relative_to(root).as_posix():sha256_file(p) for p in files},
        code_sha256=sha256_file(Path(__file__)),
        limitations=["Quarantine is a risk-control policy, not evidence that remaining sequences are experimentally confirmed.",
                     "CD-HIT clusters are heuristic; cross-cluster nearest-neighbor checks remain necessary.",
                     "Do not reuse the old 1224-sequence generator split; rebuild from this pinned release."])
    write_json(output / "release_manifest.json", manifest)
    print({k:manifest[k] for k in ("status", "total_supported_gvpa", "training_eligible", "quarantined")})
    return manifest


if __name__ == "__main__":
    prepare()
