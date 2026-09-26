"""Reproduce the 14-record review from a pinned official NCBI snapshot."""
from __future__ import annotations
import csv
import json
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.adjudication import parse_genbank_proteins, adjudicate_conflict, apply_adjudication
from gv_eval.io import read_fasta, write_fasta, write_tsv, write_json, sha256_file


def read_table(path):
    with path.open(encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream,delimiter="\t"))


def main():
    result = ROOT / "results/gv02_03_v2/family"
    ref = ROOT / "data/processed/gv02_03_v2/reference"
    output = result / "adjudication"
    snapshot = ROOT / "references/conflict_review/ncbi_14_proteins.gb"
    conflicts_path = result / "unresolved_negative_controls.tsv"
    conflicts = read_table(conflicts_path)
    official = parse_genbank_proteins(snapshot.read_text(encoding="utf-8"))
    decisions = []
    requested = set()
    for row in conflicts:
        accession = row["source_accessions"].split(";")[0]
        requested.add(accession)
        decisions.append(adjudicate_conflict(row, official[accession]))
    if requested != set(official):
        raise ValueError("Official snapshot does not match the conflict set")
    write_tsv(output / "conflict_adjudication.tsv", decisions, list(decisions[0]))
    indexed = {r["sequence_sha256"]:r for r in decisions}
    for source, name in [(result / "control_classification.tsv", "controls_reviewed_view.tsv"),
                         (ref.parent / "candidate_family_classification.tsv", "candidates_reviewed_view.tsv")]:
        rows = [apply_adjudication(row,indexed) for row in read_table(source)]
        write_tsv(output / name, rows, list(rows[0]))
    natural = {r.identifier:r for r in read_fasta(result / "natural_unique.fasta")}
    supported = list(read_fasta(ref / "high_confidence_gvpa.fasta"))
    ids = {r.identifier for r in supported}
    additions, ambiguous = [], []
    for decision in decisions:
        record = natural[decision["sequence_id"]]
        if decision["biological_identity_resolved"]:
            if record.identifier not in ids:
                supported.append(record)
                additions.append(record)
        else:
            ambiguous.append(record)
    write_fasta(additions, output / "ncbi_resolved_gvpa.fasta")
    write_fasta(ambiguous, output / "annotation_conflict_ambiguous.fasta")
    write_fasta(supported, output / "gvpa_reference_with_ncbi_corrections.fasta")
    summary = dict(reviewed_records=len(decisions), exact_sequence_matches=len(decisions),
        dispositions=dict(Counter(r["disposition"] for r in decisions)),
        resolved_annotation_changes=len(additions), biological_identity_still_ambiguous=len(ambiguous),
        supported_reference_with_corrections=len(supported),
        scope="Annotation adjudication and operational abstention; not a retrained or independently validated classifier.",
        classifier_errors_not_erased=True, final_scientific_status="partial_resolution_eight_ambiguous",
        member_b="Keep the conservative provisional release; do not restore quarantined clusters solely because six labels were corrected.")
    write_json(output / "summary.json",summary)
    inputs=[snapshot,conflicts_path,ref / "high_confidence_gvpa.fasta",result / "natural_unique.fasta",
            result / "control_classification.tsv",ref.parent / "candidate_family_classification.tsv"]
    outputs=[p for p in output.iterdir() if p.is_file() and p.name!="manifest.json"]
    write_json(output / "manifest.json",dict(inputs={p.relative_to(ROOT).as_posix():sha256_file(p) for p in inputs},
        outputs={p.relative_to(ROOT).as_posix():sha256_file(p) for p in outputs},
        code={p.relative_to(ROOT).as_posix():sha256_file(p) for p in [Path(__file__).resolve(),ROOT / "src/gv_eval/adjudication.py"]},
        retrieved_on="2026-09-22", source_url=(ROOT / "references/conflict_review/request_url.txt").read_text()))
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
