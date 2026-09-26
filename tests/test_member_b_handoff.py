import importlib.util
import csv
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location("handoff", Path(__file__).resolve().parents[1] / "experiments/prepare_member_b_handoff.py")
handoff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(handoff)
from gv_eval.io import sha256_file, read_fasta


def test_quarantine_includes_mixed_clusters_and_negative_nearest_match():
    seeds = [dict(reference_id=i, cluster_id=c, annotated_family=f) for i,c,f in
             [("a","mixed","gvpa"),("j","mixed","gvpj"),("a2","near","gvpa"),
              ("j2","negative","gvpj"),("a3","clean","gvpa")]]
    _, reasons = handoff.quarantine_reasons(seeds, [dict(sequence_id="j2",similarity_best_reference="a2")])
    assert reasons["mixed"] == {"mixed_annotation_cluster"}
    assert reasons["near"] == {"nearest_gvpa_to_unresolved_negative_cluster"}
    assert reasons["negative"] == {"unresolved_negative_cluster"}
    assert "clean" not in reasons


def test_release_and_development_split_are_traceable_and_isolated():
    root = Path(__file__).resolve().parents[1]
    directory = root / "data/processed/gv02_03_v2/member_b_handoff"
    release = json.loads((directory / "release_manifest.json").read_text(encoding="utf-8"))
    assert release["status"] == "provisional_not_scientifically_final"
    for section in ("inputs", "outputs"):
        for path, digest in release[section].items():
            assert sha256_file(root / path) == digest
    assert release["code_sha256"] == sha256_file(root / "experiments/prepare_member_b_handoff.py")
    split = json.loads((directory / "development_split/split_manifest.json").read_text(encoding="utf-8"))
    assert split["release_sha256"] == sha256_file(directory / "release_manifest.json")
    for path, digest in split["outputs"].items():
        assert sha256_file(root / path) == digest
    with (directory / "development_split/sequence_manifest.tsv").open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    for key in ("sequence_id", "sequence_sha256", "cluster_id"):
        seen = {}
        for row in rows:
            assert seen.setdefault(row[key], row["split"]) == row["split"]
    eligible = {r.identifier for r in read_fasta(directory / "gvpa_training_eligible.fasta")}
    rejected = {r.identifier for r in read_fasta(directory / "quarantined_gvpa.fasta")}
    assert eligible.isdisjoint(rejected)
    assert {r["sequence_id"] for r in rows} == eligible
    assert len(eligible) + len(rejected) == release["total_supported_gvpa"]
