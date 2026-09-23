import csv
import json
from copy import deepcopy

import numpy as np
import pytest

from gv_eval import c_handoff as c
from gv_eval.io import FastaRecord, sequence_sha256, sha256_file, write_fasta, write_json, write_tsv


SEQ = "ACDEFGHIKLMNPQRSTVWY"


def record(identifier, sequence=SEQ):
    return FastaRecord(identifier, identifier, sequence)


def metadata(records):
    return [dict(sequence_id=r.identifier, sequence_sha256=sequence_sha256(r.sequence),
                 sequence_length=str(len(r.sequence)), raw_token_length=str(len(r.sequence) + 1),
                 terminated_by_eos="true", hit_generation_cap="false", stop_reason="eos")
            for r in records]


def manifest(count=1):
    return dict(candidate_count=count, terminated_by_eos=count, hit_generation_cap=0,
                generation_parameters=dict(minimum_sequence_length=1, maximum_sequence_length=40))


def test_metadata_join_uses_ids_not_row_order():
    records = [record("a"), record("b", SEQ + "A")]
    audited = c.audit_generation(records, metadata(records)[::-1], manifest(2))
    assert list(audited) == ["a", "b"]
    assert audited["b"]["terminated_by_eos"] is True


@pytest.mark.parametrize("field,value", [
    ("sequence_sha256", "wrong"), ("sequence_length", "2"), ("raw_token_length", "20"),
    ("terminated_by_eos", "yes"), ("hit_generation_cap", "true"), ("stop_reason", "length_cap"),
])
def test_bad_generation_fields_rejected(field, value):
    rows = metadata([record("a")])
    rows[0][field] = value
    with pytest.raises(ValueError):
        c.audit_generation([record("a")], rows, manifest())


@pytest.mark.parametrize("kind", ["missing", "duplicate", "extra", "count", "eos_count", "empty", "nonstandard"])
def test_generation_cardinality_and_sequences_rejected(kind):
    records = [record("a")]
    rows, gen = metadata(records), manifest()
    if kind == "missing": rows = []
    if kind == "duplicate": rows *= 2
    if kind == "extra": rows += metadata([record("b")])
    if kind == "count": gen["candidate_count"] = 2
    if kind == "eos_count": gen["terminated_by_eos"] = 0
    if kind == "empty": records = [record("a", "")]
    if kind == "nonstandard": records = [record("a", "X")]
    with pytest.raises(ValueError):
        c.audit_generation(records, rows, gen)


def test_cap_requires_exact_limit_and_eos_at_limit_is_not_cap():
    records = [record("a", SEQ * 2)]
    rows = metadata(records)
    assert c.audit_generation(records, rows, manifest())["a"]["hit_generation_cap"] is False
    rows[0].update(terminated_by_eos="false", hit_generation_cap="true", stop_reason="length_cap", raw_token_length="40")
    gen = manifest()
    gen.update(terminated_by_eos=0, hit_generation_cap=1)
    assert c.audit_generation(records, rows, gen)["a"]["hit_generation_cap"] is True
    gen["generation_parameters"]["maximum_sequence_length"] = 41
    with pytest.raises(ValueError, match="cap"):
        c.audit_generation(records, rows, gen)


@pytest.fixture
def fixture(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    records = [record("main", SEQ + "AAAAAA"), record("ambiguous", "W" * 26), record("excluded", SEQ)]
    sources = dict(candidates="candidates.fasta", generation_metadata="metadata.tsv", generation_manifest="generation.json",
                   training_fasta="training.fasta", reference_fasta="reference.fasta", b_config="b.yaml",
                   candidate_qc="qc.tsv", family_classification="family.tsv", pool_audit="pools.tsv",
                   main_supported_gvpa="main.fasta", ambiguous_exploration="ambiguous.fasta", excluded="excluded.fasta")
    write_fasta(records, root / sources["candidates"])
    write_fasta([record("train", SEQ + "C")], root / sources["training_fasta"])
    write_fasta([record("ref", SEQ)], root / sources["reference_fasta"])
    rows = metadata(records)
    write_tsv(root / sources["generation_metadata"], rows, list(rows[0]))
    gen = manifest(3)
    gen["output_sha256"] = {sources[k]: sha256_file(root / sources[k]) for k in ("candidates", "generation_metadata")}
    write_json(root / sources["generation_manifest"], gen)
    quality = dict(minimum_length=1, maximum_length=40, minimum_unique_residues=1,
                   maximum_single_residue_fraction=1, generation_length_cap=None,
                   exact_training_match_action="exclude", exact_reference_match_action="exclude")
    write_json(root / sources["b_config"], {"quality": quality})
    qc = [dict(sequence_id=r.identifier, sequence_sha256=sequence_sha256(r.sequence),
               qc_pass=str(r.identifier != "excluded"), qc_reasons="exact_reference_match" if r.identifier == "excluded" else "") for r in records]
    write_tsv(root / sources["candidate_qc"], qc, list(qc[0]))
    family = [dict(sequence_id=r.identifier, sequence_sha256=sequence_sha256(r.sequence), sequence_length=len(r.sequence),
                   family_status=s, pf00741_status="pass") for r, s in zip(records, ["supported_gvpa", "gvpa_gvpj_ambiguous", "supported_gvpa"])]
    write_tsv(root / sources["family_classification"], family, list(family[0]))
    pools = [dict(sequence_id=r.identifier, pool=p) for r, p in zip(records, ["main_supported_gvpa", "ambiguous_exploration", "excluded"])]
    write_tsv(root / sources["pool_audit"], pools, list(pools[0]))
    for r, row in zip(records, pools): write_fasta([r], root / sources[row["pool"]])
    config = dict(inputs={k: dict(path=v, sha256=sha256_file(root / v)) for k, v in sources.items()},
                  expected_counts=dict(candidates=3, references=1, main_supported_gvpa=1, ambiguous_exploration=1, excluded=1),
                  patterns=dict(hydrophobic_run_minimum=18, homopolymer_minimum=6, tandem_repeat_minimum_length=12), alignment={})
    config_path = tmp_path / "config.json"
    write_json(config_path, config)
    return root, config_path, config, tmp_path / "output"


def load_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_runner_keeps_pools_warnings_unresolved_and_hashes(fixture):
    root, config_path, config, output = fixture
    before = {k: sha256_file(root / v["path"]) for k, v in config["inputs"].items()}
    summary = c.run_handoff(root, config_path, output)
    assert summary["candidate_count"] == 3
    assert summary["qc_pass"] == 2
    assert summary["nearest_reference"]["unresolved"] == 2
    assert summary["not_assessed"] == ["additional_domains", "transmembrane_topology"]
    rows = {r["sequence_id"]: r for r in load_rows(output / "candidate_qc_extended.tsv")}
    assert rows["main"]["qc_pass"] == "True"
    assert "homopolymer_warning" in rows["main"]["qc_warnings"]
    assert rows["main"]["qc_failures"] == ""
    assert rows["excluded"]["qc_failures"] == "exact_reference_match"
    assert rows["ambiguous"]["pool"] == "ambiguous_exploration"
    assert json.loads((output / "main_distance_ids.json").read_text()) == ["main"]
    np.testing.assert_array_equal(np.load(output / "main_distance.npy"), [[0.]])
    matches = {r["sequence_id"]: r for r in load_rows(output / "trusted_similarity.tsv")}
    assert matches["ambiguous"]["distance"] == ""
    assert matches["excluded"]["distance"] == "0.0"
    audit = json.loads((output / "qc_similarity_manifest.json").read_text())
    for name, digest in audit["outputs"].items(): assert sha256_file(output / name) == digest
    assert before == {k: sha256_file(root / v["path"]) for k, v in config["inputs"].items()}
    replay = output.with_name("replay")
    c.run_handoff(root, config_path, replay)
    for name in audit["outputs"]: assert (output / name).read_bytes() == (replay / name).read_bytes()
    with pytest.raises(ValueError, match="empty"):
        c.run_handoff(root, config_path, output)


@pytest.mark.parametrize("kind", ["hash", "baseline", "pool", "pool_sequence", "family_hash", "count", "escape", "override"])
def test_runner_rejects_inconsistent_source_without_outputs(fixture, kind):
    root, config_path, config, output = fixture
    key = None
    if kind == "hash": (root / "candidates.fasta").write_text(">wrong\nAAAA\n")
    if kind == "baseline":
        key = "candidate_qc"
        rows = load_rows(root / "qc.tsv"); rows[0]["qc_pass"] = "False"
        write_tsv(root / "qc.tsv", rows, list(rows[0]))
    if kind == "pool":
        key = "pool_audit"
        rows = load_rows(root / "pools.tsv"); rows[0]["pool"] = "excluded"
        write_tsv(root / "pools.tsv", rows, list(rows[0]))
    if kind == "pool_sequence":
        key = "main_supported_gvpa"; write_fasta([record("main")], root / "main.fasta")
    if kind == "family_hash":
        key = "family_classification"
        rows = load_rows(root / "family.tsv"); rows[0]["sequence_sha256"] = "bad"
        write_tsv(root / "family.tsv", rows, list(rows[0]))
    if kind == "count": config["expected_counts"]["candidates"] = 4
    if kind == "escape": config["inputs"]["candidates"]["path"] = "../config.json"
    if kind == "override": config["patterns"]["minimum_length"] = 100
    if key: config["inputs"][key]["sha256"] = sha256_file(root / config["inputs"][key]["path"])
    write_json(config_path, config)
    with pytest.raises(ValueError): c.run_handoff(root, config_path, output)
    assert not output.exists()
