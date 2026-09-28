from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from gv_eval.c_domains import DomainHit, run, scan_profiles, warning_reasons
from gv_eval.io import FastaRecord, sha256_file


SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 3


def _hit(accession: str, start: int, end: int) -> DomainHit:
    return DomainHit("candidate", "c1", accession, accession, start, end, 50.0, 1e-10,
                     1, 39, 39)


def _model(path: Path, gathering_cutoff: float = 1.0) -> None:
    pyhmmer = pytest.importorskip("pyhmmer")
    amino = pyhmmer.easel.Alphabet.amino()
    msa = pyhmmer.easel.TextMSA(
        name=b"mini",
        sequences=[pyhmmer.easel.TextSequence(name=f"s{i}".encode(), sequence=SEQUENCE)
                   for i in range(3)],
    )
    hmm = pyhmmer.plan7.Builder(amino).build_msa(
        msa.digitize(amino), pyhmmer.plan7.Background(amino)
    )[0]
    hmm.accession = "PF00741.1"
    hmm.cutoffs.gathering = (gathering_cutoff, gathering_cutoff)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wb") as stream:
        hmm.write(stream)


def test_nonoverlapping_non_target_domain_is_review_only():
    reasons = warning_reasons([_hit("PF00741.24", 10, 48), _hit("PF99999.1", 60, 80)], 1)
    assert reasons == ("non_target_extra_domain",)


def test_overlapping_alternative_family_is_not_called_extra_domain():
    reasons = warning_reasons([_hit("PF00741.24", 10, 48), _hit("PF12345.2", 20, 40)], 1)
    assert reasons == ("overlapping_alternative_family",)


def test_partial_overlap_is_not_called_extra_domain():
    reasons = warning_reasons([_hit("PF00741.24", 10, 48), _hit("PF12345.2", 40, 60)], 1)
    assert reasons == ("overlapping_alternative_family",)


def test_multiple_target_segments_are_reviewed_without_claiming_extra_copies():
    reasons = warning_reasons([_hit("PF00741.24", 5, 40), _hit("PF00741.24", 70, 105)], 1)
    assert reasons == ("multiple_pf00741_segments",)
    assert warning_reasons([_hit("PF00741.24", 5, 40)], 1) == ()


def test_real_pyhmmer_scan_uses_curated_threshold_and_coordinates(tmp_path: Path):
    model = tmp_path / "mini.hmm"
    _model(model)
    records = {
        "candidate": [FastaRecord("c1", "c1", SEQUENCE)],
        "reference": [FastaRecord("r1", "r1", SEQUENCE)],
    }
    hits, model_count = scan_profiles(model, records, cpus=1)
    assert model_count == 1
    assert [(h.role, h.sequence_id, h.accession, h.start, h.end) for h in hits] == [
        ("candidate", "c1", "PF00741.1", 1, 60),
        ("reference", "r1", "PF00741.1", 1, 60),
    ]
    assert [(h.hmm_start, h.hmm_end, h.hmm_length) for h in hits] == [(1, 60, 60), (1, 60, 60)]


def test_run_preserves_pool_and_writes_hashed_evidence(tmp_path: Path):
    model = tmp_path / "mini.hmm"
    _model(model)
    handoff = tmp_path / "handoff"
    handoff.mkdir()
    (handoff / "candidates.fasta").write_text(f">c1\n{SEQUENCE}\n", encoding="ascii")
    (handoff / "references.fasta").write_text(f">r1\n{SEQUENCE}\n", encoding="ascii")
    (handoff / "pools.tsv").write_text("sequence_id\tpool\nc1\tmain_supported_gvpa\n", encoding="ascii")
    config = {
        "inputs": {
            key: {"path": name, "sha256": sha256_file(handoff / name)}
            for key, name in (("candidates", "candidates.fasta"),
                              ("reference_fasta", "references.fasta"),
                              ("pool_audit", "pools.tsv"))
        },
        "expected_counts": {"candidates": 1, "references": 1},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "results"
    expected_md5 = hashlib.md5(model.read_bytes()).hexdigest()
    summary = run(config_path, handoff, model, output, cpus=1, expected_pfam_md5=expected_md5)

    assert summary["reference_max_pf00741_hit_segments"] == 1
    assert summary["candidate_review_status_counts"] == {"no_extra_domain_evidence": 1}
    with (output / "candidate_domains.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert [(r["sequence_id"], r["pool"], r["review_status"]) for r in rows] == [
        ("c1", "main_supported_gvpa", "no_extra_domain_evidence")
    ]
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    for name, digest in manifest["output_sha256"].items():
        assert sha256_file(output / name) == digest
    with pytest.raises(FileExistsError, match="not empty"):
        run(config_path, handoff, model, output, cpus=1, expected_pfam_md5=expected_md5)


def test_database_hash_mismatch_prevents_scanning(tmp_path: Path):
    model = tmp_path / "mini.hmm"
    _model(model)
    with pytest.raises(ValueError, match="Pfam MD5 mismatch"):
        scan_profiles(model, {"candidate": [FastaRecord("c1", "c1", SEQUENCE)]},
                      cpus=1, expected_pfam_md5="0" * 32)


def test_gzip_model_respects_high_gathering_cutoff(tmp_path: Path):
    model = tmp_path / "mini.hmm.gz"
    _model(model, gathering_cutoff=1000.0)
    hits, count = scan_profiles(model, {"candidate": [FastaRecord("c1", "c1", SEQUENCE)]}, cpus=1)
    assert count == 1
    assert hits == []


@pytest.mark.parametrize("result_dir", [
    "gv02_03_c_domains_member_a_v1", "gv02_03_c_tm_member_a_v1",
])
def test_generated_result_bytes_are_preserved_by_git(result_dir: str):
    repo = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        ["git", "check-attr", "text", "--", f"results/{result_dir}/summary.json"],
        cwd=repo, check=True, capture_output=True, text=True,
    )
    assert result.stdout.strip().endswith("text: unset")
