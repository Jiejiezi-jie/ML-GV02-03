"""Check the published 1,000-sequence batch and its QC evidence."""
import json
from pathlib import Path

import pandas as pd

from gv_eval.io import sha256_file

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data/processed/gv02_03_v2/vae_member_a_v1"
RESULTS = ROOT / "results/gv02_03_c_handoff_member_a_v1"


def test_qc_run_counts_and_filtering_are_auditable():
    qc = pd.read_csv(PROCESSED / "candidate_qc.tsv", sep="\t")
    summary = json.loads((PROCESSED / "qc_summary.json").read_text())
    assert len(qc) == 1000
    assert qc.sequence_id.is_unique
    assert qc.qc_pass.sum() == summary["pass_count"] == 1000
    assert summary["failure_count"] == 0
    assert qc.composition_outlier_warning.sum() == 174
    assert summary["generation_cap_count"] == 28


def test_qc_manifest_hashes_match_current_artifacts():
    manifest = json.loads((PROCESSED / "qc_manifest.json").read_text())
    for section in ("inputs", "outputs"):
        for relative, expected in manifest[section].items():
            assert sha256_file(ROOT / relative) == expected


def test_global_alignment_manifest_and_missing_pair_counts_match():
    manifest = json.loads((RESULTS / "qc_similarity_manifest.json").read_text())
    for relative, expected in manifest["outputs"].items():
        assert sha256_file(RESULTS / relative) == expected
    summary = json.loads((RESULTS / "qc_similarity_summary.json").read_text())
    assert summary["baseline_qc_reproduced"] is True
    main = summary["matrices"]["main_supported_gvpa"]
    assert main["resolved_pairs"] + main["unresolved_pairs"] == 147 * 146 // 2
    assert main["unresolved_pairs"] == 3377
