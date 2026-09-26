"""Snapshot/compare stable scientific outputs across complete pipeline runs."""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.io import sha256_file, write_json


def stable_outputs():
    result = ROOT / "results/gv02_03_v2/family"
    ref = ROOT / "data/processed/gv02_03_v2/reference"
    paths = [ref.parent / "candidate_family_classification.tsv", result / "summary.json",
             result / "control_classification.tsv", result / "reference_out_of_fold_scores.tsv"]
    paths += list(ref.glob("*.fasta")) + list(ref.glob("*.tsv"))
    paths += list((ROOT / "models/family_profiles").glob("*.hmm"))
    return {p.relative_to(ROOT).as_posix():sha256_file(p) for p in sorted(paths)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["snapshot", "compare"])
    args = parser.parse_args()
    directory = ROOT / "results/gv02_03_v2/family"
    baseline = directory / "reproducibility_baseline.json"
    if args.mode == "snapshot":
        write_json(baseline, stable_outputs())
    else:
        old = json.loads(baseline.read_text(encoding="utf-8"))
        current = stable_outputs()
        changed = [p for p in sorted(set(old) | set(current)) if old.get(p) != current.get(p)]
        report = dict(equal=not changed, compared_outputs=len(current), changed=changed,
                      run_1=old, run_2=current,
                      scope="Scientific FASTA/TSV/HMM/summary outputs; tool timing logs and provenance manifests excluded. Input text line endings were normalized to original Git LF bytes between runs.")
        write_json(directory / "reproducibility_check.json", report)
        print(json.dumps(dict(equal=not changed, compared_outputs=len(current), changed=changed)))
        sys.exit(bool(changed))
