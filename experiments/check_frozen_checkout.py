#!/usr/bin/env python3
"""Check a separate Git source snapshot using the selected Python environment."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.io import sha256_file, write_json
from gv_eval.reproducibility import compare_artifacts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=ROOT / "results/reproductions/checkout_acceptance.json")
    args = parser.parse_args()
    if args.report.exists():
        parser.error("Report already exists; choose a new --report path")
    with tempfile.TemporaryDirectory(prefix="gv-frozen-checkout-", ignore_cleanup_errors=True) as directory:
        sandbox = Path(directory)
        checkout = sandbox / "checkout"
        subprocess.run(["git", "clone", "--quiet", "--no-hardlinks", str(ROOT), str(checkout)], check=True)
        changes = subprocess.check_output(["git", "diff", "--name-only", "-z", "HEAD"], cwd=ROOT).decode().split("\0")
        untracked = subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "-z"], cwd=ROOT).decode().split("\0")
        snapshot = {name for name in changes + untracked if name}
        for name in snapshot:
            source, target = ROOT / name, checkout / name
            if source.is_file():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
            elif target.is_file():
                target.unlink()  # Mirror pending deletions only in this disposable clone.
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env["PYTHONNOUSERSITE"] = "1"
        env["PYTHONUTF8"] = "1"
        env["MPLCONFIGDIR"] = str(sandbox / "matplotlib")
        env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        output = sandbox / "batch"
        subprocess.run([sys.executable, "experiments/run_full_experiment.py", "--output", str(output)],
                       cwd=checkout, env=env, check=True)
        subprocess.run([sys.executable, "-m", "pytest", "tests", "-q"], cwd=checkout, env=env, check=True)
        published = ROOT / "results/gv02_03_v2/d_evaluation_member_a_v1"
        current = json.loads((published / "manifest.json").read_text(encoding="utf-8"))
        repeated = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        for key in ("source_commits", "config", "config_sha256", "input_sha256"):
            if current[key] != repeated[key]:
                raise ValueError(f"Isolated checkout input provenance differs: {key}")
        comparison = compare_artifacts(published, output, current["output_sha256"], repeated["output_sha256"])
        if not comparison["passed"]:
            raise ValueError(f"Isolated checkout scientific outputs differ: {comparison['failed']}")
        report = {"passed": True, "comparison": comparison, "source_commits": repeated["source_commits"],
                  "python": sys.version.split()[0],
                  "isolation": "Separate Git checkout; selected interpreter and its installed dependencies; PYTHONPATH and user-site packages disabled. This check does not install dependencies.",
                  "snapshot_sha256": {name: sha256_file(ROOT / name) for name in sorted(snapshot) if (ROOT / name).is_file()},
                  "deleted_paths": sorted(name for name in snapshot if not (ROOT / name).exists()),
                  "scope": "Frozen B/C artifacts through evaluation, selection and the full repository test suite."}
        write_json(args.report, report)
    print(f"Isolated checkout acceptance passed: {args.report}")


if __name__ == "__main__":
    main()
