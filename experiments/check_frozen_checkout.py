#!/usr/bin/env python3
"""Validate a source snapshot in an isolated checkout and fresh virtualenv."""
from __future__ import annotations

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


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="gv-frozen-checkout-") as directory:
        sandbox = Path(directory)
        checkout = sandbox / "checkout"
        subprocess.run(["git", "clone", "--quiet", "--no-hardlinks", str(ROOT), str(checkout)], check=True)
        # Include the exact pending patch under review, without committing or changing ROOT.
        changes = subprocess.check_output(["git", "diff", "--name-only", "-z", "HEAD"], cwd=ROOT).decode().split("\0")
        untracked = subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "-z"], cwd=ROOT).decode().split("\0")
        snapshot = {name for name in changes + untracked
                    if name and name != "results/gv02_03_v2/frozen_checkout_acceptance.json"}
        for name in snapshot:
            source = ROOT / name
            if source.is_file():
                target = checkout / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
        subprocess.run([sys.executable, "-m", "venv", "--system-site-packages", str(sandbox / "venv")], check=True)
        python = sandbox / "venv/bin/python"
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env["PYTHONNOUSERSITE"] = "1"
        env["MPLCONFIGDIR"] = str(sandbox / "matplotlib")
        env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        output = sandbox / "batch"
        subprocess.run([str(python), "experiments/run_full_experiment.py", "--output", str(output)],
                       cwd=checkout, env=env, check=True)
        subprocess.run([str(python), "-m", "pytest", "-q"], cwd=checkout, env=env, check=True)
        current = json.loads((ROOT / "results/gv02_03_v2/d_evaluation_member_a_v1/manifest.json").read_text())
        repeated = json.loads((output / "manifest.json").read_text())
        if current != repeated:
            raise ValueError("Isolated checkout differs from the published evaluation")
        report = {"passed": True, "compared_artifact_count": len(repeated["output_sha256"]),
                  "source_commits": repeated["source_commits"],
                  "python": sys.version.split()[0],
                  "isolation": "Separate Git checkout and new virtualenv; no PYTHONPATH or user-site packages. Dependencies inherited from the selected interpreter via system-site-packages; not a fresh package installation.",
                  "snapshot_sha256": {name: sha256_file(ROOT / name) for name in sorted(snapshot)
                                      if (ROOT / name).is_file()},
                  "scope": "Frozen B/C artifacts through evaluation, selection and all repository tests; no training or sampling replay."}
        write_json(ROOT / "results/gv02_03_v2/frozen_checkout_acceptance.json", report)
    print("Isolated frozen checkout acceptance passed.")


if __name__ == "__main__":
    main()
