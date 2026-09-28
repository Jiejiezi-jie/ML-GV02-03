#!/usr/bin/env python3
"""Combine existing candidate QC and nearest-reference runs into a review report."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.qc_report import run_qc_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/qc_report.yaml")
    parser.add_argument("--output-dir", required=True, help="New or empty directory; relative to repository root")
    args = parser.parse_args()
    print(json.dumps(run_qc_report(ROOT, args.config, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
