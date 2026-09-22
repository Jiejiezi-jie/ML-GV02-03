#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.family_pipeline import run_family_pipeline


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Member A: audit references and classify Gvp families")
    parser.add_argument("--config", default="configs/family_classification.yaml")
    parser.add_argument("--candidates", help="Override candidate FASTA, relative to project root or absolute")
    args = parser.parse_args()
    run_family_pipeline(ROOT, args.config, args.candidates)
