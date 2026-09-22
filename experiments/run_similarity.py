#!/usr/bin/env python3
"""Compute auditable global distances for one explicitly selected FASTA pool."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.similarity import run_similarity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/similarity.yaml")
    args = parser.parse_args()
    print(json.dumps(run_similarity(ROOT, args.config), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
