#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.generation_data import prepare_generator_data  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare cluster-isolated GvpA data for the V2 sequence generator"
    )
    parser.add_argument("--config", default="configs/generation.yaml")
    args = parser.parse_args()
    manifest = prepare_generator_data(ROOT, args.config)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
