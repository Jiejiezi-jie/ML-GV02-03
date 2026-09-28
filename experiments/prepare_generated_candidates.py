#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.generated_candidates import (  # noqa: E402
    build_generated_candidate_pools,
    prepare_generated_qc,
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Audit generated candidates and build family-gated pools")
    parser.add_argument("stage", choices=("qc", "pools"))
    parser.add_argument("--config", default="configs/vae_member_a_v1_candidates.yaml")
    args = parser.parse_args()
    result = (
        prepare_generated_qc(ROOT, args.config)
        if args.stage == "qc"
        else build_generated_candidate_pools(ROOT, args.config)
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
