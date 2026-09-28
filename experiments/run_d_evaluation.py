#!/usr/bin/env python3
"""Run downstream evaluation using committed B/C artifacts; no training required."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.d_evaluation import run_evaluation, verify_reproducibility

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/d_evaluation.json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-existing", action="store_true", help="Verify published results and run an independent repeat")
    args = parser.parse_args()
    if args.verify_existing:
        config = json.loads((ROOT / args.config).read_text())
        summary = verify_reproducibility(ROOT, ROOT / args.config, args.output or ROOT / config["output_dir"])
    else:
        summary = run_evaluation(ROOT, ROOT / args.config, args.output)
    print(json.dumps({k: v for k, v in summary.items() if k != "strategy_summary"}, ensure_ascii=False, indent=2))
