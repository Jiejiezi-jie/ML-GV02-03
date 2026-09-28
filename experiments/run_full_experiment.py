#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

def main() -> None:
    parser = argparse.ArgumentParser(description="Run the reproducible GV02-03 evaluation experiment")
    parser.add_argument("--workflow", choices=["frozen-v2", "legacy"],
                        help="Default: frozen-v2; an explicit YAML config retains the legacy workflow")
    parser.add_argument("--config")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-existing", action="store_true")
    args = parser.parse_args()
    workflow = args.workflow or ("legacy" if args.config and Path(args.config).suffix in (".yaml", ".yml") else "frozen-v2")
    if workflow == "frozen-v2":
        from gv_eval.release import run_release

        summary = run_release(ROOT, ROOT / (args.config or "configs/d_evaluation.json"),
                              output=args.output, verify_existing=args.verify_existing)
    else:
        if args.output or args.verify_existing:
            parser.error("--output and --verify-existing apply only to frozen-v2")
        from gv_eval.pipeline import run_pipeline

        summary = run_pipeline(ROOT, args.config or "configs/gv02_03.yaml")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
