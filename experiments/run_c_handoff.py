#!/usr/bin/env python3
"""Audit frozen B delivery and run C QC/global comparisons without modifying B."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.c_handoff import run_handoff


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff-root", required=True, help="Read-only B checkout/archive root")
    parser.add_argument("--config", default="configs/c_handoff_member_a_v1.json")
    parser.add_argument("--output-dir", required=True, help="New directory outside the B handoff root")
    args = parser.parse_args()
    summary = run_handoff(args.handoff_root, ROOT / args.config, ROOT / args.output_dir,
                          progress=lambda message: print(message, flush=True))
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
