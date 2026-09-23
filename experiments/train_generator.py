#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.generation import train_sequence_vae  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train the V2 GRU sequence VAE with validation and early stopping"
    )
    parser.add_argument("--config", default="configs/generation.yaml")
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda, or cuda:N")
    parser.add_argument("--resume", default=None, help="Checkpoint used for exact resume")
    parser.add_argument("--maximum-epochs", type=int, default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument(
        "--skip-test-evaluation",
        action="store_true",
        help="For smoke checks only: do not load or evaluate the held-out test split",
    )
    parser.add_argument(
        "--allow-provisional-data",
        action="store_true",
        help="Development only: permit the current not-yet-verified GvpA/GvpJ split",
    )
    args = parser.parse_args()
    artifacts = train_sequence_vae(
        ROOT,
        args.config,
        device_name=args.device,
        allow_provisional_data=args.allow_provisional_data,
        resume_checkpoint=args.resume,
        maximum_epochs=args.maximum_epochs,
        output_dir=args.output_dir,
        evaluate_test=not args.skip_test_evaluation,
    )
    print(
        json.dumps(
            {
                "best_checkpoint": str(artifacts.best_checkpoint),
                "last_checkpoint": str(artifacts.last_checkpoint),
                "history": str(artifacts.history_path),
                "curve": str(artifacts.curve_path),
                "manifest": str(artifacts.manifest_path),
                "best_epoch": artifacts.best_epoch,
                "stopped_epoch": artifacts.stopped_epoch,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
