#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gv_eval.generation import generate_candidates  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate auditable GvpA candidates from a trained sequence VAE"
    )
    parser.add_argument("--config", default="configs/generation.yaml")
    parser.add_argument(
        "--checkpoint", default="models/generator/sequence_vae/best.pt"
    )
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda, or cuda:N")
    parser.add_argument("--candidate-count", type=int, default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument(
        "--allow-provisional-data",
        action="store_true",
        help="Development only: permit a checkpoint trained on provisional data",
    )
    args = parser.parse_args()
    artifacts = generate_candidates(
        ROOT,
        args.config,
        args.checkpoint,
        device_name=args.device,
        allow_provisional_data=args.allow_provisional_data,
        candidate_count=args.candidate_count,
        output_dir=args.output_dir,
    )
    print(
        json.dumps(
            {
                "candidates_fasta": str(artifacts.candidates_fasta),
                "metadata": str(artifacts.metadata_path),
                "manifest": str(artifacts.manifest_path),
                "candidate_count": artifacts.candidate_count,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
