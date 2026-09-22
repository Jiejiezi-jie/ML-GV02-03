from __future__ import annotations

import csv
import importlib.metadata
import math
import os
import platform
import random
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
import yaml
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

from .io import (
    FastaRecord,
    read_fasta,
    sequence_sha256,
    sha256_file,
    write_fasta,
    write_json,
    write_tsv,
)
from .vae import (
    ProteinBatchCollator,
    ProteinSequenceDataset,
    ProteinVocabulary,
    SequenceVAE,
    compute_vae_loss,
    kl_beta_for_epoch,
)


@dataclass(frozen=True)
class TrainingArtifacts:
    best_checkpoint: Path
    last_checkpoint: Path
    history_path: Path
    curve_path: Path
    manifest_path: Path
    best_epoch: int
    stopped_epoch: int


@dataclass(frozen=True)
class GenerationArtifacts:
    candidates_fasta: Path
    metadata_path: Path
    manifest_path: Path
    candidate_count: int


def load_generation_config(root: str | Path, config_path: str | Path) -> tuple[Path, dict[str, Any]]:
    root = Path(root).resolve()
    path = Path(config_path)
    if not path.is_absolute():
        path = root / path
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError(f"Generation config is not a mapping: {path}")
    required_sections = {"seed", "model", "data", "vocabulary", "generation", "outputs"}
    missing = sorted(required_sections - set(config))
    if missing:
        raise ValueError(f"Generation config is missing sections: {missing}")
    if config["model"].get("architecture") != "gru":
        raise ValueError("Only the frozen GRU sequence-VAE architecture is supported")
    return path, config


def resolve_device(requested: str) -> torch.device:
    if requested == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(requested)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(f"CUDA device requested but unavailable: {requested}")
    return device


def seed_everything(seed: int) -> None:
    os.environ.setdefault("PYTHONHASHSEED", str(seed))
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)
    if torch.backends.cudnn.is_available():
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True


def build_sequence_vae(config: dict[str, Any], vocabulary: ProteinVocabulary) -> SequenceVAE:
    model_config = config["model"]
    return SequenceVAE(
        vocab_size=len(vocabulary.tokens),
        pad_id=vocabulary.pad_id,
        bos_id=vocabulary.bos_id,
        eos_id=vocabulary.eos_id,
        unk_id=vocabulary.unk_id,
        embedding_dim=int(model_config["embedding_dim"]),
        hidden_dim=int(model_config["hidden_dim"]),
        latent_dim=int(model_config["latent_dim"]),
        encoder_layers=int(model_config.get("encoder_layers", 1)),
        decoder_layers=int(model_config.get("decoder_layers", 1)),
        dropout=float(model_config.get("dropout", 0.0)),
    )


def _loader(
    dataset: ProteinSequenceDataset,
    *,
    batch_size: int,
    shuffle: bool,
    pad_id: int,
    seed: int,
    num_workers: int,
) -> DataLoader:
    generator = torch.Generator()
    generator.manual_seed(seed)
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        collate_fn=ProteinBatchCollator(pad_id),
        generator=generator,
        persistent_workers=bool(num_workers),
    )


def _run_epoch(
    model: SequenceVAE,
    loader: DataLoader,
    *,
    device: torch.device,
    beta: float,
    optimizer: Adam | None,
    gradient_clip_norm: float,
) -> dict[str, float]:
    training = optimizer is not None
    model.train(training)
    totals = {"total": 0.0, "reconstruction": 0.0, "kl": 0.0}
    examples = 0
    tokens = 0
    correct_tokens = 0
    eos_tokens = 0
    correct_eos = 0
    for batch in loader:
        batch_tokens = batch["tokens"]
        batch_lengths = batch["lengths"]
        if not isinstance(batch_tokens, torch.Tensor) or not isinstance(batch_lengths, torch.Tensor):
            raise TypeError("Unexpected batch tensor types")
        batch_tokens = batch_tokens.to(device)
        batch_lengths = batch_lengths.to(device)
        if training:
            optimizer.zero_grad(set_to_none=True)
        with torch.set_grad_enabled(training):
            output = model(batch_tokens, batch_lengths, sample_latent=training)
            loss = compute_vae_loss(output, pad_id=model.pad_id, beta=beta)
            if training:
                loss.total.backward()
                nn.utils.clip_grad_norm_(model.parameters(), gradient_clip_norm)
                optimizer.step()
        batch_size = batch_tokens.size(0)
        totals["total"] += float(loss.total.detach().cpu()) * batch_size
        totals["reconstruction"] += float(loss.reconstruction.detach().cpu()) * loss.token_count
        totals["kl"] += float(loss.kl.detach().cpu()) * batch_size
        examples += batch_size
        tokens += loss.token_count
        predictions = output.logits.detach().argmax(dim=-1)
        target_mask = output.targets.ne(model.pad_id)
        eos_mask = output.targets.eq(model.eos_id)
        correct_tokens += int(
            predictions.eq(output.targets).logical_and(target_mask).sum().cpu()
        )
        eos_tokens += int(eos_mask.sum().cpu())
        correct_eos += int(predictions.eq(model.eos_id).logical_and(eos_mask).sum().cpu())
    if examples == 0 or tokens == 0:
        raise RuntimeError("Data loader produced no usable batches")
    reconstruction = totals["reconstruction"] / tokens
    return {
        "loss": totals["total"] / examples,
        "reconstruction": reconstruction,
        "kl": totals["kl"] / examples,
        "perplexity": math.exp(min(reconstruction, 50.0)),
        "token_accuracy": correct_tokens / tokens,
        "eos_accuracy": correct_eos / eos_tokens if eos_tokens else 0.0,
        "examples": float(examples),
        "tokens": float(tokens),
    }


def _capture_rng_state() -> dict[str, Any]:
    state: dict[str, Any] = {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch_cpu": torch.get_rng_state(),
    }
    if torch.cuda.is_available():
        state["torch_cuda"] = torch.cuda.get_rng_state_all()
    return state


def _restore_rng_state(state: dict[str, Any]) -> None:
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.set_rng_state(state["torch_cpu"])
    if torch.cuda.is_available() and "torch_cuda" in state:
        torch.cuda.set_rng_state_all(state["torch_cuda"])


def _save_checkpoint(
    path: Path,
    *,
    model: SequenceVAE,
    optimizer: Adam,
    epoch: int,
    best_epoch: int,
    best_validation_loss: float,
    epochs_without_improvement: int,
    history: list[dict[str, Any]],
    config: dict[str, Any],
    config_sha256: str,
    vocabulary_path: Path,
    data_hashes: dict[str, str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "format_version": 1,
            "epoch": epoch,
            "best_epoch": best_epoch,
            "best_validation_loss": best_validation_loss,
            "epochs_without_improvement": epochs_without_improvement,
            "model_config": model.model_config(),
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "history": history,
            "config": config,
            "config_sha256": config_sha256,
            "vocabulary_sha256": sha256_file(vocabulary_path),
            "data_sha256": data_hashes,
            "rng_state": _capture_rng_state(),
        },
        path,
    )


def _load_torch_file(path: Path, device: torch.device) -> dict[str, Any]:
    try:
        return torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=device)


def _write_history(path: Path, history: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epoch",
        "kl_beta",
        "train_loss",
        "train_reconstruction",
        "train_kl",
        "train_perplexity",
        "train_token_accuracy",
        "train_eos_accuracy",
        "validation_loss",
        "validation_reconstruction",
        "validation_kl",
        "validation_perplexity",
        "validation_token_accuracy",
        "validation_eos_accuracy",
        "epoch_seconds",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(history)


def _write_training_curve(path: Path, history: list[dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    epochs = [row["epoch"] for row in history]
    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    axis.plot(epochs, [row["train_loss"] for row in history], label="train")
    axis.plot(epochs, [row["validation_loss"] for row in history], label="validation")
    axis.set_xlabel("Epoch")
    axis.set_ylabel("ELBO loss")
    axis.set_title("Sequence VAE training curve")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def _dependency_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {
        "python": platform.python_version(),
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
    }
    for distribution in ("numpy", "PyYAML", "biopython"):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def train_sequence_vae(
    root: str | Path,
    config_path: str | Path,
    *,
    device_name: str = "auto",
    allow_provisional_data: bool = False,
    resume_checkpoint: str | Path | None = None,
    maximum_epochs: int | None = None,
    output_dir: str | Path | None = None,
) -> TrainingArtifacts:
    root = Path(root).resolve()
    resolved_config_path, config = load_generation_config(root, config_path)
    config_sha256 = sha256_file(resolved_config_path)
    verification_status = str(config["data"]["verification_status"])
    if "pending" in verification_status and not allow_provisional_data:
        raise RuntimeError(
            "Training data is provisional; pass allow_provisional_data only for development runs"
        )
    seed = int(config["seed"])
    seed_everything(seed)
    device = resolve_device(device_name)

    output_config = config["outputs"]
    vocabulary_path = root / output_config["vocabulary"]
    vocabulary = ProteinVocabulary.from_json(vocabulary_path)
    data_dir = root / output_config["data_dir"]
    maximum_sequence_length = int(config["model"]["maximum_sequence_length"])
    datasets = {
        split: ProteinSequenceDataset(
            data_dir / f"{split}.fasta",
            vocabulary,
            maximum_sequence_length=maximum_sequence_length,
        )
        for split in ("train", "validation", "test")
    }
    data_hashes = {
        f"{split}.fasta": sha256_file(data_dir / f"{split}.fasta")
        for split in datasets
    }

    training_config = config["model"]["training"]
    epochs = int(maximum_epochs or training_config["epochs"])
    if epochs <= 0:
        raise ValueError("Training epochs must be positive")
    batch_size = int(training_config["batch_size"])
    num_workers = int(training_config.get("num_workers", 0))
    gradient_clip_norm = float(training_config.get("gradient_clip_norm", 1.0))
    patience = int(training_config["early_stopping_patience"])

    model = build_sequence_vae(config, vocabulary).to(device)
    optimizer = Adam(model.parameters(), lr=float(training_config["learning_rate"]))
    history: list[dict[str, Any]] = []
    start_epoch = 1
    best_epoch = 0
    best_validation_loss = float("inf")
    epochs_without_improvement = 0

    if output_dir is None:
        model_dir = root / output_config["model_dir"]
    else:
        model_dir = Path(output_dir).resolve()
    model_dir.mkdir(parents=True, exist_ok=True)
    best_checkpoint = model_dir / "best.pt"
    last_checkpoint = model_dir / "last.pt"
    history_path = model_dir / "training_history.tsv"
    curve_path = model_dir / "training_curve.png"
    manifest_path = model_dir / "training_manifest.json"

    if resume_checkpoint is not None:
        resume_path = Path(resume_checkpoint)
        if not resume_path.is_absolute():
            resume_path = root / resume_path
        checkpoint = _load_torch_file(resume_path, device)
        if checkpoint["config_sha256"] != config_sha256:
            raise ValueError("Resume checkpoint configuration does not match")
        if checkpoint["model_config"] != model.model_config():
            raise ValueError("Resume checkpoint model configuration does not match")
        if checkpoint["vocabulary_sha256"] != sha256_file(vocabulary_path):
            raise ValueError("Resume checkpoint vocabulary does not match")
        if checkpoint["data_sha256"] != data_hashes:
            raise ValueError("Resume checkpoint training data does not match")
        model.load_state_dict(checkpoint["model_state_dict"], strict=True)
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        history = list(checkpoint["history"])
        start_epoch = int(checkpoint["epoch"]) + 1
        best_epoch = int(checkpoint["best_epoch"])
        best_validation_loss = float(checkpoint["best_validation_loss"])
        epochs_without_improvement = int(checkpoint["epochs_without_improvement"])
        _restore_rng_state(checkpoint["rng_state"])

    stopped_epoch = start_epoch - 1
    for epoch in range(start_epoch, epochs + 1):
        started = time.perf_counter()
        beta = kl_beta_for_epoch(
            epoch,
            target_beta=float(training_config["kl_beta"]),
            warmup_epochs=int(training_config["kl_warmup_epochs"]),
        )
        train_loader = _loader(
            datasets["train"],
            batch_size=batch_size,
            shuffle=True,
            pad_id=vocabulary.pad_id,
            seed=seed + epoch,
            num_workers=num_workers,
        )
        validation_loader = _loader(
            datasets["validation"],
            batch_size=batch_size,
            shuffle=False,
            pad_id=vocabulary.pad_id,
            seed=seed,
            num_workers=num_workers,
        )
        train_metrics = _run_epoch(
            model,
            train_loader,
            device=device,
            beta=beta,
            optimizer=optimizer,
            gradient_clip_norm=gradient_clip_norm,
        )
        validation_metrics = _run_epoch(
            model,
            validation_loader,
            device=device,
            # Keep the validation objective comparable across epochs while the
            # training objective follows the KL warm-up schedule.
            beta=float(training_config["kl_beta"]),
            optimizer=None,
            gradient_clip_norm=gradient_clip_norm,
        )
        row = {
            "epoch": epoch,
            "kl_beta": beta,
            "train_loss": train_metrics["loss"],
            "train_reconstruction": train_metrics["reconstruction"],
            "train_kl": train_metrics["kl"],
            "train_perplexity": train_metrics["perplexity"],
            "train_token_accuracy": train_metrics["token_accuracy"],
            "train_eos_accuracy": train_metrics["eos_accuracy"],
            "validation_loss": validation_metrics["loss"],
            "validation_reconstruction": validation_metrics["reconstruction"],
            "validation_kl": validation_metrics["kl"],
            "validation_perplexity": validation_metrics["perplexity"],
            "validation_token_accuracy": validation_metrics["token_accuracy"],
            "validation_eos_accuracy": validation_metrics["eos_accuracy"],
            "epoch_seconds": time.perf_counter() - started,
        }
        history.append(row)
        stopped_epoch = epoch
        improved = validation_metrics["loss"] < best_validation_loss
        if improved:
            best_epoch = epoch
            best_validation_loss = validation_metrics["loss"]
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1
        _save_checkpoint(
            last_checkpoint,
            model=model,
            optimizer=optimizer,
            epoch=epoch,
            best_epoch=best_epoch,
            best_validation_loss=best_validation_loss,
            epochs_without_improvement=epochs_without_improvement,
            history=history,
            config=config,
            config_sha256=config_sha256,
            vocabulary_path=vocabulary_path,
            data_hashes=data_hashes,
        )
        if improved:
            best_checkpoint.write_bytes(last_checkpoint.read_bytes())
        _write_history(history_path, history)
        if patience > 0 and epochs_without_improvement >= patience:
            break

    if not best_checkpoint.exists():
        raise RuntimeError("Training did not create a best checkpoint")
    best_state = _load_torch_file(best_checkpoint, device)
    model.load_state_dict(best_state["model_state_dict"], strict=True)
    test_metrics = _run_epoch(
        model,
        _loader(
            datasets["test"],
            batch_size=batch_size,
            shuffle=False,
            pad_id=vocabulary.pad_id,
            seed=seed,
            num_workers=num_workers,
        ),
        device=device,
        beta=float(training_config["kl_beta"]),
        optimizer=None,
        gradient_clip_norm=gradient_clip_norm,
    )
    _write_training_curve(curve_path, history)
    manifest = {
        "format_version": 1,
        "status": "development_provisional" if "pending" in verification_status else "formal",
        "warning": (
            "This model used provisional nominal GvpA data and is not a formal V2 model."
            if "pending" in verification_status
            else None
        ),
        "seed": seed,
        "device": str(device),
        "model_config": model.model_config(),
        "training_config": training_config,
        "data_verification_status": verification_status,
        "split_sizes": {name: len(dataset) for name, dataset in datasets.items()},
        "input_sha256": {
            resolved_config_path.relative_to(root).as_posix(): config_sha256,
            vocabulary_path.relative_to(root).as_posix(): sha256_file(vocabulary_path),
            **{(data_dir / name).relative_to(root).as_posix(): digest for name, digest in data_hashes.items()},
        },
        "best_epoch": best_epoch,
        "stopped_epoch": stopped_epoch,
        "best_validation_loss": best_validation_loss,
        "test_metrics": test_metrics,
        "posterior_collapse_diagnostic": {
            "best_epoch_validation_kl": history[best_epoch - 1]["validation_kl"],
            "warning_threshold": float(training_config.get("posterior_collapse_kl_threshold", 0.01)),
            "warning": history[best_epoch - 1]["validation_kl"]
            < float(training_config.get("posterior_collapse_kl_threshold", 0.01)),
        },
        "dependencies": _dependency_versions(),
        "artifacts_sha256": {
            "best.pt": sha256_file(best_checkpoint),
            "last.pt": sha256_file(last_checkpoint),
            "training_history.tsv": sha256_file(history_path),
            "training_curve.png": sha256_file(curve_path),
        },
    }
    write_json(manifest_path, manifest)
    return TrainingArtifacts(
        best_checkpoint=best_checkpoint,
        last_checkpoint=last_checkpoint,
        history_path=history_path,
        curve_path=curve_path,
        manifest_path=manifest_path,
        best_epoch=best_epoch,
        stopped_epoch=stopped_epoch,
    )


def generate_candidates(
    root: str | Path,
    config_path: str | Path,
    checkpoint_path: str | Path,
    *,
    device_name: str = "auto",
    allow_provisional_data: bool = False,
    candidate_count: int | None = None,
    output_dir: str | Path | None = None,
) -> GenerationArtifacts:
    root = Path(root).resolve()
    resolved_config_path, config = load_generation_config(root, config_path)
    device = resolve_device(device_name)
    seed = int(config["seed"])
    seed_everything(seed)
    output_config = config["outputs"]
    vocabulary_path = root / output_config["vocabulary"]
    vocabulary = ProteinVocabulary.from_json(vocabulary_path)
    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.is_absolute():
        checkpoint_path = root / checkpoint_path
    checkpoint = _load_torch_file(checkpoint_path, device)
    if checkpoint["config_sha256"] != sha256_file(resolved_config_path):
        raise ValueError("Checkpoint configuration hash does not match current config")
    if checkpoint["vocabulary_sha256"] != sha256_file(vocabulary_path):
        raise ValueError("Checkpoint vocabulary hash does not match configured vocabulary")
    checkpoint_status = str(checkpoint["config"]["data"]["verification_status"])
    if "pending" in checkpoint_status and not allow_provisional_data:
        raise RuntimeError(
            "Checkpoint used provisional data; pass allow_provisional_data only for development output"
        )
    model = SequenceVAE(**checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.eval()

    generation_config = config["generation"]
    count = int(candidate_count or generation_config["candidate_count"])
    if count <= 0:
        raise ValueError("candidate_count must be positive")
    batch_size = int(generation_config.get("batch_size", 64))
    maximum_length = int(config["model"]["maximum_sequence_length"])
    minimum_length = int(generation_config.get("minimum_sequence_length", 1))
    temperature = float(generation_config["sampling_temperature"])
    top_k = int(generation_config["top_k"])
    top_p = float(generation_config["top_p"])
    generator = torch.Generator(device=str(device))
    generator.manual_seed(seed)

    checkpoint_digest = sha256_file(checkpoint_path)
    records: list[FastaRecord] = []
    metadata: list[dict[str, object]] = []
    for start in range(0, count, batch_size):
        current_size = min(batch_size, count - start)
        generated = model.generate(
            num_samples=current_size,
            maximum_sequence_length=maximum_length,
            minimum_sequence_length=minimum_length,
            temperature=temperature,
            top_k=top_k,
            top_p=top_p,
            generator=generator,
        )
        for offset, result in enumerate(generated):
            latent_id = start + offset
            sequence_id = f"vae_candidate_{latent_id + 1:06d}"
            sequence = vocabulary.decode(result.token_ids)
            records.append(
                FastaRecord(
                    sequence_id,
                    f"generator=sequence_vae seed={seed} latent_id={latent_id}",
                    sequence,
                )
            )
            metadata.append(
                {
                    "sequence_id": sequence_id,
                    "generator_type": "sequence_vae",
                    "checkpoint_sha256": checkpoint_digest,
                    "generation_seed": seed,
                    "temperature": temperature,
                    "top_k": top_k,
                    "top_p": top_p,
                    "latent_id": latent_id,
                    "terminated_by_eos": str(result.terminated_by_eos).lower(),
                    "hit_generation_cap": str(result.hit_generation_cap).lower(),
                    "sequence_length": len(sequence),
                    "sequence_sha256": sequence_sha256(sequence),
                }
            )

    if output_dir is None:
        destination = root / output_config["candidate_output_dir"]
    else:
        destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    candidates_fasta = destination / "vae_candidates.fasta"
    metadata_path = destination / "generation_metadata.tsv"
    manifest_path = destination / "generation_manifest.json"
    write_fasta(records, candidates_fasta)
    write_tsv(metadata_path, metadata, list(metadata[0]))

    sequences = [record.sequence for record in records]
    sequence_counts = Counter(sequences)
    lengths = [len(sequence) for sequence in sequences]
    training_path = root / output_config["data_dir"] / "train.fasta"
    if not training_path.exists():
        raise FileNotFoundError(f"Configured training FASTA is missing: {training_path}")
    training_sequences = {record.sequence for record in read_fasta(training_path)}
    reference_sequences: set[str] | None = None
    reference_path_value = config["data"].get("reference_fasta")
    reference_path: Path | None = None
    if reference_path_value:
        reference_path = root / reference_path_value
        if not reference_path.exists():
            raise FileNotFoundError(f"Configured reference FASTA is missing: {reference_path}")
        reference_sequences = {record.sequence for record in read_fasta(reference_path)}
    diagnostic_input_sha256 = {
        training_path.relative_to(root).as_posix(): sha256_file(training_path)
    }
    if reference_path is not None:
        diagnostic_input_sha256[reference_path.relative_to(root).as_posix()] = sha256_file(
            reference_path
        )
    generation_diagnostics = {
        "unique_sequence_count": len(sequence_counts),
        "unique_sequence_rate": len(sequence_counts) / count,
        "duplicate_candidate_count": sum(
            occurrence_count - 1 for occurrence_count in sequence_counts.values()
        ),
        "duplicate_sequence_groups": sum(
            occurrence_count > 1 for occurrence_count in sequence_counts.values()
        ),
        "empty_sequence_count": sum(not sequence for sequence in sequences),
        "minimum_length": min(lengths),
        "mean_length": sum(lengths) / count,
        "maximum_length": max(lengths),
        "nonstandard_sequence_count": sum(
            bool(set(sequence) - set(vocabulary.amino_acids)) for sequence in sequences
        ),
        "exact_training_match_count": sum(
            sequence in training_sequences for sequence in sequences
        ),
        "reference_match_checked": reference_sequences is not None,
        "exact_reference_match_count": (
            sum(sequence in reference_sequences for sequence in sequences)
            if reference_sequences is not None
            else None
        ),
    }
    write_json(
        manifest_path,
        {
            "format_version": 1,
            "status": "development_provisional" if "pending" in checkpoint_status else "formal",
            "candidate_count": count,
            "seed": seed,
            "device": str(device),
            "generation_parameters": {
                "minimum_sequence_length": minimum_length,
                "maximum_sequence_length": maximum_length,
                "temperature": temperature,
                "top_k": top_k,
                "top_p": top_p,
                "batch_size": batch_size,
            },
            "input_sha256": {
                resolved_config_path.relative_to(root).as_posix(): sha256_file(resolved_config_path),
                vocabulary_path.relative_to(root).as_posix(): sha256_file(vocabulary_path),
                str(checkpoint_path): checkpoint_digest,
            },
            "output_sha256": {
                candidates_fasta.name: sha256_file(candidates_fasta),
                metadata_path.name: sha256_file(metadata_path),
            },
            "diagnostic_input_sha256": diagnostic_input_sha256,
            "generation_diagnostics": generation_diagnostics,
            "terminated_by_eos": sum(row["terminated_by_eos"] == "true" for row in metadata),
            "hit_generation_cap": sum(row["hit_generation_cap"] == "true" for row in metadata),
            "dependencies": _dependency_versions(),
        },
    )
    return GenerationArtifacts(candidates_fasta, metadata_path, manifest_path, count)
