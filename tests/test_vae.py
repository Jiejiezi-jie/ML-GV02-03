from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
import torch
import yaml

from gv_eval.generation import generate_candidates, train_sequence_vae
from gv_eval.vae import (
    ProteinBatchCollator,
    ProteinSequenceDataset,
    ProteinVocabulary,
    SequenceVAE,
    VAEOutput,
    compute_vae_loss,
    kl_beta_for_epoch,
)


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def vocabulary() -> ProteinVocabulary:
    return ProteinVocabulary.from_json(
        ROOT / "models/generator/sequence_vae/vocabulary.json"
    )


def small_model(vocabulary: ProteinVocabulary) -> SequenceVAE:
    torch.manual_seed(42)
    return SequenceVAE(
        vocab_size=len(vocabulary.tokens),
        pad_id=vocabulary.pad_id,
        bos_id=vocabulary.bos_id,
        eos_id=vocabulary.eos_id,
        unk_id=vocabulary.unk_id,
        embedding_dim=8,
        hidden_dim=12,
        latent_dim=4,
    )


def test_vocabulary_and_batch_use_distinct_bos_eos_and_pad(tmp_path, vocabulary):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">a\nACDE\n>b\nACDEFG\n", encoding="utf-8")
    dataset = ProteinSequenceDataset(fasta, vocabulary, maximum_sequence_length=10)
    batch = ProteinBatchCollator(vocabulary.pad_id)([dataset[0], dataset[1]])

    assert vocabulary.pad_id == 0
    assert vocabulary.bos_id == 1
    assert vocabulary.eos_id == 2
    assert batch["tokens"][0, 0].item() == vocabulary.bos_id
    assert batch["tokens"][0, 5].item() == vocabulary.eos_id
    assert batch["tokens"][0, 6].item() == vocabulary.pad_id
    assert vocabulary.decode(dataset[0]["tokens"].tolist()) == "ACDE"


def test_forward_targets_and_loss_learn_eos_but_ignore_pad(tmp_path, vocabulary):
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(">a\nACDE\n>b\nACDEFG\n", encoding="utf-8")
    dataset = ProteinSequenceDataset(fasta, vocabulary, maximum_sequence_length=10)
    batch = ProteinBatchCollator(vocabulary.pad_id)([dataset[0], dataset[1]])
    model = small_model(vocabulary)
    output = model(batch["tokens"], batch["lengths"], sample_latent=False)
    loss = compute_vae_loss(output, pad_id=vocabulary.pad_id, beta=0.1)

    assert output.logits.shape[:2] == output.targets.shape
    assert output.targets[0, 0].item() == vocabulary.token_to_id["A"]
    assert vocabulary.eos_id in output.targets[0].tolist()
    assert loss.token_count == (4 + 1) + (6 + 1)
    assert torch.isfinite(loss.total)
    loss.total.backward()
    assert model.output_projection.weight.grad is not None

    modified_logits = output.logits.detach().clone()
    pad_positions = output.targets.eq(vocabulary.pad_id)
    modified_logits[pad_positions] = 1000.0
    modified = VAEOutput(
        logits=modified_logits,
        targets=output.targets,
        mu=output.mu.detach(),
        logvar=output.logvar.detach(),
        latent=output.latent.detach(),
    )
    modified_loss = compute_vae_loss(modified, pad_id=vocabulary.pad_id, beta=0.1)
    assert modified_loss.reconstruction.item() == pytest.approx(
        loss.reconstruction.item(), abs=1e-7
    )


def test_generation_tracks_eos_per_sequence(monkeypatch, vocabulary):
    model = small_model(vocabulary)
    amino_acid_id = vocabulary.token_to_id["A"]
    scheduled = [
        torch.tensor([vocabulary.eos_id, amino_acid_id]),
        torch.tensor([amino_acid_id, vocabulary.eos_id]),
    ]

    def fake_sample(*args, **kwargs):
        return scheduled.pop(0).to(args[0].device)

    monkeypatch.setattr("gv_eval.vae.sample_from_logits", fake_sample)
    generator = torch.Generator().manual_seed(42)
    generated = model.generate(
        num_samples=2,
        maximum_sequence_length=5,
        minimum_sequence_length=0,
        temperature=1.0,
        top_k=0,
        top_p=1.0,
        generator=generator,
    )

    assert generated[0].terminated_by_eos is True
    assert generated[0].token_ids == ()
    assert generated[1].terminated_by_eos is True
    assert generated[1].token_ids == (amino_acid_id,)
    assert scheduled == []


def test_generation_is_reproducible_for_fixed_seed(vocabulary):
    model = small_model(vocabulary).eval()
    arguments = {
        "num_samples": 3,
        "maximum_sequence_length": 8,
        "minimum_sequence_length": 3,
        "temperature": 1.0,
        "top_k": 0,
        "top_p": 1.0,
    }
    first = model.generate(**arguments, generator=torch.Generator().manual_seed(7))
    second = model.generate(**arguments, generator=torch.Generator().manual_seed(7))
    assert first == second


def test_kl_warmup_reaches_frozen_beta():
    assert kl_beta_for_epoch(1, target_beta=0.1, warmup_epochs=20) == pytest.approx(0.005)
    assert kl_beta_for_epoch(20, target_beta=0.1, warmup_epochs=20) == pytest.approx(0.1)
    assert kl_beta_for_epoch(30, target_beta=0.1, warmup_epochs=20) == pytest.approx(0.1)


def _write_tiny_training_project(root: Path) -> Path:
    data_dir = root / "data"
    data_dir.mkdir(parents=True)
    (data_dir / "train.fasta").write_text(
        ">t1\nACDEFGHIK\n>t2\nLMNPQRSTV\n>t3\nACDFGHIKL\n>t4\nMNPQRSTVW\n",
        encoding="utf-8",
    )
    (data_dir / "validation.fasta").write_text(
        ">v1\nACDEFGHIL\n>v2\nLMNPQRSTW\n", encoding="utf-8"
    )
    (data_dir / "test.fasta").write_text(
        ">e1\nACDEFGHIM\n>e2\nLMNPQRSTY\n", encoding="utf-8"
    )
    vocabulary_source = ROOT / "models/generator/sequence_vae/vocabulary.json"
    vocabulary_path = root / "models/vocabulary.json"
    vocabulary_path.parent.mkdir(parents=True)
    vocabulary_path.write_text(vocabulary_source.read_text(encoding="utf-8"), encoding="utf-8")
    config = {
        "version": 1,
        "seed": 42,
        "model": {
            "type": "sequence_vae",
            "architecture": "gru",
            "embedding_dim": 8,
            "hidden_dim": 12,
            "latent_dim": 4,
            "encoder_layers": 1,
            "decoder_layers": 1,
            "dropout": 0.0,
            "maximum_sequence_length": 20,
            "training": {
                "batch_size": 2,
                "epochs": 2,
                "learning_rate": 0.001,
                "kl_beta": 0.1,
                "kl_warmup_epochs": 2,
                "early_stopping_patience": 5,
                "gradient_clip_norm": 1.0,
                "num_workers": 0,
                "posterior_collapse_kl_threshold": 0.01,
            },
        },
        "data": {"verification_status": "pending_competitive_gvpa_gvpj_verification"},
        "vocabulary": {
            "special_tokens": ["<PAD>", "<BOS>", "<EOS>", "<UNK>"],
            "amino_acids": "ACDEFGHIKLMNPQRSTVWY",
        },
        "generation": {
            "candidate_count": 3,
            "batch_size": 2,
            "minimum_sequence_length": 3,
            "sampling_temperature": 1.0,
            "top_k": 0,
            "top_p": 1.0,
        },
        "outputs": {
            "data_dir": "data",
            "vocabulary": "models/vocabulary.json",
            "model_dir": "models/vae",
            "candidate_output_dir": "generated",
        },
    }
    config_path = root / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def test_tiny_training_resume_and_generation_are_auditable(tmp_path):
    config_path = _write_tiny_training_project(tmp_path)
    with pytest.raises(RuntimeError, match="provisional"):
        train_sequence_vae(tmp_path, config_path, device_name="cpu", maximum_epochs=1)

    first = train_sequence_vae(
        tmp_path,
        config_path,
        device_name="cpu",
        allow_provisional_data=True,
        maximum_epochs=1,
    )
    assert first.best_checkpoint.exists()
    assert first.last_checkpoint.exists()
    assert first.stopped_epoch == 1

    resumed = train_sequence_vae(
        tmp_path,
        config_path,
        device_name="cpu",
        allow_provisional_data=True,
        resume_checkpoint=first.last_checkpoint,
        maximum_epochs=2,
    )
    assert resumed.stopped_epoch == 2
    assert len(resumed.history_path.read_text(encoding="utf-8").splitlines()) == 3

    generated = generate_candidates(
        tmp_path,
        config_path,
        resumed.best_checkpoint,
        device_name="cpu",
        allow_provisional_data=True,
        candidate_count=3,
    )
    with generated.metadata_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 3
    assert len({row["sequence_id"] for row in rows}) == 3
    assert all(3 <= int(row["sequence_length"]) <= 20 for row in rows)
    assert all(row["checkpoint_sha256"] for row in rows)
    manifest = json.loads(generated.manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "development_provisional"
    assert manifest["candidate_count"] == 3
