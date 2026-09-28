from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import torch
from torch import Tensor, nn
from torch.nn import functional as F
from torch.nn.utils.rnn import pack_padded_sequence
from torch.utils.data import Dataset

from .io import FastaRecord, read_fasta


PAD_TOKEN = "<PAD>"
BOS_TOKEN = "<BOS>"
EOS_TOKEN = "<EOS>"
UNK_TOKEN = "<UNK>"
REQUIRED_SPECIAL_TOKENS = (PAD_TOKEN, BOS_TOKEN, EOS_TOKEN, UNK_TOKEN)


@dataclass(frozen=True)
class ProteinVocabulary:
    """Auditable token mapping shared by training and generation."""

    tokens: tuple[str, ...]
    token_to_id: dict[str, int]
    special_tokens: tuple[str, ...]
    amino_acids: tuple[str, ...]

    @classmethod
    def from_json(cls, path: str | Path) -> "ProteinVocabulary":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        tokens = tuple(raw["tokens"])
        token_to_id = {str(token): int(index) for token, index in raw["token_to_id"].items()}
        special_tokens = tuple(raw["special_tokens"])
        amino_acids = tuple(raw["amino_acids"])
        if len(tokens) != len(set(tokens)):
            raise ValueError("Vocabulary contains duplicate tokens")
        if token_to_id != {token: index for index, token in enumerate(tokens)}:
            raise ValueError("Vocabulary token_to_id does not match token order")
        if tuple(special_tokens) != REQUIRED_SPECIAL_TOKENS:
            raise ValueError(
                f"Special token order must be {REQUIRED_SPECIAL_TOKENS}, got {special_tokens}"
            )
        if set(special_tokens) & set(amino_acids):
            raise ValueError("Special tokens and amino-acid tokens overlap")
        if set(tokens) != set(special_tokens) | set(amino_acids):
            raise ValueError("Vocabulary tokens do not match declared token groups")
        return cls(tokens, token_to_id, special_tokens, amino_acids)

    @property
    def pad_id(self) -> int:
        return self.token_to_id[PAD_TOKEN]

    @property
    def bos_id(self) -> int:
        return self.token_to_id[BOS_TOKEN]

    @property
    def eos_id(self) -> int:
        return self.token_to_id[EOS_TOKEN]

    @property
    def unk_id(self) -> int:
        return self.token_to_id[UNK_TOKEN]

    def encode(self, sequence: str, *, reject_unknown: bool = True) -> list[int]:
        sequence = sequence.strip().upper()
        unknown = sorted(set(sequence) - set(self.amino_acids))
        if reject_unknown and unknown:
            raise ValueError(f"Sequence contains tokens outside the vocabulary: {unknown}")
        encoded = [self.bos_id]
        encoded.extend(self.token_to_id.get(residue, self.unk_id) for residue in sequence)
        encoded.append(self.eos_id)
        return encoded

    def decode(self, token_ids: Sequence[int], *, stop_at_eos: bool = True) -> str:
        residues: list[str] = []
        for raw_token_id in token_ids:
            token_id = int(raw_token_id)
            if token_id < 0 or token_id >= len(self.tokens):
                raise ValueError(f"Token ID outside vocabulary: {token_id}")
            token = self.tokens[token_id]
            if token == EOS_TOKEN and stop_at_eos:
                break
            if token in (PAD_TOKEN, BOS_TOKEN, EOS_TOKEN):
                continue
            residues.append("X" if token == UNK_TOKEN else token)
        return "".join(residues)


class ProteinSequenceDataset(Dataset):
    """FASTA dataset whose samples include BOS/EOS and never include padding."""

    def __init__(
        self,
        fasta_path: str | Path,
        vocabulary: ProteinVocabulary,
        *,
        maximum_sequence_length: int,
    ) -> None:
        self.fasta_path = Path(fasta_path)
        self.vocabulary = vocabulary
        self.maximum_sequence_length = int(maximum_sequence_length)
        self.records = list(read_fasta(self.fasta_path))
        if not self.records:
            raise ValueError(f"Training FASTA is empty: {self.fasta_path}")
        identifiers = [record.identifier for record in self.records]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError(f"Training FASTA contains duplicate identifiers: {self.fasta_path}")
        for record in self.records:
            if not record.sequence:
                raise ValueError(f"Empty sequence: {record.identifier}")
            if len(record.sequence) > self.maximum_sequence_length:
                raise ValueError(
                    f"{record.identifier} has length {len(record.sequence)}, exceeding "
                    f"maximum_sequence_length={self.maximum_sequence_length}"
                )
            self.vocabulary.encode(record.sequence, reject_unknown=True)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> dict[str, object]:
        record: FastaRecord = self.records[index]
        token_ids = self.vocabulary.encode(record.sequence, reject_unknown=True)
        return {
            "sequence_id": record.identifier,
            "tokens": torch.tensor(token_ids, dtype=torch.long),
            "length": len(token_ids),
        }


class ProteinBatchCollator:
    def __init__(self, pad_id: int) -> None:
        self.pad_id = int(pad_id)

    def __call__(self, samples: Sequence[dict[str, object]]) -> dict[str, object]:
        if not samples:
            raise ValueError("Cannot collate an empty batch")
        lengths = torch.tensor([int(sample["length"]) for sample in samples], dtype=torch.long)
        tokens = torch.full(
            (len(samples), int(lengths.max().item())), self.pad_id, dtype=torch.long
        )
        for row, sample in enumerate(samples):
            sequence_tokens = sample["tokens"]
            if not isinstance(sequence_tokens, Tensor):
                raise TypeError("Dataset sample tokens must be a Tensor")
            tokens[row, : sequence_tokens.numel()] = sequence_tokens
        return {
            "sequence_ids": [str(sample["sequence_id"]) for sample in samples],
            "tokens": tokens,
            "lengths": lengths,
            "attention_mask": tokens.ne(self.pad_id),
        }


@dataclass(frozen=True)
class VAEOutput:
    logits: Tensor
    targets: Tensor
    mu: Tensor
    logvar: Tensor
    latent: Tensor


@dataclass(frozen=True)
class VAELoss:
    total: Tensor
    reconstruction: Tensor
    kl: Tensor
    token_count: int
    beta: float


@dataclass(frozen=True)
class GeneratedProtein:
    token_ids: tuple[int, ...]
    terminated_by_eos: bool
    hit_generation_cap: bool


class SequenceVAE(nn.Module):
    """Unconditional GRU sequence VAE for the provisional GvpA corpus.

    Padding-aware packed encoding prevents padded tokens from changing the latent
    representation. During generation each sequence tracks EOS independently.
    """

    def __init__(
        self,
        *,
        vocab_size: int,
        pad_id: int,
        bos_id: int,
        eos_id: int,
        unk_id: int,
        embedding_dim: int,
        hidden_dim: int,
        latent_dim: int,
        encoder_layers: int = 1,
        decoder_layers: int = 1,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        if min(vocab_size, embedding_dim, hidden_dim, latent_dim, encoder_layers, decoder_layers) <= 0:
            raise ValueError("Model dimensions and layer counts must be positive")
        if not 0.0 <= dropout < 1.0:
            raise ValueError("dropout must be in [0, 1)")
        self.vocab_size = int(vocab_size)
        self.pad_id = int(pad_id)
        self.bos_id = int(bos_id)
        self.eos_id = int(eos_id)
        self.unk_id = int(unk_id)
        self.embedding_dim = int(embedding_dim)
        self.hidden_dim = int(hidden_dim)
        self.latent_dim = int(latent_dim)
        self.encoder_layers = int(encoder_layers)
        self.decoder_layers = int(decoder_layers)
        self.dropout = float(dropout)

        self.embedding = nn.Embedding(vocab_size, embedding_dim, padding_idx=pad_id)
        self.encoder_gru = nn.GRU(
            embedding_dim,
            hidden_dim,
            num_layers=encoder_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if encoder_layers > 1 else 0.0,
        )
        self.mu_projection = nn.Linear(hidden_dim * 2, latent_dim)
        self.logvar_projection = nn.Linear(hidden_dim * 2, latent_dim)
        self.latent_to_hidden = nn.Linear(latent_dim, decoder_layers * hidden_dim)
        self.decoder_gru = nn.GRU(
            embedding_dim + latent_dim,
            hidden_dim,
            num_layers=decoder_layers,
            batch_first=True,
            dropout=dropout if decoder_layers > 1 else 0.0,
        )
        self.output_projection = nn.Linear(hidden_dim, vocab_size)

    def model_config(self) -> dict[str, int | float]:
        return {
            "vocab_size": self.vocab_size,
            "pad_id": self.pad_id,
            "bos_id": self.bos_id,
            "eos_id": self.eos_id,
            "unk_id": self.unk_id,
            "embedding_dim": self.embedding_dim,
            "hidden_dim": self.hidden_dim,
            "latent_dim": self.latent_dim,
            "encoder_layers": self.encoder_layers,
            "decoder_layers": self.decoder_layers,
            "dropout": self.dropout,
        }

    def encode(self, tokens: Tensor, lengths: Tensor) -> tuple[Tensor, Tensor]:
        if tokens.ndim != 2 or lengths.ndim != 1 or tokens.size(0) != lengths.size(0):
            raise ValueError("tokens must be [batch, length] and lengths must be [batch]")
        embedded = self.embedding(tokens)
        packed = pack_padded_sequence(
            embedded,
            lengths.detach().cpu(),
            batch_first=True,
            enforce_sorted=False,
        )
        _, hidden = self.encoder_gru(packed)
        hidden = hidden.view(self.encoder_layers, 2, tokens.size(0), self.hidden_dim)
        final = torch.cat((hidden[-1, 0], hidden[-1, 1]), dim=-1)
        return self.mu_projection(final), self.logvar_projection(final)

    @staticmethod
    def reparameterize(mu: Tensor, logvar: Tensor) -> Tensor:
        return mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)

    def _initial_decoder_hidden(self, latent: Tensor) -> Tensor:
        hidden = torch.tanh(self.latent_to_hidden(latent))
        return hidden.view(latent.size(0), self.decoder_layers, self.hidden_dim).transpose(0, 1).contiguous()

    def decode_teacher_forced(self, decoder_tokens: Tensor, latent: Tensor) -> Tensor:
        embedded = self.embedding(decoder_tokens)
        repeated_latent = latent.unsqueeze(1).expand(-1, decoder_tokens.size(1), -1)
        decoder_input = torch.cat((embedded, repeated_latent), dim=-1)
        decoded, _ = self.decoder_gru(decoder_input, self._initial_decoder_hidden(latent))
        return self.output_projection(decoded)

    def forward(
        self, tokens: Tensor, lengths: Tensor, *, sample_latent: bool = True
    ) -> VAEOutput:
        if tokens.size(1) < 2:
            raise ValueError("Each sequence must contain at least BOS and EOS")
        mu, logvar = self.encode(tokens, lengths)
        latent = self.reparameterize(mu, logvar) if sample_latent else mu
        decoder_tokens = tokens[:, :-1]
        targets = tokens[:, 1:]
        logits = self.decode_teacher_forced(decoder_tokens, latent)
        return VAEOutput(logits=logits, targets=targets, mu=mu, logvar=logvar, latent=latent)

    @torch.no_grad()
    def generate(
        self,
        *,
        num_samples: int,
        maximum_sequence_length: int,
        minimum_sequence_length: int,
        temperature: float,
        top_k: int,
        top_p: float,
        generator: torch.Generator,
    ) -> list[GeneratedProtein]:
        if num_samples <= 0:
            raise ValueError("num_samples must be positive")
        if not 0 <= minimum_sequence_length <= maximum_sequence_length:
            raise ValueError("Generation length bounds are invalid")
        device = next(self.parameters()).device
        latent = torch.randn(
            (num_samples, self.latent_dim), device=device, generator=generator
        )
        hidden = self._initial_decoder_hidden(latent)
        previous = torch.full(
            (num_samples, 1), self.bos_id, dtype=torch.long, device=device
        )
        active = torch.ones(num_samples, dtype=torch.bool, device=device)
        terminated = torch.zeros(num_samples, dtype=torch.bool, device=device)
        generated: list[list[int]] = [[] for _ in range(num_samples)]

        for position in range(maximum_sequence_length + 1):
            embedded = self.embedding(previous)
            decoder_input = torch.cat((embedded, latent.unsqueeze(1)), dim=-1)
            decoded, hidden = self.decoder_gru(decoder_input, hidden)
            logits = self.output_projection(decoded[:, -1, :])
            logits[:, [self.pad_id, self.bos_id, self.unk_id]] = -torch.inf
            if position < minimum_sequence_length:
                logits[:, self.eos_id] = -torch.inf
            if position == maximum_sequence_length:
                break
            next_token = sample_from_logits(
                logits,
                temperature=temperature,
                top_k=top_k,
                top_p=top_p,
                generator=generator,
            )
            next_token = torch.where(
                active, next_token, torch.full_like(next_token, self.eos_id)
            )
            for row in range(num_samples):
                if not bool(active[row]):
                    continue
                token_id = int(next_token[row].item())
                if token_id == self.eos_id:
                    terminated[row] = True
                    active[row] = False
                else:
                    generated[row].append(token_id)
            if not bool(active.any()):
                break
            previous = next_token.unsqueeze(1)

        return [
            GeneratedProtein(
                token_ids=tuple(token_ids),
                terminated_by_eos=bool(terminated[index].item()),
                hit_generation_cap=not bool(terminated[index].item()),
            )
            for index, token_ids in enumerate(generated)
        ]


def compute_vae_loss(output: VAEOutput, *, pad_id: int, beta: float) -> VAELoss:
    if beta < 0:
        raise ValueError("KL beta cannot be negative")
    flat_targets = output.targets.reshape(-1)
    token_count = int(flat_targets.ne(pad_id).sum().item())
    if token_count == 0:
        raise ValueError("Batch contains no reconstruction targets")
    reconstruction_sum = F.cross_entropy(
        output.logits.reshape(-1, output.logits.size(-1)),
        flat_targets,
        ignore_index=pad_id,
        reduction="sum",
    )
    reconstruction = reconstruction_sum / token_count
    kl = -0.5 * torch.sum(
        1 + output.logvar - output.mu.pow(2) - output.logvar.exp(), dim=1
    ).mean()
    return VAELoss(
        total=reconstruction + float(beta) * kl,
        reconstruction=reconstruction,
        kl=kl,
        token_count=token_count,
        beta=float(beta),
    )


def kl_beta_for_epoch(epoch: int, *, target_beta: float, warmup_epochs: int) -> float:
    if epoch < 1:
        raise ValueError("epoch is one-based and must be positive")
    if target_beta < 0 or warmup_epochs < 0:
        raise ValueError("KL parameters cannot be negative")
    if warmup_epochs == 0:
        return float(target_beta)
    return float(target_beta) * min(1.0, epoch / warmup_epochs)


def sample_from_logits(
    logits: Tensor,
    *,
    temperature: float,
    top_k: int,
    top_p: float,
    generator: torch.Generator,
) -> Tensor:
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    if top_k < 0:
        raise ValueError("top_k cannot be negative")
    if not 0 < top_p <= 1:
        raise ValueError("top_p must be in (0, 1]")
    filtered = logits / temperature
    if 0 < top_k < filtered.size(-1):
        threshold = torch.topk(filtered, top_k, dim=-1).values[:, -1:]
        filtered = filtered.masked_fill(filtered < threshold, -torch.inf)
    if top_p < 1.0:
        sorted_logits, sorted_indices = torch.sort(filtered, descending=True, dim=-1)
        cumulative = torch.softmax(sorted_logits, dim=-1).cumsum(dim=-1)
        remove = cumulative > top_p
        remove[:, 1:] = remove[:, :-1].clone()
        remove[:, 0] = False
        sorted_logits = sorted_logits.masked_fill(remove, -torch.inf)
        filtered = torch.full_like(filtered, -torch.inf).scatter(
            1, sorted_indices, sorted_logits
        )
    probabilities = torch.softmax(filtered, dim=-1)
    return torch.multinomial(probabilities, 1, generator=generator).squeeze(1)
