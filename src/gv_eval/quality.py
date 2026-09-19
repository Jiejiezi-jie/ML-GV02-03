from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
import json
from pathlib import Path

from .io import FastaRecord, STANDARD_AA, read_fasta


def _match_index(records: Iterable[FastaRecord] | None, source: str) -> dict[str, list[str]] | None:
    """Index validated full sequences; retain all source IDs in stable order."""
    if records is None:
        return None
    index: dict[str, list[str]] = {}
    identifiers: set[str] = set()
    for record in records:
        sequence = record.sequence.upper()
        if not sequence or set(sequence) - STANDARD_AA:
            raise ValueError(f"Invalid protein sequence in {source}: {record.identifier}")
        if not record.identifier or record.identifier in identifiers:
            raise ValueError(f"Empty or duplicate identifier in {source}: {record.identifier}")
        identifiers.add(record.identifier)
        index.setdefault(sequence, []).append(record.identifier)
    if not index:
        raise ValueError(f"Empty {source} match source")
    return {sequence: sorted(ids) for sequence, ids in index.items()}


def assess_candidate_files(
    records: Iterable[FastaRecord],
    *,
    training_fasta: str | Path | None = None,
    reference_fasta: str | Path | None = None,
    **options,
) -> dict[str, dict[str, bool | int | str]]:
    """Read optional comparison files, failing on missing or invalid sources."""
    training = list(read_fasta(training_fasta)) if training_fasta is not None else None
    reference = list(read_fasta(reference_fasta)) if reference_fasta is not None else None
    return assess_candidates(records, training_records=training, reference_records=reference, **options)


def assess_candidates(
    records: Iterable[FastaRecord],
    *,
    minimum_length: int,
    maximum_length: int,
    minimum_unique_residues: int,
    maximum_single_residue_fraction: float,
    generation_length_cap: int | None = None,
    training_records: Iterable[FastaRecord] | None = None,
    reference_records: Iterable[FastaRecord] | None = None,
    exact_training_match_action: str = "warn",
    exact_reference_match_action: str = "warn",
) -> dict[str, dict[str, bool | int | str]]:
    """Return auditable quality-control results for each candidate sequence."""
    if minimum_length < 1 or maximum_length < minimum_length:
        raise ValueError("Candidate length bounds must satisfy 1 <= minimum <= maximum")
    if minimum_unique_residues < 1:
        raise ValueError("minimum_unique_residues must be positive")
    if not 0 < maximum_single_residue_fraction <= 1:
        raise ValueError("maximum_single_residue_fraction must be in (0, 1]")
    if generation_length_cap is not None and generation_length_cap < 1:
        raise ValueError("generation_length_cap must be positive when provided")

    actions = {"training": exact_training_match_action, "reference": exact_reference_match_action}
    for source, action in actions.items():
        if action not in ("warn", "exclude"):
            raise ValueError(f"exact_{source}_match_action must be 'warn' or 'exclude'")
    indexes = {
        "training": _match_index(training_records, "training"),
        "reference": _match_index(reference_records, "reference"),
    }
    for source, index in indexes.items():
        if index is None and actions[source] == "exclude":
            raise ValueError(f"Cannot exclude {source} matches without a comparison source")

    rows = list(records)
    if len({record.identifier for record in rows}) != len(rows):
        raise ValueError("Candidate FASTA contains duplicate identifiers")
    sequence_counts = Counter(record.sequence.upper() for record in rows)
    result: dict[str, dict[str, bool | int | str]] = {}
    for record in rows:
        sequence = record.sequence.upper()
        length = len(sequence)
        standard_aa_only = bool(sequence) and not (set(sequence) - STANDARD_AA)
        duplicate_candidate = sequence_counts[sequence] > 1
        counts = Counter(sequence)
        maximum_fraction = max(counts.values(), default=0) / length if length else 0.0
        low_complexity_warning = bool(sequence) and (
            len(counts) < minimum_unique_residues
            or maximum_fraction > maximum_single_residue_fraction
        )
        generation_cap_warning = generation_length_cap is not None and length >= generation_length_cap

        failures: list[str] = []
        if not sequence:
            failures.append("empty_sequence")
        elif not standard_aa_only:
            failures.append("nonstandard_amino_acids")
        if length < minimum_length:
            failures.append("length_below_minimum")
        if length > maximum_length:
            failures.append("length_above_maximum")
        if duplicate_candidate:
            failures.append("duplicate_candidate")
        if low_complexity_warning:
            failures.append("low_complexity")

        match_fields: dict[str, bool | str] = {}
        warnings: list[str] = []
        for source, index in indexes.items():
            matched_ids = index.get(sequence, []) if index is not None else []
            match_fields[f"{source}_match_checked"] = index is not None
            match_fields[f"exact_{source}_match"] = bool(matched_ids)
            match_fields[f"{source}_match_ids"] = json.dumps(matched_ids, ensure_ascii=False)
            if matched_ids:
                reason = f"exact_{source}_match"
                if actions[source] == "exclude":
                    failures.append(reason)
                else:
                    warnings.append(f"{reason}_warning")

        reasons = [*failures, *warnings]
        if generation_cap_warning:
            reasons.append("generation_cap_warning")
        result[record.identifier] = {
            "length": length,
            "standard_aa_only": standard_aa_only,
            "duplicate_candidate": duplicate_candidate,
            "generation_cap_warning": generation_cap_warning,
            "low_complexity_warning": low_complexity_warning,
            **match_fields,
            "qc_pass": not failures,
            "qc_reasons": ";".join(reasons),
        }
    return result


def passing_candidates(
    records: Iterable[FastaRecord], results: dict[str, dict[str, bool | int | str]]
) -> list[FastaRecord]:
    """Keep QC-passing records in their original deterministic order."""
    return [record for record in records if results[record.identifier]["qc_pass"]]
