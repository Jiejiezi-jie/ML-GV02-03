from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

from .io import FastaRecord, STANDARD_AA


def assess_candidates(
    records: Iterable[FastaRecord],
    *,
    minimum_length: int,
    maximum_length: int,
    minimum_unique_residues: int,
    maximum_single_residue_fraction: float,
    generation_length_cap: int | None = None,
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

        reasons = [*failures]
        if generation_cap_warning:
            reasons.append("generation_cap_warning")
        result[record.identifier] = {
            "length": length,
            "standard_aa_only": standard_aa_only,
            "duplicate_candidate": duplicate_candidate,
            "generation_cap_warning": generation_cap_warning,
            "low_complexity_warning": low_complexity_warning,
            "qc_pass": not failures,
            "qc_reasons": ";".join(reasons),
        }
    return result


def passing_candidates(
    records: Iterable[FastaRecord], results: dict[str, dict[str, bool | int | str]]
) -> list[FastaRecord]:
    """Keep QC-passing records in their original deterministic order."""
    return [record for record in records if results[record.identifier]["qc_pass"]]
