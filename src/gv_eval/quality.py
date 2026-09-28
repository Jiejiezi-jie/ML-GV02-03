from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
import json
import math
from pathlib import Path

from .io import FastaRecord, STANDARD_AA, read_fasta
from .sequence_patterns import (
    hydrophobic_segments, homopolymer_segments, tandem_repeat_segments,
    validate_pattern_options,
)


AA_ORDER = tuple(sorted(STANDARD_AA))


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


def _composition(sequence: str) -> tuple[float, ...]:
    counts = Counter(sequence)
    return tuple(counts[residue] / len(sequence) for residue in AA_ORDER)


def _jensen_shannon(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    midpoint = tuple((a + b) / 2.0 for a, b in zip(left, right, strict=True))

    def divergence(values: tuple[float, ...]) -> float:
        return sum(
            value * math.log2(value / middle)
            for value, middle in zip(values, midpoint, strict=True)
            if value > 0
        )

    return (divergence(left) + divergence(right)) / 2.0


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _composition_model(
    records: list[FastaRecord] | None,
    quantile: float | None,
) -> tuple[tuple[float, ...], float] | None:
    if quantile is None:
        return None
    if not 0 < quantile < 1:
        raise ValueError("composition_outlier_quantile must be in (0, 1)")
    if records is None:
        raise ValueError("Composition outlier detection requires reference_records")
    if len(records) < 2:
        raise ValueError("Composition outlier detection requires at least two reference records")

    profiles = [_composition(record.sequence.upper()) for record in records]
    mean_profile = tuple(
        sum(profile[index] for profile in profiles) / len(profiles)
        for index in range(len(AA_ORDER))
    )
    reference_divergences = [_jensen_shannon(profile, mean_profile) for profile in profiles]
    return mean_profile, _quantile(reference_divergences, quantile)


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
    composition_outlier_quantile: float | None = None,
    composition_outlier_action: str = "warn",
    hydrophobic_run_minimum: int | None = None,
    homopolymer_minimum: int | None = None,
    tandem_repeat_minimum_length: int | None = None,
    tandem_repeat_minimum_copies: int = 3,
    tandem_repeat_maximum_motif_length: int = 6,
) -> dict[str, dict[str, object]]:
    """Return auditable quality-control results for each candidate sequence."""
    validate_pattern_options(
        hydrophobic_run_minimum=hydrophobic_run_minimum,
        homopolymer_minimum=homopolymer_minimum,
        tandem_repeat_minimum_length=tandem_repeat_minimum_length,
        tandem_repeat_minimum_copies=tandem_repeat_minimum_copies,
        tandem_repeat_maximum_motif_length=tandem_repeat_maximum_motif_length,
    )
    if minimum_length < 1 or maximum_length < minimum_length:
        raise ValueError("Candidate length bounds must satisfy 1 <= minimum <= maximum")
    if minimum_unique_residues < 1:
        raise ValueError("minimum_unique_residues must be positive")
    if not 0 < maximum_single_residue_fraction <= 1:
        raise ValueError("maximum_single_residue_fraction must be in (0, 1]")
    if generation_length_cap is not None and generation_length_cap < 1:
        raise ValueError("generation_length_cap must be positive when provided")

    training_rows = list(training_records) if training_records is not None else None
    reference_rows = list(reference_records) if reference_records is not None else None
    actions = {"training": exact_training_match_action, "reference": exact_reference_match_action}
    for source, action in actions.items():
        if action not in ("warn", "exclude"):
            raise ValueError(f"exact_{source}_match_action must be 'warn' or 'exclude'")
    if composition_outlier_action not in ("warn", "exclude"):
        raise ValueError("composition_outlier_action must be 'warn' or 'exclude'")
    indexes = {
        "training": _match_index(training_rows, "training"),
        "reference": _match_index(reference_rows, "reference"),
    }
    for source, index in indexes.items():
        if index is None and actions[source] == "exclude":
            raise ValueError(f"Cannot exclude {source} matches without a comparison source")
    composition_model = _composition_model(reference_rows, composition_outlier_quantile)

    rows = list(records)
    if len({record.identifier for record in rows}) != len(rows):
        raise ValueError("Candidate FASTA contains duplicate identifiers")
    sequence_counts = Counter(record.sequence.upper() for record in rows)
    result: dict[str, dict[str, object]] = {}
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

        match_fields: dict[str, object] = {}
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

        composition_checked = composition_model is not None and standard_aa_only
        composition_divergence: float | None = None
        composition_threshold: float | None = None
        composition_outlier_warning = False
        if composition_checked:
            mean_profile, composition_threshold = composition_model
            composition_divergence = _jensen_shannon(_composition(sequence), mean_profile)
            composition_outlier_warning = composition_divergence > composition_threshold
            if composition_outlier_warning:
                if composition_outlier_action == "exclude":
                    failures.append("composition_outlier")
                else:
                    warnings.append("composition_outlier_warning")

        if generation_cap_warning:
            warnings.append("generation_cap_warning")
        pattern_fields: dict[str, object] = {}
        checks = (
            ("hydrophobic", hydrophobic_run_minimum, hydrophobic_segments,
             "hydrophobic_run_warning"),
            ("homopolymer", homopolymer_minimum, homopolymer_segments,
             "homopolymer_warning"),
            ("tandem_repeat", tandem_repeat_minimum_length,
             lambda seq, minimum: tandem_repeat_segments(
                 seq, minimum, tandem_repeat_minimum_copies, tandem_repeat_maximum_motif_length),
             "tandem_repeat_warning"),
        )
        for name, minimum, detector, warning in checks:
            checked = minimum is not None and standard_aa_only
            segments = detector(sequence, minimum) if checked else []
            pattern_fields[f"{name}_checked"] = checked
            pattern_fields[f"{name}_warning"] = bool(segments)
            pattern_fields[f"{name}_segments"] = json.dumps(segments, ensure_ascii=False)
            if segments:
                warnings.append(warning)
        reasons = [*failures, *warnings]
        result[record.identifier] = {
            "length": length,
            "standard_aa_only": standard_aa_only,
            "duplicate_candidate": duplicate_candidate,
            "generation_cap_warning": generation_cap_warning,
            "low_complexity_warning": low_complexity_warning,
            "composition_outlier_checked": composition_checked,
            "composition_js_divergence": composition_divergence,
            "composition_js_threshold": composition_threshold,
            "composition_outlier_warning": composition_outlier_warning,
            **match_fields,
            **pattern_fields,
            "qc_pass": not failures,
            "qc_reasons": ";".join(reasons),
            "qc_warnings": ";".join(warnings),
        }
    return result


def passing_candidates(
    records: Iterable[FastaRecord], results: dict[str, dict[str, object]]
) -> list[FastaRecord]:
    """Keep QC-passing records in their original deterministic order."""
    return [record for record in records if results[record.identifier]["qc_pass"]]
