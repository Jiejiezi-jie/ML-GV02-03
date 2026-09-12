from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable

from .io import FastaRecord, STANDARD_AA, header_has_label, read_fasta, sequence_sha256


def build_clean_reference(
    sources: Iterable[str | Path],
    *,
    required_label: str,
    excluded_labels: Iterable[str],
    excluded_header_terms: Iterable[str],
    minimum_length: int,
    maximum_length: int,
) -> tuple[list[FastaRecord], list[dict], dict]:
    """Build a traceable, exact-deduplicated nominal GvpA reference set.

    Inclusion requires an explicit GvpA label in the original header.  This keeps
    unlabeled and GvpJ records out of the reference used to define conservation.
    Every source record is retained in the audit mapping, including duplicates.
    """

    excluded_labels = list(excluded_labels)
    excluded_terms = [term.lower() for term in excluded_header_terms]
    accepted_by_sequence: dict[str, FastaRecord] = {}
    canonical_id_by_sequence: dict[str, str] = {}
    mapping: list[dict] = []
    reason_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    total = 0

    for source in sources:
        source = Path(source)
        for record in read_fasta(source):
            total += 1
            header_lower = record.description.lower()
            reasons: list[str] = []
            if not header_has_label(record.description, required_label):
                reasons.append("missing_required_gvpa_label")
            if any(header_has_label(record.description, label) for label in excluded_labels):
                reasons.append("excluded_family_label")
            if any(term in header_lower for term in excluded_terms):
                reasons.append("excluded_header_term")
            invalid = sorted(set(record.sequence) - STANDARD_AA)
            if invalid:
                reasons.append("nonstandard_residue")
            if not minimum_length <= len(record.sequence) <= maximum_length:
                reasons.append("length_outside_range")

            accepted = not reasons
            duplicate = accepted and record.sequence in accepted_by_sequence
            if accepted and not duplicate:
                stable_id = f"ref_{len(accepted_by_sequence) + 1:06d}"
                clean = FastaRecord(stable_id, record.description, record.sequence)
                accepted_by_sequence[record.sequence] = clean
                canonical_id_by_sequence[record.sequence] = stable_id
                source_counts[source.as_posix()] += 1
            elif duplicate:
                reasons.append("exact_duplicate")

            for reason in reasons or ["accepted_unique"]:
                reason_counts[reason] += 1
            mapping.append(
                {
                    "reference_id": canonical_id_by_sequence.get(record.sequence, ""),
                    "source_path": source.as_posix(),
                    "source_id": record.identifier,
                    "source_header": record.description,
                    "sequence_sha256": sequence_sha256(record.sequence),
                    "length": len(record.sequence),
                    "accepted": accepted,
                    "unique_representative": accepted and not duplicate,
                    "status": ";".join(reasons) if reasons else "accepted_unique",
                }
            )

    clean_records = list(accepted_by_sequence.values())
    if not clean_records:
        raise ValueError("Reference cleaning produced no eligible sequence")
    summary = {
        "source_records": total,
        "accepted_source_records_including_duplicates": sum(
            row["accepted"] for row in mapping
        ),
        "clean_unique_sequences": len(clean_records),
        "reason_counts": dict(sorted(reason_counts.items())),
        "unique_sequences_first_seen_by_source": dict(sorted(source_counts.items())),
    }
    return clean_records, mapping, summary

