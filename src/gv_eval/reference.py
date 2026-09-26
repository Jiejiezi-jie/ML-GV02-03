from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable
import re
import json

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


def audit_family_sources(sources: Iterable[str | Path]):
    """Inventory local database annotations, never infer family from filenames.

    Return unique sequence records plus a row for EVERY source occurrence.
    Conflicting annotations (sequence or accession) quarantine all occurrences.
    These are provisional seeds; the caller must add independent domain and
    out-of-cluster evidence before calling a reference high confidence.
    """
    rows = []
    unique = {}
    sequence_labels = {}
    accession_labels = {}
    for source in sources:
        for index, record in enumerate(read_fasta(source), 1):
            digest = sequence_sha256(record.sequence)
            identifier = "nat_" + digest[:20]
            unique.setdefault(identifier, FastaRecord(identifier, "", record.sequence))
            labels = sorted(set(re.findall(
                r"(?i)(?<![a-z0-9])gvp[a-z](?![a-z0-9])", record.description.lower()
            )))
            accession = re.split(r"[|.]", record.identifier)[0]
            sequence_labels.setdefault(digest, set()).update(labels)
            accession_labels.setdefault(accession, set()).update(labels)
            reasons = []
            if not labels:
                reasons.append("missing_family_annotation")
            if not record.sequence or set(record.sequence) - STANDARD_AA:
                reasons.append("empty_or_nonstandard_sequence")
            if re.search(r"(?i)partial|fragment|predicted|rescued|probable", record.description):
                reasons.append("incomplete_or_provisional_annotation")
            family = labels[0] if len(labels) == 1 else ""
            bounds = (50, 180) if family == "gvpa" else (40, 1000)
            if not bounds[0] <= len(record.sequence) <= bounds[1]:
                reasons.append("length_outside_seed_range")
            rows.append(dict(reference_id=identifier, source_path=Path(source).as_posix(),
                             source_record_index=index, source_id=record.identifier,
                             source_accession=accession, source_header=record.description,
                             sequence_sha256=digest, length=len(record.sequence),
                             annotated_family=family, annotation_reasons=";".join(reasons)))
    for row in rows:
        labels = sequence_labels[row["sequence_sha256"]] | accession_labels[row["source_accession"]]
        if len(labels) > 1:
            row["annotation_reasons"] = ";".join(filter(None, [
                row["annotation_reasons"], "conflicting_family_annotation"]))
        row["seed_eligible"] = not row["annotation_reasons"]
    return list(unique.values()), rows


def prepare_reviewed_references(snapshot: str | Path, output: str | Path):
    """Extract explicit A/J recommended names from a pinned UniProt JSON snapshot.

    Reviewed does not mean experimentally proven. Preserve evidence codes,
    fragment/probable status and record version in a separate source manifest.
    """
    from .io import write_fasta
    payload = json.loads(Path(snapshot).read_text(encoding="utf-8"))
    records, rows = [], []
    for entry in payload["results"]:
        if entry["entryType"] != "UniProtKB reviewed (Swiss-Prot)":
            raise ValueError("Non-reviewed entry in curated snapshot")
        description = entry["proteinDescription"]
        name = description["recommendedName"]["fullName"]["value"]
        match = re.search(r"(?i)gas vesicle protein ([AJ])(?:\d+)?$", name)
        if not match:
            raise ValueError(f"Unrecognized reviewed family: {name}")
        family = "gvp" + match[1].lower()
        accession = entry["primaryAccession"]
        flag = description.get("flag", "")
        record = FastaRecord(accession, f"{family} reviewed {name} {flag}", entry["sequence"]["value"])
        records.append(record)
        rows.append(dict(accession=accession, family=family, recommended_name=name,
                         organism=entry["organism"]["scientificName"], flag=flag,
                         sequence_sha256=sequence_sha256(record.sequence),
                         evidence_codes=";".join(sorted(set(re.findall(r'ECO:\d+', json.dumps(entry))))),
                         entry_version=entry["entryAudit"]["entryVersion"],
                         last_updated=entry["entryAudit"]["lastAnnotationUpdateDate"],
                         source_url=f"https://www.uniprot.org/uniprotkb/{accession}/entry"))
    write_fasta(records, output)
    return records, rows

