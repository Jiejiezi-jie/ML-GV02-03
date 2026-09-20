from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:
    path = Path(path)
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8-sig") as stream:
        for line_number, raw in enumerate(stream, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _record(header, parts, path, line_number)
                header, parts = line[1:].strip(), []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before header")
                parts.append("".join(line.split()).upper())
    if header is not None:
        yield _record(header, parts, path, None)


def _record(
    header: str, parts: Sequence[str], path: Path, line_number: int | None
) -> FastaRecord:
    sequence = "".join(parts)
    if not header:
        raise ValueError(f"{path}:{line_number or 'EOF'}: empty FASTA header")
    identifier = header.split()[0]
    if not identifier:
        raise ValueError(f"{path}:{line_number or 'EOF'}: empty FASTA identifier")
    return FastaRecord(identifier, header, sequence)


def write_fasta(
    records: Iterable[FastaRecord], path: str | Path, width: int = 80
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for record in records:
            stream.write(f">{record.identifier}")
            if record.description and record.description != record.identifier:
                stream.write(f" {record.description}")
            stream.write("\n")
            for start in range(0, len(record.sequence), width):
                stream.write(record.sequence[start : start + width] + "\n")


def read_alignment(path: str | Path) -> dict[str, str]:
    rows = list(read_fasta(path))
    if not rows:
        raise ValueError(f"Alignment is empty: {path}")
    lengths = {len(row.sequence) for row in rows}
    if len(lengths) != 1:
        raise ValueError(f"Alignment rows have unequal lengths: {path}")
    if len({row.identifier for row in rows}) != len(rows):
        raise ValueError(f"Alignment contains duplicate identifiers: {path}")
    return {row.identifier: row.sequence for row in rows}


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def write_json(path: str | Path, value: object) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_tsv(path: str | Path, rows: Iterable[dict], fieldnames: Sequence[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def header_has_label(header: str, label: str) -> bool:
    return bool(re.search(rf"(?i)(?:^|[^a-z0-9]){re.escape(label)}(?:[^a-z0-9]|$)", header))
