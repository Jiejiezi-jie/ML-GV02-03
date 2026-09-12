from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .io import write_json


@dataclass(frozen=True)
class DomainHit:
    target: str
    target_length: int
    model: str
    model_length: int
    sequence_evalue: float
    sequence_score: float
    domain_evalue: float
    domain_score: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int

    @property
    def model_coverage(self) -> float:
        return (self.hmm_to - self.hmm_from + 1) / self.model_length


@dataclass(frozen=True)
class BlastHit:
    query: str
    target: str
    pident: float
    alignment_length: int
    query_length: int
    target_length: int
    evalue: float
    bitscore: float

    @property
    def effective_identity(self) -> float:
        """Identity adjusted for unaligned tails, on a 0--1 scale."""

        denominator = max(self.query_length, self.target_length)
        return (self.pident / 100.0) * self.alignment_length / denominator

    @property
    def query_coverage(self) -> float:
        return self.alignment_length / self.query_length

    @property
    def target_coverage(self) -> float:
        return self.alignment_length / self.target_length


def require_tools(names: Iterable[str]) -> dict[str, str]:
    resolved: dict[str, str] = {}
    for name in names:
        path = shutil.which(name)
        if path is None:
            raise RuntimeError(f"Required executable is missing: {name}")
        resolved[name] = path
    return resolved


def tool_version(name: str) -> str:
    executable = Path(name).name
    commands = {
        "hmmsearch": [name, "-h"],
        "mafft": [name, "--version"],
        "blastp": [name, "-version"],
        "makeblastdb": [name, "-version"],
        "cd-hit": [name, "-h"],
    }
    result = subprocess.run(
        commands.get(executable, [name, "--version"]),
        check=False,
        capture_output=True,
        text=True,
    )
    text = (result.stdout or result.stderr).strip().splitlines()
    return text[0].strip() if text else "unknown"


def run_command(
    args: list[str],
    *,
    stdout_path: str | Path | None = None,
    stderr_path: str | Path | None = None,
) -> None:
    stdout_stream = None
    stderr_stream = None
    try:
        if stdout_path is not None:
            stdout_path = Path(stdout_path)
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_stream = stdout_path.open("w", encoding="utf-8", newline="\n")
        if stderr_path is not None:
            stderr_path = Path(stderr_path)
            stderr_path.parent.mkdir(parents=True, exist_ok=True)
            stderr_stream = stderr_path.open("w", encoding="utf-8", newline="\n")
        subprocess.run(
            args,
            check=True,
            text=True,
            stdout=stdout_stream,
            stderr=stderr_stream,
        )
    finally:
        if stdout_stream is not None:
            stdout_stream.close()
        if stderr_stream is not None:
            stderr_stream.close()


def download_pfam_hmm(url: str, accession: str, target: str | Path) -> dict:
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        text = target.read_text(encoding="utf-8")
        if f"ACC   {accession}" not in text:
            raise ValueError(f"Existing HMM does not contain {accession}: {target}")
        return {
            "accession": accession,
            "url": url,
            "path": target.as_posix(),
            "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
            "cached": True,
        }
    request = urllib.request.Request(url, headers={"User-Agent": "GV02-03-course-project/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = response.read()
        content_type = response.headers.get("Content-Type", "")
    if payload.startswith(b"\x1f\x8b") or "gzip" in content_type:
        payload = gzip.decompress(payload)
    text = payload.decode("utf-8")
    if not text.startswith("HMMER3/") or f"ACC   {accession}" not in text:
        raise ValueError("Downloaded payload is not the requested HMMER profile")
    target.write_text(text, encoding="utf-8", newline="\n")
    metadata = {
        "accession": accession,
        "url": url,
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "path": target.as_posix(),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "cached": False,
    }
    write_json(target.with_suffix(target.suffix + ".metadata.json"), metadata)
    return metadata


def parse_hmm_metadata(path: str | Path) -> dict:
    metadata: dict[str, object] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line == "//":
            break
        if line.startswith("NAME"):
            metadata["name"] = line.split(maxsplit=1)[1]
        elif line.startswith("ACC"):
            metadata["accession"] = line.split(maxsplit=1)[1]
        elif line.startswith("LENG"):
            metadata["length"] = int(line.split()[1])
        elif line.startswith("GA"):
            values = line.replace(";", "").split()[1:]
            metadata["ga_sequence"] = float(values[0])
            metadata["ga_domain"] = float(values[1] if len(values) > 1 else values[0])
    for required in ("name", "accession", "length"):
        if required not in metadata:
            raise ValueError(f"HMM metadata is missing {required}: {path}")
    return metadata


def parse_domtbl(path: str | Path) -> dict[str, DomainHit]:
    best: dict[str, DomainHit] = {}
    with Path(path).open(encoding="utf-8") as stream:
        for raw in stream:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError(f"Malformed HMMER domtblout row: {raw.rstrip()}")
            hit = DomainHit(
                target=fields[0],
                target_length=int(fields[2]),
                model=fields[3],
                model_length=int(fields[5]),
                sequence_evalue=float(fields[6]),
                sequence_score=float(fields[7]),
                domain_evalue=float(fields[12]),
                domain_score=float(fields[13]),
                hmm_from=int(fields[15]),
                hmm_to=int(fields[16]),
                ali_from=int(fields[17]),
                ali_to=int(fields[18]),
            )
            previous = best.get(hit.target)
            if previous is None or (hit.domain_score, -hit.domain_evalue) > (
                previous.domain_score,
                -previous.domain_evalue,
            ):
                best[hit.target] = hit
    return best


def parse_blast(path: str | Path) -> list[BlastHit]:
    hits: list[BlastHit] = []
    with Path(path).open(encoding="utf-8") as stream:
        for raw in stream:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 12:
                raise ValueError(f"Malformed BLAST row: {raw.rstrip()}")
            hits.append(
                BlastHit(
                    query=fields[0],
                    target=fields[1],
                    pident=float(fields[2]),
                    alignment_length=int(fields[3]),
                    query_length=int(fields[4]),
                    target_length=int(fields[5]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11]),
                )
            )
    return hits
