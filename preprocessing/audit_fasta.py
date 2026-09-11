"""Read-only FASTA inventory. Uses only the Python standard library."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
from statistics import mean, median


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path):
    header, parts = None, []
    with Path(path).open(encoding="utf-8-sig") as stream:
        for line_number, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError("{}:{}: sequence before FASTA header".format(path, line_number))
                parts.append("".join(line.split()))
    if header is not None:
        yield header, "".join(parts)


def audit(path, root):
    records = list(read_fasta(path))
    sequences = [sequence for _, sequence in records]
    lengths = [len(sequence) for sequence in sequences]
    ids = [header.split()[0] if header.split() else "" for header, _ in records]
    labels = Counter()
    for header, _ in records:
        matches = set(re.findall(r"\bgvp([acfgijklmnsuvwyz])\b", header, flags=re.I))
        labels.update({"Gvp" + suffix.upper() for suffix in matches} or {"unlabelled"})
    return {
        "path": path.relative_to(root).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "records": len(records),
        "unique_sequences": len(set(sequences)),
        "duplicate_sequence_records": len(sequences) - len(set(sequences)),
        "duplicate_id_records": len(ids) - len(set(ids)),
        "empty_sequences": sum(not sequence for sequence in sequences),
        "min_length": min(lengths) if lengths else None,
        "max_length": max(lengths) if lengths else None,
        "mean_length": mean(lengths) if lengths else None,
        "median_length": median(lengths) if lengths else None,
        "nonstandard_characters": sorted(set("".join(sequences)) - STANDARD_AA),
        "header_family_labels": dict(sorted(labels.items())),
        "partial_or_fragment_headers": sum(bool(re.search(r"\b(partial|fragment)\b", h, re.I)) for h, _ in records),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    inventory = [audit(path, root) for path in sorted((root / "data").rglob("*.fasta"))]
    print(json.dumps({"note": "Header labels are annotations, not verified family assignments. Sequences are not aligned, filtered or modified.", "files": inventory}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
