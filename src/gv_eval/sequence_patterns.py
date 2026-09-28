"""Conservative, exact sequence-pattern alerts, not topology/domain predictions."""
from __future__ import annotations

from itertools import groupby


HYDROPHOBIC_RESIDUES = frozenset("AVILMFWY")


def validate_pattern_options(**options: int | None) -> None:
    optional = {"hydrophobic_run_minimum", "homopolymer_minimum", "tandem_repeat_minimum_length"}
    for name, value in options.items():
        if value is None and name in optional:
            continue
        lower = 1 if name in optional else 2
        if type(value) is not int or value < lower:
            raise ValueError(f"{name} must be an integer >= {lower}" +
                             (" or null" if name in optional else ""))


def hydrophobic_segments(sequence: str, minimum: int) -> list[dict[str, int]]:
    segments = []
    offset = 0
    for hydrophobic, group in groupby(sequence, key=lambda aa: aa in HYDROPHOBIC_RESIDUES):
        length = sum(1 for _ in group)
        if hydrophobic and length >= minimum:
            segments.append({"start": offset + 1, "end": offset + length, "length": length})
        offset += length
    return segments


def homopolymer_segments(sequence: str, minimum: int) -> list[dict[str, object]]:
    segments = []
    offset = 0
    for residue, group in groupby(sequence):
        length = sum(1 for _ in group)
        if length >= minimum:
            segments.append({"start": offset + 1, "end": offset + length,
                             "length": length, "residue": residue})
        offset += length
    return segments


def tandem_repeat_segments(
    sequence: str, minimum_length: int, minimum_copies: int, maximum_motif_length: int,
) -> list[dict[str, object]]:
    """Report primitive exact repeats; leftmost phase, complete copies only."""
    segments = []
    size = len(sequence)
    for period in range(2, min(maximum_motif_length, size // minimum_copies) + 1):
        for start in range(size - period * minimum_copies + 1):
            # Same periodic run at a shifted phase: retain its leftmost start.
            if start and sequence[start - 1] == sequence[start + period - 1]:
                continue
            motif = sequence[start:start + period]
            if any(period % divisor == 0 and motif == motif[:divisor] * (period // divisor)
                   for divisor in range(1, period)):
                continue
            end = start + period
            while sequence[end:end + period] == motif:
                end += period
            length = end - start
            copies = length // period
            if copies >= minimum_copies and length >= minimum_length:
                segments.append({"start": start + 1, "end": end, "length": length,
                                 "motif": motif, "copies": copies})
    return sorted(segments, key=lambda row: (row["start"], row["end"], row["motif"]))
