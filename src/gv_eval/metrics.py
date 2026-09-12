from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .external import BlastHit, DomainHit


@dataclass(frozen=True)
class ConservedSite:
    alignment_position: int
    consensus: str
    consensus_fraction: float
    gap_fraction: float


def find_conserved_sites(
    reference_alignment: dict[str, str],
    maximum_gap_fraction: float,
    minimum_consensus_fraction: float,
) -> list[ConservedSite]:
    if not reference_alignment:
        raise ValueError("Reference alignment is empty")
    rows = list(reference_alignment.values())
    width = len(rows[0])
    sites: list[ConservedSite] = []
    for column in range(width):
        residues = [row[column].upper() for row in rows]
        gap_fraction = sum(value in "-." for value in residues) / len(residues)
        nongap = [value for value in residues if value not in "-."]
        if not nongap or gap_fraction > maximum_gap_fraction:
            continue
        consensus, count = Counter(nongap).most_common(1)[0]
        consensus_fraction = count / len(nongap)
        if consensus_fraction >= minimum_consensus_fraction:
            sites.append(
                ConservedSite(
                    alignment_position=column + 1,
                    consensus=consensus,
                    consensus_fraction=consensus_fraction,
                    gap_fraction=gap_fraction,
                )
            )
    if not sites:
        raise ValueError("No conserved sites satisfy the configured thresholds")
    return sites


def conservation_scores(
    combined_alignment: dict[str, str], sites: Iterable[ConservedSite], ids: Iterable[str]
) -> dict[str, dict[str, float | int]]:
    sites = list(sites)
    result: dict[str, dict[str, float | int]] = {}
    for identifier in ids:
        aligned = combined_alignment.get(identifier)
        if aligned is None:
            raise ValueError(f"Candidate is missing from combined alignment: {identifier}")
        residues = [aligned[site.alignment_position - 1].upper() for site in sites]
        covered = sum(value not in "-." for value in residues)
        matches = sum(value == site.consensus for value, site in zip(residues, sites))
        result[identifier] = {
            "conserved_sites_total": len(sites),
            "conserved_sites_covered": covered,
            "conserved_sites_matched": matches,
            "conservation_coverage": covered / len(sites),
            "conservation_score": matches / len(sites),
        }
    return result


def empirical_reference_score(value: float, reference_values: Iterable[float]) -> float:
    values = np.sort(np.asarray(list(reference_values), dtype=float))
    if values.size == 0:
        raise ValueError("Reference score distribution is empty")
    if value < values[0]:
        return 0.0
    return float(np.searchsorted(values, value, side="right") / values.size)


def domain_metrics(
    candidate_ids: Iterable[str],
    candidate_hits: dict[str, DomainHit],
    reference_hits: dict[str, DomainHit],
    domain_ga: float | None,
    coverage_floor: float,
    coverage_ceiling: float,
) -> tuple[dict[str, dict[str, float | bool | str]], dict[str, float]]:
    reference_scores = [hit.domain_score for hit in reference_hits.values()]
    if not reference_scores:
        raise ValueError("The Pfam model did not hit any cleaned reference sequence")
    cutoff = float(domain_ga) if domain_ga is not None else float(np.quantile(reference_scores, 0.01))
    trusted_coverages = [
        hit.model_coverage for hit in reference_hits.values() if hit.domain_score >= cutoff
    ]
    if not trusted_coverages:
        raise ValueError("No reference sequence reaches the HMM score cutoff")
    coverage_cutoff = float(
        np.clip(np.quantile(trusted_coverages, 0.05), coverage_floor, coverage_ceiling)
    )
    result: dict[str, dict[str, float | bool | str]] = {}
    for identifier in candidate_ids:
        hit = candidate_hits.get(identifier)
        if hit is None:
            result[identifier] = {
                "domain_status": "not_detected",
                "domain_score": 0.0,
                "domain_evalue": 1000.0,
                "model_coverage": 0.0,
                "domain_pass": False,
                "constraint_score": 0.0,
            }
            continue
        percentile = empirical_reference_score(hit.domain_score, reference_scores)
        result[identifier] = {
            "domain_status": "detected",
            "domain_score": hit.domain_score,
            "domain_evalue": hit.domain_evalue,
            "sequence_evalue": hit.sequence_evalue,
            "sequence_score": hit.sequence_score,
            "model_coverage": hit.model_coverage,
            "domain_pass": hit.domain_score >= cutoff and hit.model_coverage >= coverage_cutoff,
            "constraint_score": percentile * min(1.0, hit.model_coverage / coverage_cutoff),
            "hmm_from": hit.hmm_from,
            "hmm_to": hit.hmm_to,
            "alignment_from": hit.ali_from,
            "alignment_to": hit.ali_to,
        }
    return result, {
        "domain_score_cutoff": cutoff,
        "model_coverage_cutoff": coverage_cutoff,
        "reference_score_min": float(np.min(reference_scores)),
        "reference_score_median": float(np.median(reference_scores)),
        "reference_score_max": float(np.max(reference_scores)),
        "reference_hits": len(reference_scores),
    }


def closest_reference_metrics(
    candidate_ids: Iterable[str], hits: Iterable[BlastHit]
) -> dict[str, dict[str, float | str]]:
    best: dict[str, BlastHit] = {}
    for hit in hits:
        previous = best.get(hit.query)
        key = (hit.effective_identity, hit.bitscore, -hit.evalue, hit.target)
        previous_key = (
            previous.effective_identity,
            previous.bitscore,
            -previous.evalue,
            previous.target,
        ) if previous is not None else None
        if previous_key is None or key > previous_key:
            best[hit.query] = hit
    result: dict[str, dict[str, float | str]] = {}
    for identifier in candidate_ids:
        hit = best.get(identifier)
        if hit is None:
            result[identifier] = {
                "closest_reference_id": "",
                "closest_reference_pident": 0.0,
                "closest_reference_effective_identity": 0.0,
                "closest_reference_query_coverage": 0.0,
                "closest_reference_target_coverage": 0.0,
                "closest_reference_evalue": 1000.0,
                "novelty_score": 1.0,
            }
        else:
            result[identifier] = {
                "closest_reference_id": hit.target,
                "closest_reference_pident": hit.pident / 100.0,
                "closest_reference_effective_identity": hit.effective_identity,
                "closest_reference_query_coverage": hit.query_coverage,
                "closest_reference_target_coverage": hit.target_coverage,
                "closest_reference_evalue": hit.evalue,
                "novelty_score": 1.0 - hit.effective_identity,
            }
    return result


def candidate_distance_matrix(
    candidate_ids: list[str], hits: Iterable[BlastHit]
) -> tuple[np.ndarray, dict[str, dict[str, float | str]]]:
    index = {identifier: position for position, identifier in enumerate(candidate_ids)}
    identity = np.zeros((len(candidate_ids), len(candidate_ids)), dtype=float)
    np.fill_diagonal(identity, 1.0)
    for hit in hits:
        if hit.query == hit.target or hit.query not in index or hit.target not in index:
            continue
        left, right = index[hit.query], index[hit.target]
        value = hit.effective_identity
        identity[left, right] = max(identity[left, right], value)
        identity[right, left] = max(identity[right, left], value)
    distance = 1.0 - identity
    details: dict[str, dict[str, float | str]] = {}
    for position, identifier in enumerate(candidate_ids):
        others = np.delete(distance[position], position)
        nearest_position = int(
            np.argmax(np.where(np.arange(len(candidate_ids)) == position, -1.0, identity[position]))
        )
        details[identifier] = {
            "uniqueness_score": float(np.mean(others)),
            "nearest_candidate_id": candidate_ids[nearest_position],
            "nearest_candidate_identity": float(identity[position, nearest_position]),
            "nearest_candidate_distance": float(distance[position, nearest_position]),
        }
    return distance, details


def subset_diversity(distance: np.ndarray, positions: list[int]) -> float:
    if len(positions) < 2:
        return float("nan")
    values = distance[np.ix_(positions, positions)]
    upper = values[np.triu_indices(len(positions), k=1)]
    return float(np.mean(upper))

