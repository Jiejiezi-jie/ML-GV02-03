"""Conservative competitive HMM classification; PF00741 is not GvpA-specific."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

from .external import DomainHit
from .io import FastaRecord, STANDARD_AA, sequence_sha256


@dataclass(frozen=True)
class FamilyThresholds:
    pfam_sequence_score: float = 25.0
    pfam_domain_score: float = 25.0
    pfam_coverage: float = 0.8
    family_score: float = 25.0
    family_model_coverage: float = 0.6
    family_query_coverage: float = 0.6
    margin_bits: float = 10.0
    minimum_gvpa_length: int = 50
    maximum_gvpa_length: int = 180
    similarity_coverage: float = 0.70
    similarity_identity: float = 0.40
    similarity_margin: float = 0.05

    def __post_init__(self):
        for key, value in vars(self).items():
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"Non-finite threshold: {key}")
            if "coverage" in key and not 0 < value <= 1:
                raise ValueError(f"Coverage must be in (0, 1]: {key}")
        if self.margin_bits <= 0 or self.family_score < 0:
            raise ValueError("Family score must be nonnegative and margin positive")
        if not 0 < self.minimum_gvpa_length <= self.maximum_gvpa_length:
            raise ValueError("Invalid GvpA length bounds")
        if not 0 < self.similarity_identity <= 1 or not 0 < self.similarity_margin <= 1:
            raise ValueError("Similarity identity and margin must be in (0, 1]")


def pfam_pass(hit: DomainHit | None, thresholds: FamilyThresholds) -> bool:
    return bool(hit and hit.sequence_score >= thresholds.pfam_sequence_score
                and hit.domain_score >= thresholds.pfam_domain_score
                and hit.model_coverage >= thresholds.pfam_coverage)


def classify_family(record: FastaRecord, pfam_hit: DomainHit | None,
                    family_hits: Mapping[str, DomainHit], *,
                    available_profiles: set[str],
                    thresholds: FamilyThresholds | None = None) -> dict:
    """Classify one sequence using separately searched profiles.

    A missing profile is NOT equivalent to a searched profile with no hit.
    Missing runners-up have no invented score/margin. A reported short local
    motif cannot establish family identity; both coverages are mandatory.
    """
    t = thresholds or FamilyThresholds()
    passed = pfam_pass(pfam_hit, t)
    hits = sorted(((family, hit) for family, hit in family_hits.items()
                   if family in available_profiles),
                  key=lambda item: (-item[1].domain_score, item[0]))
    # Keep all reported scores for audit; only sufficiently covered alignments
    # can win, but poorly covered competitors can still make a winner ambiguous.
    covered = [(f, h) for f, h in hits
               if h.model_coverage >= t.family_model_coverage
               and (h.ali_to - h.ali_from + 1) / max(1, len(record.sequence)) >= t.family_query_coverage]
    best, best_hit = covered[0] if covered else ("", None)
    competitors = [(f, h) for f, h in hits if f != best]
    second, second_hit = competitors[0] if competitors else ("", None)
    margin = best_hit.domain_score - second_hit.domain_score if best_hit and second_hit else None
    row = dict(sequence_id=record.identifier, sequence_sha256=sequence_sha256(record.sequence),
               sequence_length=len(record.sequence), pf00741_status="pass" if passed else "fail",
               pf00741_score=pfam_hit.domain_score if pfam_hit else None,
               pf00741_sequence_score=pfam_hit.sequence_score if pfam_hit else None,
               pf00741_coverage=pfam_hit.model_coverage if pfam_hit else 0.0,
               best_family=best, best_family_score=best_hit.domain_score if best_hit else None,
               second_family=second, second_family_score=second_hit.domain_score if second_hit else None,
               family_margin=margin,
               best_model_coverage=best_hit.model_coverage if best_hit else 0.0,
               best_query_coverage=((best_hit.ali_to-best_hit.ali_from+1)/max(1,len(record.sequence))) if best_hit else 0.0,
               family_status="unsupported", classification_reason="no_qualified_subfamily_hit")
    for family in ("gvpa", "gvpj"):
        hit = family_hits.get(family)
        row[family + "_score"] = hit.domain_score if hit else None
    if not record.sequence or set(record.sequence) - STANDARD_AA:
        row["classification_reason"] = "empty_or_nonstandard_sequence"
    elif not {"gvpa", "gvpj"}.issubset(available_profiles):
        row["classification_reason"] = "missing_competitive_profile"
    elif best_hit and best_hit.domain_score >= t.family_score:
        if margin is not None and margin < t.margin_bits:
            row.update(family_status="gvpa_gvpj_ambiguous" if {best, second} <= {"gvpa", "gvpj"}
                       else "unsupported", classification_reason="insufficient_competitive_margin")
        elif best == "gvpa":
            if not passed:
                row["classification_reason"] = "gvpa_profile_but_pf00741_gate_failed"
            elif not t.minimum_gvpa_length <= len(record.sequence) <= t.maximum_gvpa_length:
                row["classification_reason"] = "gvpa_length_outside_range"
            else:
                row.update(family_status="supported_gvpa", classification_reason="pf00741_and_competitive_gvpa_support")
        else:
            row.update(family_status="other_gvp", classification_reason="competitive_" + best + "_support")
    if best_hit and second_hit is None and row["family_status"] != "unsupported":
        # hmmsearch is run at domain report threshold 0 bits. A non-hit's score
        # is NOT assumed to be zero; strong best score plus explicit no-hit is
        # the decision rule, and the numerical margin remains unavailable.
        row["classification_reason"] += ";no_reportable_runner_up"
    return row


def apply_similarity_evidence(row: dict, hits: list[dict], thresholds: FamilyThresholds) -> dict:
    """Require agreement with competitive, bidirectionally covered references.

    HMM profiles can conflate homologous A/J subfamilies. Strong profile support
    alone is insufficient. Sequence references are excluded by cluster during
    validation. Low identity/missing evidence causes abstention, never promotion.
    """
    row = dict(row)
    qualified = [h for h in hits if h["query_coverage"] >= thresholds.similarity_coverage
                 and h["target_coverage"] >= thresholds.similarity_coverage]
    best_by_family = {}
    for hit in sorted(qualified, key=lambda h: (-h["effective_identity"], h["target_id"])):
        best_by_family.setdefault(hit["family"], hit)
    ranked = list(best_by_family.values())
    best = ranked[0] if ranked else None
    second = ranked[1] if len(ranked) > 1 else None
    margin = best["effective_identity"] - second["effective_identity"] if second else None
    row.update(hmm_family_status=row["family_status"],
               similarity_best_family=best["family"] if best else "",
               similarity_best_reference=best["target_id"] if best else "",
               similarity_effective_identity=best["effective_identity"] if best else None,
               similarity_query_coverage=best["query_coverage"] if best else None,
               similarity_target_coverage=best["target_coverage"] if best else None,
               similarity_second_family=second["family"] if second else "",
               similarity_second_identity=second["effective_identity"] if second else None,
               similarity_margin=margin)
    if row["family_status"] not in ("supported_gvpa", "other_gvp"):
        return row
    reason = ""
    if not best or best["effective_identity"] < thresholds.similarity_identity:
        reason = "insufficient_covered_reference_similarity"
    elif best["family"] != row["best_family"]:
        reason = "profile_sequence_family_conflict"
    elif margin is not None and margin < thresholds.similarity_margin:
        reason = "insufficient_sequence_competitive_margin"
    if reason:
        families = {row["best_family"]}
        if best:
            families.add(best["family"])
        if second and margin < thresholds.similarity_margin:
            families.add(second["family"])
        row["family_status"] = ("gvpa_gvpj_ambiguous" if families == {"gvpa", "gvpj"}
                                else "unsupported")
        row["classification_reason"] += ";" + reason
    else:
        row["classification_reason"] += ";competitive_sequence_agreement"
    return row
