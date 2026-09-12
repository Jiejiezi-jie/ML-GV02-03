from __future__ import annotations

import numpy as np

from gv_eval.external import BlastHit
from gv_eval.metrics import (
    candidate_distance_matrix,
    conservation_scores,
    find_conserved_sites,
)


def test_effective_identity_penalizes_unaligned_tails():
    hit = BlastHit("q", "t", 100.0, 50, 100, 80, 1e-10, 70.0)
    assert hit.effective_identity == 0.5
    assert hit.query_coverage == 0.5
    assert hit.target_coverage == 0.625


def test_conserved_sites_are_defined_from_references_only():
    reference = {"r1": "AC-D", "r2": "ACED", "r3": "AC-D"}
    sites = find_conserved_sites(reference, maximum_gap_fraction=0.67, minimum_consensus_fraction=1.0)
    assert [(site.alignment_position, site.consensus) for site in sites] == [
        (1, "A"), (2, "C"), (3, "E"), (4, "D")
    ]
    scores = conservation_scores({**reference, "q": "AT-D"}, sites, ["q"])["q"]
    assert scores["conserved_sites_total"] == 4
    assert scores["conserved_sites_covered"] == 3
    assert scores["conserved_sites_matched"] == 2
    assert scores["conservation_score"] == 0.5


def test_candidate_distance_is_symmetric_and_has_zero_diagonal():
    hits = [BlastHit("a", "b", 80.0, 10, 10, 10, 1e-5, 30.0)]
    matrix, details = candidate_distance_matrix(["a", "b"], hits)
    np.testing.assert_allclose(matrix, [[0.0, 0.2], [0.2, 0.0]])
    assert abs(details["a"]["uniqueness_score"] - 0.2) < 1e-12
