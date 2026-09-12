from __future__ import annotations

from gv_eval.external import parse_blast, parse_domtbl


def test_parse_blast_tabular(tmp_path):
    target = tmp_path / "blast.tsv"
    target.write_text("q\tt\t75\t80\t100\t90\t1\t80\t1\t80\t1e-20\t100\n", encoding="utf-8")
    hit = parse_blast(target)[0]
    assert hit.query == "q"
    assert hit.target == "t"
    assert abs(hit.effective_identity - 0.6) < 1e-12


def test_parse_domtbl_keeps_best_domain(tmp_path):
    target = tmp_path / "hits.domtbl"
    rows = [
        "target - 100 PF00741 PF00741.1 70 1e-5 20 0 1 1 1e-5 1e-5 15 0 1 60 10 69 8 71 0.9 first",
        "target - 100 PF00741 PF00741.1 70 1e-8 40 0 1 1 1e-8 1e-8 35 0 2 69 5 72 4 74 0.9 second",
    ]
    target.write_text("# header\n" + "\n".join(rows) + "\n", encoding="utf-8")
    hit = parse_domtbl(target)["target"]
    assert hit.domain_score == 35.0
    assert hit.hmm_from == 2
    assert hit.hmm_to == 69

