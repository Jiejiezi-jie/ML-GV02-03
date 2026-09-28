"""Local Pfam-A architecture review of the frozen B-to-C candidate handoff.

Pfam hits are warning evidence, not replacements for B's QC or pool labels.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import platform

from .c_tm import _load_pools, _load_records, ensure_empty_output
from .io import STANDARD_AA, sequence_sha256, sha256_file, write_json, write_tsv


PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_MD5 = "7ab3c4e215d0daaea3004e37c4e24f8a"
TARGET_ACCESSION = "PF00741"
SEQUENCE_FIELDS = ("sequence_id", "pool", "sequence_sha256", "length",
                   "pfam_hit_count", "pf00741_hit_segments", "non_target_hit_count",
                   "domain_architecture", "warning_reasons", "review_status")
HIT_FIELDS = ("sequence_id", "pool", "accession", "accession_base", "model_name",
              "env_start_1based", "env_end_1based", "bit_score", "i_evalue",
              "hmm_start_1based", "hmm_end_1based", "hmm_length", "hmm_span_fraction",
              "reference_sequences_with_accession")


@dataclass(frozen=True)
class DomainHit:
    role: str
    sequence_id: str
    accession: str
    model_name: str
    start: int
    end: int
    bit_score: float
    i_evalue: float
    hmm_start: int
    hmm_end: int
    hmm_length: int

    @property
    def accession_base(self) -> str:
        return self.accession.split(".", 1)[0]


def warning_reasons(hits: list[DomainHit], reference_max_target_segments: int) -> tuple[str, ...]:
    targets = [hit for hit in hits if hit.accession_base == TARGET_ACCESSION]
    reasons = set()
    if len(targets) > reference_max_target_segments:
        reasons.add("multiple_pf00741_segments")
    for hit in hits:
        if hit.accession_base == TARGET_ACCESSION:
            continue
        overlap = max((max(0, min(hit.end, target.end) - max(hit.start, target.start) + 1)
                       for target in targets), default=0)
        if overlap > 0:
            reasons.add("overlapping_alternative_family")
        else:
            reasons.add("non_target_extra_domain")
    return tuple(sorted(reasons))


def _file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scan_profiles(pfam: Path, records: dict, cpus: int,
                  expected_pfam_md5: str | None = None) -> tuple[list[DomainHit], int]:
    """Search all supplied proteins against every Pfam HMM at curated GA cutoffs."""
    if not 1 <= cpus <= 32:
        raise ValueError("cpus must be between 1 and 32")
    if expected_pfam_md5 is not None and _file_md5(pfam) != expected_pfam_md5.lower():
        raise ValueError("Pfam MD5 mismatch")
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    names = {}
    for role in ("candidate", "reference"):
        for index, record in enumerate(records.get(role, ()), start=1):
            if not record.sequence or set(record.sequence) - STANDARD_AA:
                raise ValueError(f"Invalid protein sequence: {role}/{record.identifier}")
            key = ("c" if role == "candidate" else "r") + f"{index:06d}"
            names[key] = (role, record.identifier)
            sequences.append(pyhmmer.easel.TextSequence(
                name=key.encode("ascii"), sequence=record.sequence
            ).digitize(alphabet))
    if not sequences:
        raise ValueError("No protein sequences supplied")

    results = []
    count = 0
    opener = gzip.open if pfam.suffix == ".gz" else open
    with opener(pfam, "rb") as compressed, pyhmmer.plan7.HMMFile(compressed) as models:
        for top_hits in pyhmmer.hmmsearch(models, sequences, cpus=cpus,
                                         bit_cutoffs="gathering"):
            count += 1
            accession = top_hits.query.accession or top_hits.query.name
            model_name = top_hits.query.name
            for hit in top_hits:
                if not hit.included:
                    continue
                role, identifier = names[hit.name]
                for domain in hit.domains:
                    if domain.included:
                        results.append(DomainHit(
                            role, identifier, accession, model_name,
                            domain.env_from, domain.env_to, float(domain.score),
                            float(domain.i_evalue), domain.alignment.hmm_from,
                            domain.alignment.hmm_to, domain.alignment.hmm_length,
                        ))
            if count % 2000 == 0:
                print(f"Scanned {count} Pfam profile HMMs", flush=True)
    results.sort(key=lambda h: (0 if h.role == "candidate" else 1,
                                h.sequence_id, h.start, h.end, h.accession))
    return results, count


def _sequence_row(record, pool: str, hits: list[DomainHit], reference_max: int) -> dict:
    reasons = warning_reasons(hits, reference_max)
    return {
        "sequence_id": record.identifier,
        "pool": pool,
        "sequence_sha256": sequence_sha256(record.sequence),
        "length": len(record.sequence),
        "pfam_hit_count": len(hits),
        "pf00741_hit_segments": sum(h.accession_base == TARGET_ACCESSION for h in hits),
        "non_target_hit_count": sum(h.accession_base != TARGET_ACCESSION for h in hits),
        "domain_architecture": ";".join(f"{h.accession_base}:{h.start}-{h.end}" for h in hits),
        "warning_reasons": ";".join(reasons),
        "review_status": ("manual_review" if reasons else
                          "no_pfam_hit" if not hits else "no_extra_domain_evidence"),
    }


def _hit_row(hit: DomainHit, pool: str, reference_counts: Counter) -> dict:
    return {
        "sequence_id": hit.sequence_id,
        "pool": pool,
        "accession": hit.accession,
        "accession_base": hit.accession_base,
        "model_name": hit.model_name,
        "env_start_1based": hit.start,
        "env_end_1based": hit.end,
        "bit_score": hit.bit_score,
        "i_evalue": hit.i_evalue,
        "hmm_start_1based": hit.hmm_start,
        "hmm_end_1based": hit.hmm_end,
        "hmm_length": hit.hmm_length,
        "hmm_span_fraction": (hit.hmm_end - hit.hmm_start + 1) / hit.hmm_length,
        "reference_sequences_with_accession": reference_counts[hit.accession_base],
    }


def run(config_path: Path, handoff_root: Path, pfam: Path, output: Path, cpus: int = 4,
        expected_pfam_md5: str = PFAM_MD5) -> dict:
    ensure_empty_output(output)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    records = _load_records(handoff_root, config)
    pools = _load_pools(handoff_root, config)
    if set(pools) != {record.identifier for record in records["candidate"]}:
        raise ValueError("Candidate/pool identifiers disagree")
    hits, model_count = scan_profiles(pfam, records, cpus, expected_pfam_md5)
    import pyhmmer

    by_record = defaultdict(list)
    for hit in hits:
        by_record[(hit.role, hit.sequence_id)].append(hit)
    reference_max = max(sum(hit.accession_base == TARGET_ACCESSION
                            for hit in by_record[("reference", record.identifier)])
                        for record in records["reference"])
    if reference_max == 0:
        raise ValueError("No PF00741 hit in frozen reference set; cannot calibrate architecture")
    reference_counts = Counter()
    for record in records["reference"]:
        reference_counts.update({h.accession_base for h in by_record[("reference", record.identifier)]})

    candidate_rows = [_sequence_row(record, pools[record.identifier],
                                    by_record[("candidate", record.identifier)], reference_max)
                      for record in records["candidate"]]
    reference_rows = [_sequence_row(record, "reference",
                                    by_record[("reference", record.identifier)], reference_max)
                      for record in records["reference"]]
    candidate_hits = [_hit_row(hit, pools[hit.sequence_id], reference_counts)
                      for hit in hits if hit.role == "candidate"]
    reference_hits = [_hit_row(hit, "reference", reference_counts)
                      for hit in hits if hit.role == "reference"]
    outputs = {
        "candidate_domains.tsv": (candidate_rows, SEQUENCE_FIELDS),
        "reference_domains.tsv": (reference_rows, SEQUENCE_FIELDS),
        "candidate_domain_hits.tsv": (candidate_hits, HIT_FIELDS),
        "reference_domain_hits.tsv": (reference_hits, HIT_FIELDS),
    }
    for name, (rows, fields) in outputs.items():
        write_tsv(output / name, rows, fields)
    summary = {
        "method": "Local PyHMMER hmmsearch against all Pfam-A profiles with curated GA sequence and domain thresholds",
        "pfam_source": PFAM_URL,
        "pfam_md5": expected_pfam_md5,
        "pfam_profile_count": model_count,
        "candidate_count": len(candidate_rows),
        "reference_count": len(reference_rows),
        "reference_max_pf00741_hit_segments": reference_max,
        "reference_non_target_accessions": dict(sorted((key, value) for key, value in reference_counts.items()
                                                      if key != TARGET_ACCESSION)),
        "candidate_review_status_counts": dict(Counter(row["review_status"] for row in candidate_rows)),
        "candidate_review_status_by_pool": {
            pool: dict(Counter(row["review_status"] for row in candidate_rows if row["pool"] == pool))
            for pool in sorted(set(pools.values()))
        },
        "warning_reason_counts": dict(Counter(reason for row in candidate_rows
                                              for reason in filter(None, row["warning_reasons"].split(";")))),
        "decision_rule": "Non-target Pfam hit with no overlap to PF00741, any overlapping alternative-family hit, or PF00741 hit-segment count exceeding reference maximum triggers manual review only; B pool/QC unchanged.",
        "limitations": "Multiple PF00741 hit segments may reflect repeat-like sequence or split/overlapping HMMER alignments, not proven extra domain copies. Pfam GA significance is family-specific, but absence of a hit does not establish absence of a domain. Reference set remains provisional.",
    }
    write_json(output / "summary.json", summary)
    manifest = {
        "input_sha256": {key: config["inputs"][key]["sha256"] for key in
                         ("candidates", "reference_fasta", "pool_audit")},
        "config_sha256": sha256_file(config_path),
        "pfam_source": PFAM_URL,
        "pfam_md5": _file_md5(pfam),
        "pfam_sha256": sha256_file(pfam),
        "pyhmmer_version": pyhmmer.__version__,
        "python_version": platform.python_version(),
        "cpus": cpus,
        "output_sha256": {name: sha256_file(output / name) for name in (*outputs, "summary.json")},
    }
    write_json(output / "manifest.json", manifest)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("config", "handoff-root", "pfam", "output"):
        parser.add_argument(f"--{arg}", type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=4)
    args = parser.parse_args()
    run(args.config.resolve(), args.handoff_root.resolve(), args.pfam.resolve(),
        args.output.resolve(), args.cpus)


if __name__ == "__main__":
    main()
