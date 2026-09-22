"""Evidence-bound adjudication of pinned NCBI records, never model-driven relabeling."""
from __future__ import annotations

import re
from .io import sequence_sha256


def parse_genbank_proteins(text):
    records = {}
    for block in text.split("//"):
        if not block.strip():
            continue
        version = re.search(r"^VERSION\s+(\S+)", block, re.M)
        definition = re.search(r"^DEFINITION\s+(.+?)(?=\nACCESSION)", block, re.M | re.S)
        origin = block.split("ORIGIN", 1)
        if not version or not definition or len(origin) != 2:
            raise ValueError("Incomplete GenBank protein record")
        sequence = re.sub(r"[^a-zA-Z]", "", origin[1]).upper()
        accession = version[1]
        if accession in records:
            raise ValueError(f"Duplicate NCBI accession: {accession}")
        def field(label):
            match = re.search(re.escape(label) + r"\s*::\s*([^\n]+)", block)
            return match[1].strip() if match else ""
        records[accession] = dict(accession=accession, definition=" ".join(definition[1].split()),
            sequence=sequence, sequence_sha256=sequence_sha256(sequence),
            name_evidence_category=field("Evidence Category"),
            name_evidence_accession=field("Evidence Accession"),
            name_evidence_source=field("Evidence Source"),
            source_identifier=field("Source Identifier"),
            gene_names=";".join(re.findall(r'/gene="([^"]+)"',block)),
            region_names=";".join(re.findall(r'/region_name="([^"]+)"',block)))
    return records


def adjudicate_conflict(row, official):
    if row["sequence_sha256"] != official["sequence_sha256"]:
        raise ValueError(f"NCBI/local sequence mismatch: {official['accession']}")
    if official["accession"] not in row["source_accessions"].split(";"):
        raise ValueError("Accession mismatch")
    name = official["definition"].lower()
    if re.search(r"\bgvpa\b", name) and "gvpa" in official["gene_names"].lower().split(";"):
        disposition = "resolved_stale_local_gvpj_annotation"
        status = "supported_gvpa"
        reason = "current_ncbi_gvpa_annotation_exact_sequence_match"
    elif re.search(r"\bgvpj\b", name) and official["source_identifier"].startswith("PF00741."):
        disposition = "adjudicated_ambiguous_not_gold_gvpj"
        status = "gvpa_gvpj_ambiguous"
        reason = "ncbi_gvpj_name_from_nonexclusive_pf00741_conflicts_with_competitive_gvpa_evidence"
    else:
        disposition = "unresolved_requires_additional_evidence"
        status = "unsupported"
        reason = "no_sufficient_official_adjudication_evidence"
    return dict(sequence_id=row["sequence_id"], sequence_sha256=row["sequence_sha256"],
        accession=official["accession"], original_local_family=row["expected_family"],
        original_prediction=row["family_status"], ncbi_definition=official["definition"],
        name_evidence_category=official["name_evidence_category"],
        name_evidence_accession=official["name_evidence_accession"],
        source_identifier=official["source_identifier"], region_names=official["region_names"],
        disposition=disposition, operational_family_status=status, reason=reason,
        biological_identity_resolved=disposition == "resolved_stale_local_gvpj_annotation",
        independent_experimental_confirmation=False,
        source_url="https://www.ncbi.nlm.nih.gov/protein/" + official["accession"])


def apply_adjudication(row, decisions):
    """Attach exact-sequence external review; retain blind prediction unchanged.

    This is a reviewed operational view, not a new classifier or held-out test.
    Novel homologues are not automatically assigned the reviewed sequence label.
    """
    result = dict(row)
    decision = decisions.get(row["sequence_sha256"])
    result["blind_family_status"] = row["family_status"]
    result["operational_family_status"] = row["family_status"]
    result["adjudication_disposition"] = "not_reviewed"
    result["adjudication_reason"] = ""
    if decision:
        result["operational_family_status"] = decision["operational_family_status"]
        result["adjudication_disposition"] = decision["disposition"]
        result["adjudication_reason"] = decision["reason"]
    return result
