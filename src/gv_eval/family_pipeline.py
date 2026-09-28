"""Reproducible member-A workflow (MAFFT, HMMER and CD-HIT required)."""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import json
import random
import subprocess

import yaml

from .classification import FamilyThresholds, classify_family, pfam_pass, apply_similarity_evidence
from .external import parse_domtbl, parse_hmm_metadata, require_tools
from .generation_data import parse_cdhit_clusters
from .io import (FastaRecord, STANDARD_AA, read_fasta, sha256_file,
                 write_fasta, write_json, write_tsv)
from .reference import audit_family_sources, prepare_reviewed_references
from .io import sequence_sha256


def table(path, rows):
    if not rows:
        raise ValueError(f"Cannot write an untyped empty audit table: {path}")
    write_tsv(path, rows, list(rows[0]))


def command(args, log, output=None):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as err:
        if output:
            with output.open("w", encoding="utf-8") as out:
                subprocess.run([str(x) for x in args], check=True, stdout=out, stderr=err)
        else:
            subprocess.run([str(x) for x in args], check=True, stdout=err, stderr=err)


def search(model, fasta, target):
    command(["hmmsearch", "--cpu", "1", "--seed", "42", "--max", "-T", "0",
             "--domT", "0", "--domtblout", target, model, fasta], target.with_suffix(".log"))
    return parse_domtbl(target)


def build_profiles(records_by_family, directory, minimum):
    directory.mkdir(parents=True, exist_ok=True)
    models, manifest = {}, []
    for family, records in sorted(records_by_family.items()):
        if len(records) < minimum:
            continue
        fasta = directory / (family + ".seeds.fasta")
        msa = directory / (family + ".alignment.fasta")
        model = directory / (family + ".hmm")
        write_fasta(sorted(records, key=lambda r: r.identifier), fasta)
        command(["mafft", "--thread", "1", "--auto", fasta], directory / (family + ".mafft.log"), msa)
        command(["hmmbuild", "--cpu", "1", "--seed", "42", "-n", family, model, msa],
                directory / (family + ".hmmbuild.log"))
        # HMMER embeds wall-clock DATE; strip it so equivalent reruns hash identically.
        text = model.read_text(encoding="utf-8")
        with model.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write("\n".join(line for line in text.splitlines() if not line.startswith("DATE ")) + "\n")
        models[family] = model
        manifest.append(dict(family=family, sequence_count=len(records),
                             training_ids=[r.identifier for r in sorted(records, key=lambda r:r.identifier)],
                             hmm_sha256=sha256_file(model), alignment_sha256=sha256_file(msa),
                             seeds_sha256=sha256_file(fasta)))
    return models, manifest


def classify_batch(records, models, pfam, directory, thresholds, similarity_references, annotations):
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "queries.fasta"
    valid = [r for r in records if r.sequence and not set(r.sequence) - STANDARD_AA]
    write_fasta(valid, fasta)
    pf_hits = search(pfam, fasta, directory / "pf00741.domtbl") if valid else {}
    family_hits = {f: search(p, fasta, directory / (f + ".domtbl")) if valid else {}
                   for f, p in models.items()}
    rows = [classify_family(r, pf_hits.get(r.identifier),
                           {f: hits[r.identifier] for f, hits in family_hits.items() if r.identifier in hits},
                           available_profiles=set(models), thresholds=thresholds) for r in records]
    reference_fasta = directory / "similarity_references.fasta"
    write_fasta(similarity_references, reference_fasta)
    blast_output = directory / "competitive_similarity.tsv"
    hits = defaultdict(list)
    if valid and similarity_references:
        database = directory / "blastdb" / "references"
        database.parent.mkdir(parents=True, exist_ok=True)
        command(["makeblastdb", "-in", reference_fasta, "-dbtype", "prot", "-out", database],
                directory / "makeblastdb.log")
        command(["blastp", "-query", fasta, "-db", database, "-seg", "yes",
                 "-evalue", "0.001", "-max_target_seqs", str(max(500, len(similarity_references))),
                 "-num_threads", "4", "-outfmt",
                 "6 qseqid sseqid nident qlen slen qstart qend sstart send evalue bitscore",
                 "-out", blast_output], directory / "blastp.log")
        for line in blast_output.read_text(encoding="utf-8").splitlines():
            q, target, identical, qlen, tlen, qs, qe, ts, te, evalue, score = line.split("\t")
            hits[q].append(dict(target_id=target, family=annotations[target],
                effective_identity=int(identical)/max(int(qlen), int(tlen)),
                query_coverage=(int(qe)-int(qs)+1)/int(qlen),
                target_coverage=(int(te)-int(ts)+1)/int(tlen)))
    return [apply_similarity_evidence(row, hits[row["sequence_id"]], thresholds) for row in rows]


def run_family_pipeline(root: Path, config_path: str, candidates_override: str | None = None):
    config_file = root / config_path
    config = yaml.safe_load(config_file.read_text(encoding="utf-8"))
    t = FamilyThresholds(**config["thresholds"])
    if config["folds"] < 2 or config["minimum_profile_sequences"] < 2:
        raise ValueError("At least two folds and two profile sequences are required")
    if not 0.6 <= config["cluster_identity"] <= 0.7 or not 0 < config["cluster_coverage"] <= 1:
        raise ValueError("CD-HIT word size 4 requires identity 0.6..0.7; coverage must be in (0,1]")
    tools = require_tools(["mafft", "hmmbuild", "hmmsearch", "cd-hit", "blastp", "makeblastdb"])
    refdir, modeldir, resultdir = [root / config[k] for k in ("reference_dir", "profile_dir", "result_dir")]
    for directory in (refdir, modeldir, resultdir):
        directory.mkdir(parents=True, exist_ok=True)
    reviewed, reviewed_manifest = prepare_reviewed_references(
        root / config["reviewed_snapshot"], root / config["reviewed_fasta"])
    table(refdir / "reviewed_source_manifest.tsv", reviewed_manifest)
    reviewed_ids = {"nat_" + sequence_sha256(r.sequence)[:20] for r in reviewed
                    if "fragment" not in r.description.lower() and "probable" not in r.description.lower()}
    # Work on relative source paths to keep manifests portable across OSes.
    sources = [root / p for p in config["sources"]]
    records, audit = audit_family_sources(sources)
    for row in audit:
        row["source_path"] = Path(row["source_path"]).relative_to(root).as_posix()
    by_id = {r.identifier: r for r in records}
    annotations = {row["reference_id"]: row["annotated_family"] for row in audit if row["seed_eligible"]}
    pfam = root / config["pfam_model"]
    pfmeta = parse_hmm_metadata(pfam)
    if pfmeta["accession"].split(".")[0] != "PF00741":
        raise ValueError("Expected PF00741 model")
    valid = [r for r in records if r.sequence and not set(r.sequence) - STANDARD_AA]
    natural_fasta = resultdir / "natural_unique.fasta"
    write_fasta(valid, natural_fasta)
    pf_hits = search(pfam, natural_fasta, resultdir / "natural_pf00741.domtbl")
    seeds = [r for r in records if r.identifier in annotations
             and (annotations[r.identifier] not in ("gvpa", "gvpj") or pfam_pass(pf_hits.get(r.identifier), t))]
    if not seeds:
        raise ValueError("No eligible seeds")
    seed_fasta = refdir / "provisional_seeds.fasta"
    write_fasta(sorted(seeds, key=lambda r:r.identifier), seed_fasta)
    cluster_fasta = refdir / "seed_cluster_representatives.fasta"
    command(["cd-hit", "-i", seed_fasta, "-o", cluster_fasta, "-c", str(config["cluster_identity"]),
             "-n", "4", "-aS", str(config["cluster_coverage"]), "-aL", str(config["cluster_coverage"]),
             "-g", "1", "-d", "0", "-T", "1", "-M", "2000"], resultdir / "clustering.log")
    clusters = parse_cdhit_clusters(str(cluster_fasta) + ".clstr")
    # Assign whole clusters to folds, stratified by annotation where possible.
    family_clusters = defaultdict(list)
    for members in clusters:
        family_clusters[tuple(sorted({annotations[i] for i in members}))].append(sorted(members))
    assignments, cluster_ids, representatives = {}, {}, set()
    for families, grouped in sorted(family_clusters.items()):
        grouped.sort(key=lambda members: members[0])
        for index, members in enumerate(grouped):
            fold = index % config["folds"]
            cluster_id = "cluster_" + members[0][4:]
            for identifier in members:
                assignments[identifier] = fold
                cluster_ids[identifier] = cluster_id
            # A multi-label cluster can contribute one seed per actual family,
            # but ALL its members are excluded together in held-out evaluation.
            for family in families:
                eligible = [i for i in members if annotations[i] == family
                            and (family not in ("gvpa", "gvpj") or i in reviewed_ids)]
                if eligible:
                    representatives.add(eligible[0])
    seed_manifest = [dict(reference_id=r.identifier, annotated_family=annotations[r.identifier],
                          cluster_id=cluster_ids[r.identifier], fold=assignments[r.identifier],
                          reviewed_anchor=r.identifier in reviewed_ids,
                          profile_representative=r.identifier in representatives) for r in seeds]
    table(refdir / "seed_manifest.tsv", seed_manifest)
    out_of_fold, control_rows, fold_manifests = {}, [], []
    for fold in range(config["folds"]):
        print(f"Building and checking fold {fold + 1}/{config['folds']}", flush=True)
        train = defaultdict(list)
        for identifier in sorted(representatives):
            if assignments[identifier] != fold:
                train[annotations[identifier]].append(by_id[identifier])
        directory = resultdir / f"fold_{fold}"
        models, metadata = build_profiles(train, directory / "profiles", config["minimum_profile_sequences"])
        heldout = [r for r in seeds if assignments[r.identifier] == fold]
        similarity_references = [r for r in seeds if assignments[r.identifier] != fold]
        rows = classify_batch(heldout, models, pfam, directory / "natural", t,
                              similarity_references, annotations)
        out_of_fold.update({row["sequence_id"]: row for row in rows})
        for row in rows:
            family = annotations[row["sequence_id"]]
            control_rows.append(dict(**row, control_type="natural_heldout", expected_family=family,
                                     fold=fold, parent_id=row["sequence_id"]))
        synthetic, meta = [], {}
        rng = random.Random(config["seed"] + fold)
        sampled = Counter()
        for record in heldout:
            family = annotations[record.identifier]
            if record.identifier not in representatives or sampled[family] >= 5:
                continue
            sampled[family] += 1
            shuffled = list(record.sequence)
            rng.shuffle(shuffled)
            for kind, sequence in (("shuffled", "".join(shuffled)),
                                   ("truncated", record.sequence[:max(1, len(record.sequence)//3)]),
                                   ("exact_duplicate", record.sequence)):
                identifier = f"{kind}_{record.identifier}"
                synthetic.append(FastaRecord(identifier, "", sequence))
                meta[identifier] = (kind, family, record.identifier)
        for row in classify_batch(synthetic, models, pfam, directory / "synthetic", t,
                                  similarity_references, annotations):
            kind, family, parent = meta[row["sequence_id"]]
            control_rows.append(dict(**row, control_type=kind, expected_family=family, fold=fold, parent_id=parent))
        fold_manifests.append(dict(fold=fold, profiles=metadata, heldout_ids=[r.identifier for r in heldout]))
    # Promote only records supported without their own homology cluster.
    accepted, statuses = defaultdict(list), {}
    for record in records:
        identifier = record.identifier
        evidence = out_of_fold.get(identifier)
        family = annotations.get(identifier)
        category = "possible_gvpa" if family == "gvpa" else "excluded_or_unresolved"
        if evidence and evidence["best_family"] == family:
            if evidence["family_status"] == "supported_gvpa" and family == "gvpa":
                category = "high_confidence_gvpa"
            elif evidence["family_status"] == "other_gvp":
                category = "gvpj_negative_control" if family == "gvpj" else "other_gvp_negative_control"
        statuses[identifier] = category
        if category in ("high_confidence_gvpa", "gvpj_negative_control", "other_gvp_negative_control"):
            accepted[family].append(record)
    for category in ("high_confidence_gvpa", "possible_gvpa", "gvpj_negative_control",
                     "other_gvp_negative_control", "excluded_or_unresolved"):
        write_fasta([r for r in records if statuses[r.identifier] == category], refdir / (category + ".fasta"))
    for row in audit:
        identifier = row["reference_id"]
        evidence = out_of_fold.get(identifier)
        row["reference_status"] = statuses[identifier] if row["seed_eligible"] else "excluded_or_unresolved"
        row["fold"] = assignments.get(identifier, "")
        row["cluster_id"] = cluster_ids.get(identifier, "")
        row["pf00741_pass"] = pfam_pass(pf_hits.get(identifier), t)
        row["pf00741_score"] = pf_hits[identifier].domain_score if identifier in pf_hits else None
        row["pf00741_coverage"] = pf_hits[identifier].model_coverage if identifier in pf_hits else 0.0
        row["family_status"] = evidence["family_status"] if evidence else "unsupported"
        row["best_family"] = evidence["best_family"] if evidence else ""
        row["family_margin"] = evidence["family_margin"] if evidence else None
        row["classification_reason"] = (row["annotation_reasons"] or
                                         (evidence["classification_reason"] if evidence else "not_eligible_for_profile_audit"))
    table(refdir / "reference_classification.tsv", audit)
    table(resultdir / "control_classification.tsv", control_rows)
    table(resultdir / "reference_out_of_fold_scores.tsv", list(out_of_fold.values()))
    # A/J profiles remain anchored to reviewed records; never bootstrap them
    # from the model's own predictions on the noisy nominal local references.
    final_train = defaultdict(list)
    for family, family_records in accepted.items():
        if family in ("gvpa", "gvpj"):
            continue
        seen = set()
        for record in sorted(family_records, key=lambda r:r.identifier):
            cluster = cluster_ids[record.identifier]
            if cluster not in seen:
                final_train[family].append(record)
                seen.add(cluster)
    for identifier in sorted(representatives):
        if annotations[identifier] in ("gvpa", "gvpj"):
            final_train[annotations[identifier]].append(by_id[identifier])
    models, profile_metadata = build_profiles(final_train, modeldir, config["minimum_profile_sequences"])
    if not {"gvpa", "gvpj"}.issubset(models):
        raise ValueError("Insufficient supported GvpA/GvpJ references to build final profiles")
    candidate_path = root / (candidates_override or config["candidates"])
    candidates = list(read_fasta(candidate_path))
    if not candidates or len({r.identifier for r in candidates}) != len(candidates):
        raise ValueError("Candidates must be nonempty with unique identifiers")
    candidate_rows = classify_batch(candidates, models, pfam, resultdir / "candidates", t, seeds, annotations)
    candidate_output = refdir.parent / "candidate_family_classification.tsv"
    table(candidate_output, candidate_rows)
    summary = dict(source_records=len(audit), unique_natural_sequences=len(records),
                   seed_sequences=len(seeds), seed_clusters=len(clusters),
                   reference_status_counts=dict(Counter(statuses.values())),
                   candidate_status_counts=dict(Counter(r["family_status"] for r in candidate_rows)),
                   available_profiles=sorted(models),
                   controls={}, pf00741_by_annotated_family={})
    for kind in sorted({r["control_type"] for r in control_rows}):
        subset = [r for r in control_rows if r["control_type"] == kind]
        summary["controls"][kind] = dict(total=len(subset), supported_gvpa=sum(r["family_status"] == "supported_gvpa" for r in subset),
            non_gvpa_called_gvpa=sum(r["family_status"] == "supported_gvpa" and r["expected_family"] != "gvpa" for r in subset))
    for family in sorted(set(annotations.values())):
        identifiers = [i for i, f in annotations.items() if f == family]
        summary["pf00741_by_annotated_family"][family] = dict(total=len(identifiers),
            pass_count=sum(pfam_pass(pf_hits.get(i), t) for i in identifiers))
    write_json(resultdir / "summary.json", summary)
    write_json(modeldir / "profile_manifest.json", dict(method="MAFFT + hmmbuild, cluster representative weighting",
        evidence_level="reviewed_A_J_anchors_plus_out_of_cluster_computational_support",
        sources={p: sha256_file(root / p) for p in config["sources"]},
        pfam=dict(**pfmeta, sha256=sha256_file(pfam)), profiles=profile_metadata,
        cross_validation=fold_manifests, config=config,
        limitation="Reviewed is not experimentally verified. Cluster-heldout audit is internal validation; no full-Pfam off-target scan."))
    outputs = [p for directory in (refdir, modeldir) for p in directory.iterdir()
               if p.is_file() and p.suffix in (".tsv", ".fasta", ".hmm", ".json")]
    outputs += [candidate_output, resultdir / "control_classification.tsv", resultdir / "summary.json",
                resultdir / "reference_out_of_fold_scores.tsv"]
    write_json(resultdir / "manifest.json", dict(config=config, tools=tools,
        inputs={p.relative_to(root).as_posix():sha256_file(p) for p in sources + [config_file, pfam, candidate_path,
                root / config["reviewed_snapshot"], root / config["pfam_metadata_snapshot"]]},
        outputs={p.relative_to(root).as_posix():sha256_file(p) for p in outputs},
        code={p.relative_to(root).as_posix():sha256_file(p) for p in [Path(__file__),
              root / "src/gv_eval/classification.py", root / "src/gv_eval/reference.py",
              root / "src/gv_eval/external.py"]}))
    print(json.dumps(summary, indent=2), flush=True)
    return summary
