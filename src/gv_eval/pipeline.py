from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .analysis import (
    ablation_results,
    build_strategy_results,
    correlation_tables,
    random_baseline,
    save_figures,
    threshold_robustness,
    weight_robustness,
)
from .external import (
    download_pfam_hmm,
    parse_blast,
    parse_domtbl,
    parse_hmm_metadata,
    require_tools,
    run_command,
    tool_version,
)
from .io import (
    FastaRecord,
    STANDARD_AA,
    read_alignment,
    read_fasta,
    sha256_file,
    write_fasta,
    write_json,
    write_tsv,
)
from .metrics import (
    candidate_distance_matrix,
    closest_reference_metrics,
    conservation_scores,
    domain_metrics,
    find_conserved_sites,
)
from .reference import build_clean_reference
from .selection import pareto_ranking, weighted_ranking


BLAST_COLUMNS = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"


def _path(root: Path, value: str) -> Path:
    return root / value


def _write_frame(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format="%.10g")


def _selection_fasta(selected: pd.DataFrame, candidates: dict[str, FastaRecord], path: Path) -> None:
    write_fasta([candidates[value] for value in selected["sequence_id"]], path)


def run_pipeline(root: str | Path, config_path: str | Path) -> dict:
    root = Path(root).resolve()
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = root / config_path
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    os.environ.setdefault("PYTHONHASHSEED", str(config["seed"]))

    tools = require_tools(["hmmsearch", "mafft", "blastp", "makeblastdb", "cd-hit"])
    processed = _path(root, config["outputs"]["processed_dir"])
    result = _path(root, config["outputs"]["result_dir"])
    raw = result / "raw"
    tables = result / "tables"
    figures = result / "figures"
    selections_dir = result / "selections"
    work = result / "work"
    for directory in (processed, raw, tables, figures, selections_dir, work):
        directory.mkdir(parents=True, exist_ok=True)

    candidate_path = _path(root, config["inputs"]["candidates"])
    candidates_list = list(read_fasta(candidate_path))
    candidate_ids = [row.identifier for row in candidates_list]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("Candidate FASTA contains duplicate identifiers")
    if not candidates_list or any(not row.sequence or set(row.sequence) - STANDARD_AA for row in candidates_list):
        raise ValueError("Candidate FASTA contains empty or non-standard protein sequences")
    candidates = {row.identifier: row for row in candidates_list}

    reference_paths = [_path(root, value) for value in config["inputs"]["nominal_gvpa_sources"]]
    reference, reference_mapping, reference_summary = build_clean_reference(
        reference_paths,
        required_label=config["reference"]["required_label"],
        excluded_labels=config["reference"]["excluded_labels"],
        excluded_header_terms=config["reference"]["excluded_header_terms"],
        minimum_length=config["reference"]["min_length"],
        maximum_length=config["reference"]["max_length"],
    )
    reference_clean = processed / "reference_clean.fasta"
    write_fasta(reference, reference_clean)
    write_tsv(processed / "reference_mapping.tsv", reference_mapping, list(reference_mapping[0]))
    write_json(processed / "reference_cleaning_summary.json", reference_summary)

    representative_path = processed / "reference_cdhit90.fasta"
    run_command(
        [
            tools["cd-hit"], "-i", str(reference_clean), "-o", str(representative_path),
            "-c", str(config["reference"]["cluster_identity"]), "-n",
            str(config["reference"]["cluster_word_length"]), "-d", "0", "-M", "0",
            "-T", str(config["threads"]),
        ],
        stdout_path=work / "cdhit.log",
        stderr_path=work / "cdhit.stderr.log",
    )
    representatives = list(read_fasta(representative_path))
    if not representatives:
        raise RuntimeError("CD-HIT produced no reference representative")

    reference_alignment = processed / "reference_cdhit90_alignment.fasta"
    run_command(
        [
            tools["mafft"], "--auto", "--thread",
            str(config["conservation"]["alignment_threads"]), str(representative_path),
        ],
        stdout_path=reference_alignment,
        stderr_path=work / "mafft_reference.log",
    )
    combined_alignment = processed / "reference_candidate_alignment.fasta"
    run_command(
        [
            tools["mafft"], "--add", str(candidate_path), "--keeplength", "--thread",
            str(config["conservation"]["alignment_threads"]), str(reference_alignment),
        ],
        stdout_path=combined_alignment,
        stderr_path=work / "mafft_add_candidates.log",
    )
    reference_msa = read_alignment(reference_alignment)
    combined_msa = read_alignment(combined_alignment)
    sites = find_conserved_sites(
        reference_msa,
        config["conservation"]["maximum_reference_gap_fraction"],
        config["conservation"]["minimum_consensus_fraction"],
    )
    site_rows = [site.__dict__ for site in sites]
    write_tsv(processed / "conserved_sites.tsv", site_rows, list(site_rows[0]))
    conservation = conservation_scores(combined_msa, sites, candidate_ids)

    hmm_path = _path(root, config["pfam"]["model_path"])
    download_pfam_hmm(config["pfam"]["url"], config["pfam"]["accession"], hmm_path)
    hmm_metadata = parse_hmm_metadata(hmm_path)
    candidate_domtbl = raw / "candidate_pfam.domtbl"
    reference_domtbl = raw / "reference_pfam.domtbl"
    hmm_common = [
        "--noali", "-E", str(config["pfam"]["permissive_sequence_evalue"]),
        "--domE", str(config["pfam"]["permissive_domain_evalue"]),
    ]
    run_command(
        [tools["hmmsearch"], *hmm_common, "--domtblout", str(candidate_domtbl), str(hmm_path), str(candidate_path)],
        stdout_path=raw / "candidate_pfam.txt", stderr_path=work / "hmm_candidate.log",
    )
    run_command(
        [tools["hmmsearch"], *hmm_common, "--domtblout", str(reference_domtbl), str(hmm_path), str(representative_path)],
        stdout_path=raw / "reference_pfam.txt", stderr_path=work / "hmm_reference.log",
    )
    domain, calibration = domain_metrics(
        candidate_ids,
        parse_domtbl(candidate_domtbl),
        parse_domtbl(reference_domtbl),
        hmm_metadata.get("ga_domain"),
        config["pfam"]["minimum_coverage_floor"],
        config["pfam"]["minimum_coverage_ceiling"],
    )
    write_json(processed / "domain_calibration.json", calibration)

    blastdb_dir = work / "blastdb"
    blastdb_dir.mkdir(parents=True, exist_ok=True)
    reference_db = blastdb_dir / "reference"
    candidate_db = blastdb_dir / "candidates"
    run_command(
        [tools["makeblastdb"], "-in", str(reference_clean), "-dbtype", "prot", "-out", str(reference_db)],
        stdout_path=work / "makeblastdb_reference.log", stderr_path=work / "makeblastdb_reference.stderr.log",
    )
    run_command(
        [tools["makeblastdb"], "-in", str(candidate_path), "-dbtype", "prot", "-out", str(candidate_db)],
        stdout_path=work / "makeblastdb_candidates.log", stderr_path=work / "makeblastdb_candidates.stderr.log",
    )
    reference_blast = raw / "candidate_vs_reference.tsv"
    candidate_blast = raw / "candidate_vs_candidate.tsv"
    run_command(
        [
            tools["blastp"], "-query", str(candidate_path), "-db", str(reference_db), "-out", str(reference_blast),
            "-outfmt", BLAST_COLUMNS, "-evalue", str(config["blast"]["evalue"]),
            "-max_target_seqs", str(config["blast"]["maximum_reference_targets"]),
            "-num_threads", str(config["threads"]),
        ]
    )
    run_command(
        [
            tools["blastp"], "-query", str(candidate_path), "-db", str(candidate_db), "-out", str(candidate_blast),
            "-outfmt", BLAST_COLUMNS, "-evalue", str(config["blast"]["evalue"]),
            "-max_target_seqs", str(config["blast"]["maximum_candidate_targets"]),
            "-num_threads", str(config["threads"]),
        ]
    )
    novelty = closest_reference_metrics(candidate_ids, parse_blast(reference_blast))
    distance, uniqueness = candidate_distance_matrix(candidate_ids, parse_blast(candidate_blast))
    np.save(processed / "candidate_distance.npy", distance, allow_pickle=False)
    write_json(processed / "candidate_distance_ids.json", candidate_ids)

    rows: list[dict] = []
    for record in candidates_list:
        identifier = record.identifier
        rows.append(
            {
                "sequence_id": identifier,
                "length": len(record.sequence),
                **domain[identifier],
                **conservation[identifier],
                **novelty[identifier],
                **uniqueness[identifier],
            }
        )
    frame = pd.DataFrame(rows)
    metrics = list(config["selection"]["metrics"])
    _write_frame(frame, tables / "candidate_scores.tsv")

    pearson, spearman = correlation_tables(frame, metrics)
    pearson.to_csv(tables / "correlation_pearson.tsv", sep="\t", lineterminator="\n", float_format="%.10g")
    spearman.to_csv(tables / "correlation_spearman.tsv", sep="\t", lineterminator="\n", float_format="%.10g")
    selections, strategy_rows, overlap_rows = build_strategy_results(
        frame, metrics, config["selection"]["equal_weights"], config["selection"]["budgets"], distance
    )
    for (strategy, budget), selected in selections.items():
        _write_frame(selected, selections_dir / f"{strategy}_top{budget}.tsv")
        _selection_fasta(selected, candidates, selections_dir / f"{strategy}_top{budget}.fasta")
    strategy_frame = pd.DataFrame(strategy_rows)
    overlap_frame = pd.DataFrame(overlap_rows)
    _write_frame(strategy_frame, tables / "strategy_summary.tsv")
    _write_frame(overlap_frame, tables / "strategy_overlap.tsv")

    primary_k = config["selection"]["primary_k"]
    ablation_rows, ablation_selections = ablation_results(
        frame, metrics, config["selection"]["equal_weights"], primary_k
    )
    _write_frame(pd.DataFrame(ablation_rows), tables / "ablation_summary.tsv")
    write_json(tables / "ablation_selections.json", ablation_selections)

    weight_rows, frequency = weight_robustness(
        frame, metrics, config["selection"]["equal_weights"], primary_k,
        config["robustness"]["weight_samples"], config["robustness"]["relative_weight_sigma"], config["seed"],
    )
    weight_frame = pd.DataFrame(weight_rows)
    _write_frame(weight_frame, tables / "weight_robustness.tsv")
    _write_frame(
        pd.DataFrame([{"sequence_id": key, "selection_frequency": value} for key, value in sorted(frequency.items())]),
        tables / "weight_selection_frequency.tsv",
    )
    weighted = weighted_ranking(frame, metrics, config["selection"]["equal_weights"])
    threshold_rows = threshold_robustness(
        weighted, primary_k, calibration["domain_score_cutoff"], calibration["model_coverage_cutoff"],
        config["robustness"]["domain_cutoff_multipliers"], config["robustness"]["coverage_offsets"],
    )
    threshold_frame = pd.DataFrame(threshold_rows)
    _write_frame(threshold_frame, tables / "threshold_robustness.tsv")

    random_rows, random_comparisons = random_baseline(
        frame, strategy_rows, metrics, distance, primary_k,
        config["robustness"]["random_baseline_samples"], config["seed"],
    )
    _write_frame(pd.DataFrame(random_rows), tables / "random_baseline_samples.tsv")
    _write_frame(pd.DataFrame(random_comparisons), tables / "random_baseline_comparison.tsv")

    pareto = pareto_ranking(frame, metrics)
    _write_frame(pareto, tables / "pareto_ranking.tsv")
    save_figures(frame, metrics, pearson, spearman, strategy_frame, overlap_frame, weight_frame, threshold_frame, figures)

    summary = {
        "candidate_count": len(frame),
        "clean_reference_count": len(reference),
        "reference_representative_count": len(representatives),
        "conserved_site_count": len(sites),
        "candidate_domain_hit_count": int((frame["domain_status"] == "detected").sum()),
        "candidate_domain_pass_count": int(frame["domain_pass"].sum()),
        "domain_calibration": calibration,
        "metric_summary": {
            metric: {
                "mean": float(frame[metric].mean()), "median": float(frame[metric].median()),
                "min": float(frame[metric].min()), "max": float(frame[metric].max()),
            }
            for metric in metrics
        },
        "correlation_pearson": pearson.to_dict(),
        "correlation_spearman": spearman.to_dict(),
        "pareto_front_size": int((pareto["pareto_rank"] == 0).sum()),
        "strategy_summary": strategy_rows,
        "strategy_overlap": overlap_rows,
        "ablation_summary": ablation_rows,
        "weight_robustness": {
            "samples": len(weight_frame),
            "jaccard_mean": float(weight_frame["top_k_jaccard"].mean()),
            "jaccard_min": float(weight_frame["top_k_jaccard"].min()),
            "rank_spearman_mean": float(weight_frame["rank_spearman"].mean()),
            "rank_spearman_min": float(weight_frame["rank_spearman"].min()),
        },
        "threshold_robustness": threshold_rows,
        "random_baseline_comparison": random_comparisons,
    }
    write_json(result / "summary.json", summary)

    input_files = [config_path, candidate_path, *reference_paths]
    output_files = [
        processed / "reference_clean.fasta", processed / "reference_cdhit90.fasta",
        processed / "conserved_sites.tsv", processed / "candidate_distance.npy",
        tables / "candidate_scores.tsv", tables / "strategy_summary.tsv",
        tables / "ablation_summary.tsv", tables / "weight_robustness.tsv",
        result / "summary.json",
    ]
    manifest = {
        "pipeline_version": "1.0.0",
        "seed": config["seed"],
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "tools": {name: {"path": path, "version": tool_version(path)} for name, path in tools.items()},
        "pfam": {
            **hmm_metadata,
            "url": config["pfam"]["url"],
            "sha256": sha256_file(hmm_path),
        },
        "inputs": {path.relative_to(root).as_posix(): sha256_file(path) for path in input_files},
        "outputs": {path.relative_to(root).as_posix(): sha256_file(path) for path in output_files},
        "config": config,
    }
    write_json(result / "manifest.json", manifest)
    return summary
