"""Create the M1 input audit without filtering or altering any source sequence."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re
from statistics import mean, median

from audit_fasta import audit, read_fasta, STANDARD_AA


def percentile(values, fraction):
    ordered = sorted(values)
    location = (len(ordered) - 1) * fraction
    low = int(location)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (location - low)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    candidates_path = root / 'data/candidates/t05/generated_gvp.fasta'
    source_paths = [root / 'data/raw/design/GvpA.fasta'] + sorted((root / 'data/raw/t05/gvpa').glob('*.fasta'))
    assert len(source_paths) == 4, 'Expected four preserved nominal GvpA source files'
    source_sets = {p.relative_to(root).as_posix(): {s for _, s in read_fasta(p)} for p in source_paths}
    reference_union = set().union(*source_sets.values())
    mixed = {s for _, s in read_fasta(root / 'data/raw/t05/real_gvp.fasta')}
    candidates = list(read_fasta(candidates_path))
    if not candidates:
        raise ValueError('Candidate file is empty')
    lengths = [len(s) for _, s in candidates]
    candidate_sequences = {s for _, s in candidates}
    unique_contributions = {
        name: len(sequences - set().union(*(other for other_name, other in source_sets.items() if other_name != name)))
        for name, sequences in source_sets.items()
    }
    unique_header_labels = {}
    for path in source_paths:
        name = path.relative_to(root).as_posix()
        others = set().union(*(values for other, values in source_sets.items() if other != name))
        labels = Counter()
        for header, sequence in read_fasta(path):
            if sequence not in others:
                matches = set(re.findall(r'\bgvp([acfgijklmnsuvwyz])\b', header, flags=re.I))
                labels.update({'Gvp' + suffix.upper() for suffix in matches} or {'unlabelled'})
        unique_header_labels[name] = dict(sorted(labels.items()))
    summary = {
        'stage': 'M1 input audit; no domain scan, alignment, family assignment or functional validation performed',
        'candidate_file': candidates_path.relative_to(root).as_posix(),
        'candidate_records': len(candidates),
        'candidate_unique_sequences': len(candidate_sequences),
        'candidate_min_length': min(lengths),
        'candidate_max_length': max(lengths),
        'candidate_mean_length': mean(lengths),
        'candidate_median_length': median(lengths),
        'candidate_q1_length': percentile(lengths, 0.25),
        'candidate_q3_length': percentile(lengths, 0.75),
        'candidate_standard_only_records': sum(set(s) <= STANDARD_AA for _, s in candidates),
        'candidate_length_512_records': sum(len(s) == 512 for _, s in candidates),
        'candidate_exact_matches_mixed_reference': len(candidate_sequences & mixed),
        'candidate_exact_matches_nominal_gvpa_union': len(candidate_sequences & reference_union),
        'candidate_length_counts': dict(sorted(Counter(lengths).items())),
        'nominal_gvpa_union_unique_sequences': len(reference_union),
        'nominal_gvpa_union_length_counts': dict(sorted(Counter(map(len, reference_union)).items())),
        'nominal_gvpa_source_unique_contributions': unique_contributions,
        'source_exclusive_sequence_record_header_labels': unique_header_labels,
        'files': [audit(p, root) for p in sorted((root / 'data').rglob('*.fasta'))],
    }
    out = root / 'results/input_audit'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    with (out / 'candidate_records.tsv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['sequence_id', 'length', 'standard_aa_only', 'at_generation_length_cap', 'exact_match_mixed_reference', 'exact_match_nominal_gvpa_union', 'family_status'])
        for header, sequence in candidates:
            writer.writerow([header.split()[0], len(sequence), set(sequence) <= STANDARD_AA, len(sequence) == 512, sequence in mixed, sequence in reference_union, 'not_assessed'])
    print(json.dumps({key: value for key, value in summary.items() if key not in ['files', 'candidate_length_counts', 'nominal_gvpa_union_length_counts']}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
