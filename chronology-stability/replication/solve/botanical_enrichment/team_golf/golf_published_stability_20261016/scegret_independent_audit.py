"""Independently validate scEGRET's original R execution; not a replacement run."""
from pathlib import Path
import argparse
import csv
import hashlib
import itertools
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_matrix(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        rows = list(csv.reader(handle))
    return rows[0][1:], [row[0] for row in rows[1:]], np.array([[float(v) for v in row[1:]] for row in rows[1:]])


def bins(positions, steps):
    """R cut.default numeric-break semantics: right closed, expanded endpoints."""
    p = np.asarray(positions, dtype=float)
    lo, hi = float(p.min()), float(p.max())
    width = hi-lo
    if width == 0:
        change = abs(lo) if lo else 1.
        edges = np.linspace(lo-change/1000, hi+change/1000, steps+1)
    else:
        edges = np.linspace(lo, hi, steps+1)
        edges[0] -= width/1000
        edges[-1] += width/1000
    assigned = np.searchsorted(edges, p, side='left')
    assert np.all((assigned >= 1) & (assigned <= steps))
    occupancy = np.bincount(assigned, minlength=steps+1)[1:]
    # droplevels changes internal integer factor codes, not numeric labels.
    occupied = np.flatnonzero(occupancy)+1
    factor_codes = np.searchsorted(occupied, assigned)+1
    return assigned, factor_codes, occupancy


def expected(matrix, features, observations, positions, steps):
    assigned, factor_codes, occupancy = bins(positions, steps)
    keep_bins = np.flatnonzero(occupancy > 3)+1
    out = dict(assignment=assigned, factor_codes=factor_codes, occupancy=occupancy,
               retained_bins=keep_bins, sufficient=len(keep_bins) >= 2)
    if not out['sufficient']:
        return out
    detected = (matrix != 0).mean(axis=1) > .1
    selected = matrix[detected]
    feature_names = np.asarray(features)[detected]
    denominator = np.log1p(selected.mean(axis=1))
    average = np.column_stack([np.log1p(selected[:, assigned == b].mean(axis=1)) for b in keep_bins])
    contributions = np.abs(np.diff(average, axis=1))/denominator[:, None]/len(selected)
    pairs = list(itertools.combinations(range(len(keep_bins)), 2))
    full = [float(np.mean(abs(average[:, a]-average[:, b])/denominator)) for a, b in pairs]
    out.update(features=feature_names.tolist(), feature_count=len(selected),
        adjacent=contributions.sum(axis=0), contributions=contributions,
        cellnumber=np.diff(occupancy[keep_bins-1]),
        full_pairs=[[int(keep_bins[a]), int(keep_bins[b])] for a, b in pairs],
        full_enthalpy=full,
        bin_observations={str(b):[observations[i] for i in np.flatnonzero(assigned == b)] for b in keep_bins})
    return out


def as_list(value):
    return value if isinstance(value, list) else [value]


def numeric_error(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    assert a.shape == b.shape, (a.shape, b.shape)
    return float(np.max(np.abs(a-b))) if a.size else 0.


def main(partial=None):
    inputs = json.loads((HERE/'scegret_inputs.json').read_text())
    manifest = json.loads((HERE/'scegret_source/source_manifest.json').read_text())
    for item in manifest['files']:
        assert sha(HERE/'scegret_source'/item['path']) == item['sha256']
    for path, digest in inputs['source_sha256'].items():
        assert sha(ROOT/path) == digest, path
    paths = [Path(__file__), HERE/'scegret_inputs.json', HERE/'run_scegret.mjs',
             HERE/'scegret_source/source_manifest.json', HERE/'scegret_source/R/egret_flap.R',
             HERE/'scegret_source/R/egret_fly.R', HERE/'normalize_scegret_output.py']
    expected_runs = {}
    matrix_checks = []
    for cohort in inputs['cohorts']:
        matrices = {}
        for matrix_id, filename in cohort['matrices'].items():
            path = HERE/filename; paths.append(path)
            observations, features, matrix = read_matrix(path)
            assert observations == [m['id'] for m in cohort['meta']]
            assert len(set(features)) == len(features) and np.all(matrix >= 0)
            assert np.isfinite(matrix).all() and (matrix.sum(axis=0) > 0).all()
            matrices[matrix_id] = (features, matrix)
            if matrix_id.endswith('per_million'):
                assert np.max(abs(matrix.sum(axis=0)-1e6)) < 1e-5
            else:
                assert np.array_equal(matrix, np.round(matrix))
            for candidate in cohort['candidates']:
                for steps in inputs['design']['steps'][cohort['cohort']]:
                    key = (cohort['cohort'], matrix_id, candidate['id'], steps)
                    expected_runs[key] = expected(matrix, features, observations, candidate['positions'], steps)
        for mode in ['word', 'line']:
            raw_features, raw = matrices[mode+'_raw']
            normalized_features, normalized = matrices[mode+'_per_million']
            assert raw_features == normalized_features
            error = numeric_error(raw/raw.sum(axis=0)*1e6, normalized)
            assert error < 1e-5
            matrix_checks.append(dict(cohort=cohort['cohort'], channel=mode,
                feature_count=len(raw_features), observations=raw.shape[1], per_million_max_error=error))
    out = dict(status='prepared', method='Independent Python validation of unmodified scEGRET R output',
        source_commit=manifest['commit'], matrix_checks=matrix_checks,
        expected_run_count=len(expected_runs),
        runs_with_factor_code_label_difference=sum(not np.array_equal(v['assignment'], v['factor_codes']) for v in expected_runs.values()),
        source_sha256={str(p.relative_to(HERE)):sha(p) for p in paths}, runs=[],
        exclusions='Does not audit optional Wilcoxon p-values or biological interpretation; no new curve fits or replacement scEGRET run.')
    result_path = HERE/(partial or 'scegret_results.json')
    if result_path.exists():
        snapshot = result_path.read_bytes()
        result = json.loads(snapshot)
        if result['status'] == 'complete' or partial:
            if 'published_source_manifest' in result:
                assert result['published_source_manifest'] == manifest
            assert result['runtime']['webr'] == '0.6.0'
            out['source_sha256'][result_path.name] = hashlib.sha256(snapshot).hexdigest()
            seen = set()
            for cohort in result['cohorts']:
                for run in cohort['runs']:
                    key = (cohort['cohort'], run['matrix'], run['candidate'], run['steps'])
                    assert key not in seen; seen.add(key)
                    e = expected_runs[key]
                    assert numeric_error(run['occupancy'], e['occupancy']) == 0
                    if not e['sufficient']:
                        assert run['status'] == 'insufficient_occupied_bins'
                        out['runs'].append(dict(key=list(key), status='insufficient_occupied_bins_verified'))
                        continue
                    assert run['status'] == 'complete', (key, run['status'])
                    assert run['feature_count'] == e['feature_count'] and run['features'] == e['features']
                    assert numeric_error(run['retained_bins'], e['retained_bins']) == 0
                    assert numeric_error(run['assignment'], e['assignment']) == 0, ('Original-bin assignment mismatch', key)
                    if 'factor_codes' in run:
                        assert numeric_error(run['factor_codes'], e['factor_codes']) == 0
                    assert {k:as_list(v) for k, v in run['bin_observations'].items()} == e['bin_observations']
                    adjacent = as_list(run['adjacent'])
                    adjacent_error = numeric_error([r['Relative_Step_Change'] for r in adjacent], e['adjacent'])
                    step_error = numeric_error([r['StepEnthalpy'] for r in adjacent], e['adjacent']*e['feature_count'])
                    assert numeric_error([r['CellNumber'] for r in adjacent], e['cellnumber']) == 0
                    assert numeric_error([r['Pseudostep'] for r in adjacent], np.arange(1, len(adjacent)+1)) == 0
                    assert numeric_error(run['adjacent_bin_pairs'], np.column_stack([e['retained_bins'][:-1], e['retained_bins'][1:]])) == 0
                    assert numeric_error(run['full_pairs'], e['full_pairs']) == 0
                    full_error = numeric_error(as_list(run['full_enthalpy']), e['full_enthalpy'])
                    top_error = 0.
                    for j, top in enumerate(run['top_adjacent_features']):
                        indices = [e['features'].index(name) for name in as_list(top['features'])]
                        actual = np.asarray(as_list(top['contribution']))
                        top_error = max(top_error, numeric_error(actual, e['contributions'][indices, j]))
                        cutoff = np.sort(e['contributions'][:, j])[-min(10, e['feature_count'])]
                        assert np.all(e['contributions'][indices, j] >= cutoff-1e-12)
                    assert max(adjacent_error, full_error, top_error) < 1e-10
                    assert step_error < 1e-8
                    out['runs'].append(dict(key=list(key), status='passed', retained_bins=e['retained_bins'].tolist(),
                        feature_count=e['feature_count'], adjacent_error=adjacent_error, full_pair_error=full_error,
                        top_feature_error=top_error, unnormalized_sum_error=step_error,
                        factor_codes_differ=not np.array_equal(e['assignment'], e['factor_codes'])))
            if partial:
                assert seen <= set(expected_runs)
                out['status'] = 'partial_passed'
            else:
                assert seen == set(expected_runs) and result['run_count'] == len(expected_runs)
                assert sha(result_path) == out['source_sha256'][result_path.name]
                out['status'] = 'passed'
    destination = 'scegret_partial_independent_audit.json' if partial else 'scegret_independent_audit.json'
    (HERE/destination).write_text(json.dumps(out, separators=(',', ':'), allow_nan=False)+'\n')
    print(json.dumps(dict(status=out['status'], expected=len(expected_runs), validated=len(out['runs']),
        factor_code_label_difference_runs=out['runs_with_factor_code_label_difference'])))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--partial', help='Audit currently saved cases from a named child output; does not imply all648 completed.')
    main(parser.parse_args().partial)
