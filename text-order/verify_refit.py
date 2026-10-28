"""Compare a fresh refit with the frozen release. Verification only, no fitting."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import runtime as R


def verify(output):
    audit = json.loads((output / 'refit-audit.json').read_text(encoding='utf-8'))
    if audit['status'] not in ('complete', 'smoke-test'):
        raise ValueError('Refit is unfinished')
    fresh = output / 'chronology-stability/replication/solve/botanical_enrichment/team_golf'
    checks = []
    total = 0
    for rep, files in audit['result_sha256'].items():
        folder = 'golf_published_stability_20261016' if rep == 'ZL' else 'golf_rf1b_stability_20261016'
        for name, digest in files.items():
            path = fresh / folder / name
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError('Changed refit file: ' + str(path))
            expected = json.loads((R.TEAM / folder / name).read_text(encoding='utf-8'))
            actual = json.loads(path.read_text(encoding='utf-8'))
            old = {(r['replicate'], r['start_index']): r for r in expected['runs']}
            seen = set()
            maximum = {k: 0. for k in ('positions', 'objective', 'node_mse', 'edge_energy', 'star_energy')}
            for got in actual['runs']:
                key = got['replicate'], got['start_index']
                if key in seen:
                    raise ValueError('Duplicate fit: ' + str(key))
                seen.add(key)
                saved = old[key]
                if got['converged'] != saved['converged']:
                    raise ValueError('Convergence differs: ' + name + str(key))
                for k in maximum:
                    error = float(np.max(abs(np.asarray(got[k]) - np.asarray(saved[k]))))
                    maximum[k] = max(maximum[k], error)
                    np.testing.assert_allclose(got[k], saved[k], atol=1e-10, rtol=0,
                        err_msg=name + ' ' + str(key) + ' ' + k)
            if audit['full_published_design'] and seen != set(old):
                raise ValueError('Missing fits in ' + name)
            total += len(seen)
            checks.append(dict(representation=rep, setting=name, fits=len(seen), maximum_errors=maximum))
    if len(checks) != 8 or total != audit['completed_fits'] or total != audit['planned_fits']:
        raise ValueError('Incomplete fit accounting')
    score_checks = {}
    if audit['full_published_design']:
        if total != 6400:
            raise ValueError('Full refit must contain 6400 fits')
        for unit in ('pages', 'bifolia'):
            with (R.SITE / 'path-v2' / (unit + '.csv')).open(encoding='utf-8', newline='') as f:
                expected = list(csv.DictReader(f))
            with (output / 'path-v2' / (unit + '.csv')).open(encoding='utf-8', newline='') as f:
                actual = list(csv.DictReader(f))
            if [r['id'] for r in actual] != [r['id'] for r in expected]:
                raise ValueError('Different cohort: ' + unit)
            fields = [k for k in expected[0] if k.endswith(('_score', '_rank', '_retained_fits'))]
            maximum = 0.
            for a, b in zip(actual, expected):
                for key in fields:
                    error = abs(float(a[key]) - float(b[key]))
                    maximum = max(maximum, error)
                    np.testing.assert_allclose(float(a[key]), float(b[key]), atol=1e-10, rtol=0)
            score_checks[unit] = maximum
    return dict(status='passed', full_published_design=audit['full_published_design'],
        fits_verified=total, settings=checks, score_maximum_errors=score_checks,
        note='A smoke test checks only the requested fits, not all 6400 or the consensus.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(verify(args.output), indent=2))
