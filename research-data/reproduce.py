"""Verify published transition data, weighted JS and MDS. Requires NumPy only.

Run from any directory: python path/to/research-data/reproduce.py
Reads the published files without modifying them or refitting a curve.
"""
from collections import defaultdict
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def rows(name):
    with (HERE / name).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def square(name):
    with (HERE / name).open(encoding='utf-8', newline='') as f:
        table = list(csv.reader(f))
    ids = table[0][1:]
    assert ids == [r[0] for r in table[1:]]
    return ids, np.array([r[1:] for r in table[1:]], dtype=float)


def weighted_js(x):
    totals = x.sum(axis=2)
    shares = totals / np.maximum(totals.sum(axis=1)[:, None], 1)
    p = x / np.maximum(totals[:, :, None], 1)
    result = np.zeros((len(x), len(x)))
    mass = np.zeros_like(result)
    for i in range(len(x)):
        mid = (p[i] + p) / 2
        # KL contributions with 0 log(0 / m) = 0.
        def kl(a):
            ratio = np.divide(a, mid, out=np.ones_like(mid), where=mid > 0)
            logs = np.zeros_like(mid)
            np.log2(ratio, out=logs, where=ratio > 0)
            return (a * logs).sum(axis=-1)
        js = (kl(p[i]) + kl(p)) / 2
        weight = (shares[i] + shares) / 2 * ((totals[i] >= 3) & (totals >= 3))
        result[i] = (js * weight).sum(axis=1)
        mass[i] = weight.sum(axis=1)
    np.fill_diagonal(result, 0)
    return result, mass


def main():
    manifest = json.loads((HERE / 'manifest.json').read_text(encoding='utf-8'))
    for name, record in manifest['csv_files'].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == record['sha256'], name
    profile_rows = rows('profiles.csv')
    profiles = {r['profile_id']: r for r in profile_rows}
    alphabet = [r['unit'] for r in rows('alphabet.csv')]
    index = {g: i for i, g in enumerate(alphabet)}
    counts = {}
    for mode in ('word', 'line'):
        tables = defaultdict(lambda: np.zeros((len(alphabet), len(alphabet))))
        observed = set()
        for r in rows(f'transitions-{mode}-counts.csv'):
            uid, g = r['profile_id'], r['from_unit']
            assert (uid, g) not in observed
            observed.add((uid, g))
            x = np.array([int(r['to_' + h]) for h in alphabet])
            assert x.sum() == int(r['row_total']) and np.all(x >= 0)
            tables[uid][index[g]] = x
        assert len(observed) == len(profiles) * len(alphabet)
        for r in rows(f'transitions-{mode}-probabilities.csv'):
            x = tables[r['profile_id']][index[r['from_unit']]]
            total = x.sum()
            assert int(r['row_total']) == total
            if total:
                np.testing.assert_allclose([float(r['to_' + h]) for h in alphabet], x / total, atol=1e-15, rtol=0)
            else:
                assert all(r['to_' + h] == '' for h in alphabet)
            np.testing.assert_allclose(float(r['row_weight']), total / max(tables[r['profile_id']].sum(), 1), atol=1e-15)
            assert int(r['row_eligible_ge3']) == int(total >= 3)
        for uid, r in profiles.items():
            members = r['constituent_units'].split(';')
            np.testing.assert_array_equal(tables[uid], sum((tables[u] for u in members), np.zeros_like(tables[uid])))
        counts[mode] = tables
    shared = [r['profile_id'] for r in sorted((r for r in profile_rows if r['shared_mds_index'] != ''), key=lambda r: int(r['shared_mds_index']))]
    distances = {}
    for mode in counts:
        d, mass = weighted_js(np.array([counts[mode][u] for u in shared]))
        ids, expected = square(f'pairwise/shared-{mode}-js.csv')
        assert ids == shared
        np.testing.assert_allclose(d, expected, atol=1e-12, rtol=0)
        _, expected_mass = square(f'pairwise/shared-{mode}-eligible-weight.csv')
        np.testing.assert_allclose(mass, expected_mass, atol=1e-12, rtol=0)
        distances[mode] = d
    distances['combined'] = (distances['word'] + distances['line']) / 2
    _, expected = square('pairwise/shared-combined-js.csv')
    np.testing.assert_allclose(distances['combined'], expected, atol=1e-12, rtol=0)
    for name in ('paragraph', 'all-nonlabel'):
        ids, expected = square(f'pairwise/matrix-{name}-word-js.csv')
        actual, _ = weighted_js(np.array([counts['word'][u] for u in ids]))
        np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)
    for model in rows('mds/models.csv'):
        name = model['model']
        table = rows('mds/' + name + '.csv')
        ids = [r['profile_id'] for r in table]
        cols = [k for k in table[0] if k.startswith('MDS')]
        saved = np.array([[float(r[k]) for k in cols] for r in table])
        ix = [shared.index(u) for u in ids]
        mode = name.split('-')[-1]
        d = distances[mode][np.ix_(ix, ix)]
        gram = -.5 * (d - d.mean(axis=0) - d.mean(axis=1)[:, None] + d.mean())
        values, vectors = np.linalg.eigh(gram)
        positive = values > 1e-12
        expected_gram = (vectors[:, positive] * values[positive]) @ vectors[:, positive].T
        # Gram comparison tolerates eigenvector signs and rotations at ties.
        np.testing.assert_allclose(saved @ saved.T, expected_gram, atol=1e-10, rtol=0)
        assert saved.shape[1] == int(model['positive_dimensions']) == positive.sum()
        print(name, len(ids), 'profiles,', positive.sum(), 'positive axes verified')
    print(len(manifest['csv_files']), 'CSV hashes verified. Counts, probabilities, pools, JS and MDS agree.')


if __name__ == '__main__':
    main()
