"""Frozen cohorts and fold-specific classical MDS for global topology tests.

Retains no model arrays. JS is used as the squared-distance input, as in Golf.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import numpy as np

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent
ROOT = TEAM.parents[2]
DATA = ROOT / 'apps/clustermatrix/research-data'
SEED = 20261015


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def rows(path):
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def load_data():
    previous = module('progressive_source_for_global_topology', TEAM / 'progressive_mds_20261015/common.py')
    source = previous.load_data()
    ix = np.array([i for i, m in enumerate(source['meta']) if m['kind'] == 'para'])
    meta = [dict(source['meta'][i]) for i in ix]
    blocks = np.array([m['bifolium'] for m in meta])
    assert len(meta) == 207 and len(set(blocks)) == 50
    pools = [p for p in rows(DATA / 'profiles.csv') if p['profile_id'].startswith('pool|') and p['shared_mds_id']]
    pmeta = [dict(id=p['shared_mds_id'], profile_id=p['profile_id'], group=p['golf21_code'],
                  kind='pool', bifolium=p['bifolium'], words=int(p['readable_words']),
                  section=p['section'], hand=p['lfd_hand']) for p in pools]
    pds = {}
    for mode in ['word', 'line', 'combined']:
        table = {r['profile_id']: r for r in rows(DATA / f'pairwise/shared-{mode}-js.csv')}
        pds[mode] = np.array([[float(table[a['profile_id']][b['profile_id']]) for b in pmeta] for a in pmeta])
    assert len(pmeta) == 50
    return dict(meta=meta, blocks=blocks, ds={k:d[np.ix_(ix, ix)] for k,d in source['ds'].items()},
                global_t=source['global_t'][ix], reference=source['reference'][ix],
                full_meta=source['meta'], full_ds=source['ds'],
                pool_meta=pmeta, pool_ds=pds,
                sources=source['sources'] + [TEAM / 'progressive_mds_20261015/common.py', HERE / 'common.py'])


def embed(d, reference=None):
    g = -.5 * (d - d.mean(0) - d.mean(1)[:, None] + d.mean())
    val, vec = np.linalg.eigh(g)
    ix = np.argsort(val)[::-1]
    val, vec = val[ix], vec[:, ix]
    keep = val > 1e-10
    x = vec[:, keep] * np.sqrt(val[keep])
    for j in range(x.shape[1]):
        flip = x[:, j] @ reference[:, j] < 0 if reference is not None and j < reference.shape[1] else x[np.argmax(abs(x[:, j])), j] < 0
        if flip:
            x[:, j] *= -1
    return x, val


def fold_geometry(d_train, d_test_train=None):
    """Training-only MDS and Gower projection with signed residual correction.

    Add test_offset to squared projection error in x_test. This evaluates
    signed source-distance reconstruction to nodes in the positive training
    affine span; it does not discard a test observation's orthogonal component.
    The offset is constant across all candidate models within this fold.
    """
    rowmean = d_train.mean(0)
    grand = float(d_train.mean())
    g = -.5 * (d_train - rowmean - rowmean[:, None] + grand)
    val, vec = np.linalg.eigh(g)
    ix = np.argsort(val)[::-1]
    val, vec = val[ix], vec[:, ix]
    keep = val > 1e-10
    basis, roots = vec[:, keep], np.sqrt(val[keep])
    x = basis * roots
    info = dict(positive_dimensions=int(keep.sum()),
                negative_eigenvalue_mass=float(-val[val < 0].sum() / abs(val).sum()),
                smallest_retained_eigenvalue=float(val[keep][-1]))
    if d_test_train is None:
        return dict(x_train=x, x_test=None, test_offset=None, info=info)
    cross = -.5 * (d_test_train - rowmean[None, :] - d_test_train.mean(1)[:, None] + grand)
    test = cross @ (basis / roots)
    radius = d_test_train.mean(1) - .5 * grand
    offset = radius - np.sum(test * test, axis=1)
    return dict(x_train=x, x_test=test, test_offset=offset, info=info)


def make_folds(blocks, n_splits=10, seed=SEED):
    """Return train/test row-index pairs; no physical sheet crosses a fold."""
    blocks = np.asarray(blocks)
    unique = np.array(sorted(set(blocks)))
    np.random.default_rng(seed).shuffle(unique)
    output = []
    for held in np.array_split(unique, min(n_splits, len(unique))):
        test = np.isin(blocks, held)
        output.append((np.flatnonzero(~test), np.flatnonzero(test)))
    return output


def sheet_weights(blocks):
    blocks = np.asarray(blocks)
    weights = np.array([1. / np.sum(blocks == b) for b in blocks])
    return weights / weights.sum()


def summarize_errors(errors, blocks):
    errors, blocks = np.asarray(errors), np.asarray(blocks)
    per = {b:float(errors[blocks == b].mean()) for b in sorted(set(blocks))}
    return dict(equal_bifolium_mean=float(np.mean(list(per.values()))),
                page_mean=float(errors.mean()), minimum_page_score=float(errors.min()),
                per_bifolium=per)


def hash_sources(paths):
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in dict.fromkeys(paths)}


def verify_geometry():
    rng = np.random.default_rng(SEED)
    x = rng.normal(size=(25, 7))
    d = ((x[:, None] - x[None, :]) ** 2).sum(2)
    tr, te = np.arange(15), np.arange(15, 25)
    out = fold_geometry(d[np.ix_(tr, tr)], d[np.ix_(te, tr)])
    weights = rng.dirichlet(np.ones(len(tr)), size=4)
    node = weights @ out['x_train']
    got = ((out['x_test'][:, None] - node[None, :]) ** 2).sum(2) + out['test_offset'][:, None]
    expected = ((x[te, None] - (weights @ x[tr])[None, :]) ** 2).sum(2)
    np.testing.assert_allclose(got, expected, atol=1e-10, rtol=0)
    # A truly unseen dimension must remain in held-out reconstruction error.
    x[te, -1] += 3
    x[tr, -1] = 0
    d = ((x[:, None] - x[None, :]) ** 2).sum(2)
    out = fold_geometry(d[np.ix_(tr, tr)], d[np.ix_(te, tr)])
    got = ((out['x_test'][:, None] - (weights @ out['x_train'])[None, :]) ** 2).sum(2) + out['test_offset'][:, None]
    expected = ((x[te, None] - (weights @ x[tr])[None, :]) ** 2).sum(2)
    np.testing.assert_allclose(got, expected, atol=1e-10, rtol=0)
    return float(abs(got - expected).max())


if __name__ == '__main__':
    data = load_data()
    print('Paragraphs', len(data['meta']), 'bifolia', len(set(data['blocks'])), 'pools', len(data['pool_meta']))
    print('Held-out geometry maximum Euclidean check error', verify_geometry())
