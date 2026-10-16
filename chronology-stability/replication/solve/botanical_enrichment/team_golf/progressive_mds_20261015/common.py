"""Frozen source data and nested cohorts for the progressive Golf MDS audit."""
from pathlib import Path
import csv
import json
import numpy as np

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent
ROOT = TEAM.parents[2]
DATA = ROOT / 'apps/clustermatrix/research-data'


def rows(path):
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def embed(d, reference=None):
    gram = -.5 * (d - d.mean(0) - d.mean(1)[:, None] + d.mean())
    val, vec = np.linalg.eigh(gram)
    order = np.argsort(val)[::-1]
    val, vec = val[order], vec[:, order]
    x = vec[:, val > 1e-12] * np.sqrt(val[val > 1e-12])
    for k in range(x.shape[1]):
        flip = (x[:, k] @ reference[:, k] < 0) if reference is not None and k < reference.shape[1] else x[np.argmax(abs(x[:, k])), k] < 0
        if flip:
            x[:, k] *= -1
    return x, val


def load_data():
    profiles = rows(DATA / 'profiles.csv')
    byid = {p['shared_mds_id']: p for p in profiles if p['shared_mds_id']}
    viewer = rows(DATA / 'mds/viewer-branches.csv')
    selected = [r for r in viewer if r['version'] == 'golf21' and r['branch'] == 'root' and r['view'] == 'pages']
    b_rows = {r['viewer_id']: r for r in viewer if r['version'] == 'golf21' and r['branch'] == 'b' and r['view'] == 'pages'}
    meta = [dict(id=r['viewer_id'], profile_id=byid[r['viewer_id']]['profile_id'],
                 group=r['membership_in_view'], kind=r['kind'], bifolium=r['bifolium'],
                 words=int(r['readable_words']), section=byid[r['viewer_id']]['section'],
                 hand=byid[r['viewer_id']]['lfd_hand']) for r in selected]
    sources = [DATA / 'profiles.csv', DATA / 'mds/viewer-branches.csv']
    ds = {}
    for mode in ['word', 'line', 'combined']:
        path = DATA / f'pairwise/shared-{mode}-js.csv'
        table = rows(path)
        lookup = {r['profile_id']: r for r in table}
        ds[mode] = np.array([[float(lookup[a['profile_id']][b['profile_id']]) for b in meta] for a in meta])
        sources.append(path)
    base_indices = np.array([i for i, m in enumerate(meta) if m['id'] in b_rows])
    masks = [np.array([i in base_indices for i in range(len(meta))])]
    definitions = [('b', 'B alone', 'Original B-related individual-page cohort'),
                   ('physical_pharma', '+ physical Pharma', 'Section P paragraphs (16 units)'),
                   ('all_2a', '+ remaining 2a', 'Remaining Golf 2a units (10)'),
                   ('transitional', '+ Transitional A', 'Golf 1.5a (20)'),
                   ('early_1b', '+ Early A 1b', 'Golf 1b (31)'),
                   ('early_1a', '+ Early A 1a', 'Golf 1a (39)'),
                   ('whole', '+ remaining units', 'Remaining 23 units; whole manuscript')]
    masks.append(masks[-1] | np.array([m['section'] == 'P' and m['kind'] == 'para' for m in meta]))
    for group in ['2a', '1.5a', '1b', '1a']:
        masks.append(masks[-1] | np.array([m['group'] == group for m in meta]))
    masks.append(np.ones(len(meta), dtype=bool))
    stages = [dict(key=k, label=l, addition=a, indices=np.flatnonzero(mask).tolist()) for (k, l, a), mask in zip(definitions, masks)]
    assert [len(s['indices']) for s in stages] == [84, 100, 110, 130, 161, 200, 223]
    reference = np.array([[float(r[f'MDS{k}']) for k in range(1, 5)] for r in selected])
    b_reference = np.array([[float(b_rows[meta[i]['id']][f'MDS{k}']) for k in range(1, 5)] for i in base_indices])
    seqpath = TEAM / 'nonlinear_sequence_20261010/plot_data.json'
    seq = json.loads(seqpath.read_text())
    positions = dict(zip([p['id'] for p in seq['points']], seq['models']['combined']['proposed_positions']))
    global_t = np.array([positions[m['id']] for m in meta])
    sources.append(seqpath)
    return dict(meta=meta, ds=ds, stages=stages, base_indices=base_indices,
                global_t=global_t, reference=reference, b_reference=b_reference, sources=sources)
