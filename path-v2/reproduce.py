"""Path collection v2: conditional pairwise consensus from frozen ensembles.

No trajectory refitting. Peak RAM <250 MB; generated data <5 MB.
Published copy: python path-v2/reproduce.py --site-root . --out reproduced-v2
Workspace: python apps/clustermatrix/build_path_v2.py --out tmp/path-v2
"""
from pathlib import Path
from functools import lru_cache
import argparse
import csv
import hashlib
import io
import json
import numpy as np

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
TOL = 1e-7
DEFAULT = 'v2_default'
VARIANTS = [
    dict(id=DEFAULT, name='V2 default · overall consensus', conditions=None),
    dict(id='v2_3a_3b__4a_4c', name='V2 · 3a before 3b; 4a before 4c', conditions=[-1, -1]),
    dict(id='v2_3a_3b__4c_4a', name='V2 · 3a before 3b; 4c before 4a', conditions=[-1, 1]),
    dict(id='v2_3b_3a__4a_4c', name='V2 · 3b before 3a; 4a before 4c', conditions=[1, -1]),
]
COHORTS = {'pages': 'paragraph207', 'bifolia': 'pool50'}
OUTPUTS = {'path-v2/collection.json', 'path-v2/pages.csv', 'path-v2/bifolia.csv',
    'path-v2/selected-fits.csv', 'path-v2/README.md', 'path-v2/reproduce.py', 'path-v2/.gitattributes'}


def encoded(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(',', ':'), allow_nan=False)+'\n').encode('utf-8')


def csv_bytes(rows, fields=None):
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=fields or list(rows[0]), lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return out.getvalue().encode('utf-8')


def sign(run, key):
    return run['group_readout_all']['contrasts'][key]['sign']


def select(runs, conditions):
    eligible = [r for r in runs if r['converged'] and sign(r, '1a|1b') == -1
        and (conditions is None or [sign(r, '3a|3b'), sign(r, '4a|4c')] == conditions)]
    return [min([r for r in eligible if r['replicate'] == mask],
        key=lambda r: (r['objective'], r['start_index'])) for mask in sorted({r['replicate'] for r in eligible})]


def aggregate(runs, design):
    n = len(design['meta']); masks = np.asarray(design['masks'], bool)
    after = np.zeros((n, n), dtype=np.int32); ties = np.zeros_like(after); joint = np.zeros_like(after)
    for r in runs:
        keep = masks[r['replicate']]; both = keep[:, None] & keep[None, :]
        t = np.asarray(r['positions']); delta = t[:, None]-t[None, :]
        after += both & (delta > TOL); ties += both & (np.abs(delta) <= TOL); joint += both
    assert np.array_equal(after+after.T+ties, joint)
    assert np.all(joint > 0), 'Some peers have no joint evidence; do not silently impute ranks.'
    return (after+.5*ties)/joint, joint


def source_team(site_root=None):
    root = Path(site_root)/'chronology-stability/replication' if site_root else ROOT
    return root/'solve/botanical_enrichment/team_golf'


def build(site_root=None):
    team = source_team(site_root); selected_rows = []; sources = {}; cohorts = {}
    for unit, cohort in COHORTS.items():
        settings = []
        for rep, folder in [('RF', 'golf_rf1b_stability_20261016'), ('ZL', 'golf_published_stability_20261016')]:
            index = json.loads((team/folder/'ensemble_index.json').read_text())
            for code in ['010', '050']:
                name = f'ensemble_{cohort}_point_l{code}.json'; path = team/folder/name
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                assert digest == index['result_sha256'][name], path
                sources[f'{folder}/{name}'] = digest
                d = json.loads(path.read_text()); assert d['status'] == 'complete'
                settings.append((rep, float(d['design']['penalty']), d))
        meta = [{k: m[k] for k in ['id', 'group', 'bifolium']} for m in settings[0][2]['design']['meta']]
        for _, _, d in settings:
            assert meta == [{k: m[k] for k in ['id', 'group', 'bifolium']} for m in d['design']['meta']]
        n = len(meta); groups = np.array([m['group'] for m in meta]); paths = {}
        for variant in VARIANTS:
            fractions = []; support = []; retained = np.zeros(n, dtype=int)
            for rep, penalty, d in settings:
                runs = select(d['runs'], variant['conditions'])
                fraction, joint = aggregate(runs, d['design']); fractions.append(fraction)
                retained += np.diag(joint)
                support.append(dict(representation=rep, penalty=penalty, eligible_samples=len(runs),
                    omitted_samples=100-len(runs), minimum_joint_samples=int(joint[np.triu_indices(n, 1)].min())))
                selected_rows.extend(dict(path_id=variant['id'], unit=unit, representation=rep, penalty=penalty,
                    replicate=r['replicate'], start_index=r['start_index'], start_id=r['start_id'], objective=r['objective']) for r in runs)
            # Match the inspector: each of the four settings contributes 1/4,
            # after normalization by its own jointly retained fit count.
            pairwise = np.mean(fractions, axis=0)
            np.testing.assert_allclose(pairwise+pairwise.T, 1, atol=1e-14)
            score = (pairwise.sum(1)-.5)/(n-1)
            ranks = np.empty(n, int); ranks[np.argsort(score, kind='stable')] = np.arange(1, n+1)
            medians = {g: float(np.median(score[groups == g])) for g in sorted(set(groups))}
            signature = [int(np.sign(medians[a]-medians[b])) for a, b in [('3a', '3b'), ('4a', '4c')]]
            paths[variant['id']] = dict(scores=score.tolist(), ranks=ranks.tolist(),
                retained_fits=retained.tolist(), support=support, group_score_medians=medians,
                resulting_group_order=signature,
                pairwise_after_half_tie=pairwise[np.triu_indices(n, 1)].tolist())
        assert paths[DEFAULT]['resulting_group_order'] == [1, 1]
        cohorts[unit] = dict(meta=meta, paths=paths)
    return dict(version=2, classification_version='Golf 2.1', default=DEFAULT,
        variants=VARIANTS, omitted_conditional_variant=dict(conditions=[1, 1], reason='Same ending directions as the overall default, which remains unconditioned.'),
        presets=dict(representations=['RF/basic EVA characters', 'ZL/Zattera'], representation_weights=[.5, .5],
            penalties=[.01, .05], penalty_weights=[.5, .5], scope='best eligible fit per sample',
            early_a_filter='1a before 1b', sections='all, including unassigned', page_sampling='individual pages',
            bifolium_sampling='whole pooled bifolia', tie_tolerance=TOL),
        score_definition='Mean across all other cohort members of same-fit after + half-tie frequency; four settings normalized separately then equally weighted. Lower is earlier. No min-max rescaling.',
        variant_selection='Filter full-reference group-median directions first, then minimize objective within each representation/penalty/sample. Missing eligible samples are reported, never replaced.',
        interpretation='Descriptive consensus, not a probability, elapsed time, unique chronology, or one actual fitted curve. Page and bifolium scores use different peer populations.',
        source_sha256=sources, cohorts=cohorts), selected_rows


def wide_rows(collection, unit):
    c = collection['cohorts'][unit]; out = []
    for i, meta in enumerate(c['meta']):
        row = dict(id=meta['id'], golf_2_1=meta['group'], bifolium=meta['bifolium'],
            path_collection_version=2, default_path_id=DEFAULT)
        for variant in VARIANTS:
            key = variant['id']; p = c['paths'][key]
            row.update({key+'_score': p['scores'][i], key+'_rank': p['ranks'][i], key+'_retained_fits': p['retained_fits'][i]})
        out.append(row)
    return out


@lru_cache(maxsize=1)
def collection():
    return build()


def assignment_fields(collection, page, sheet):
    row = dict(path_collection_version=2, default_path_id=DEFAULT)
    for unit, key in [('pages', page), ('bifolia', sheet)]:
        c = collection['cohorts'][unit]; lookup = {m['id']: i for i, m in enumerate(c['meta'])}; i = lookup.get(key)
        prefix = 'page' if unit == 'pages' else 'bifolium'
        row[prefix+'_path_score'] = '' if i is None else c['paths'][DEFAULT]['scores'][i]
        row[prefix+'_path_rank'] = '' if i is None else c['paths'][DEFAULT]['ranks'][i]
        for variant in VARIANTS[1:]:
            p = c['paths'][variant['id']]
            row[prefix+'_'+variant['id']+'_score'] = '' if i is None else p['scores'][i]
            row[prefix+'_'+variant['id']+'_rank'] = '' if i is None else p['ranks'][i]
    return row


def display_data(base):
    """Smooth guides at fixed consensus scores, never new ordering estimates.

    Local linear regression in the eight existing display coordinates; fixed
    bandwidth 0.06 on the score axis, equal observation weights, 160 samples.
    Both representations contribute to scores; the preserved maps are ZL maps.
    """
    c, _ = collection()
    points = base['points']; para = [i for i,p in enumerate(points) if p['kind'] == 'para']
    pools = [i for i,p in enumerate(points) if p['kind'] == 'pool']
    lookup = {p['id']: i for i,p in enumerate(points)}

    def guide(scores, xyz, use):
        t = np.array([scores[i] for i in use]); y = np.asarray(xyz)[use]
        grid = np.linspace(t.min(), t.max(), 160); delta = t[None, :]-grid[:, None]
        w = np.exp(-.5*(delta/.06)**2)
        a = w.sum(1); b = (w*delta).sum(1); d = (w*delta**2).sum(1)
        yy = w@y; dy = (w*delta)@y
        nodes = (d[:, None]*yy-b[:, None]*dy)/(a*d-b*b)[:, None]
        assert np.isfinite(nodes).all()
        return dict(nodes=nodes.tolist(), positions=scores, count=len(use),
            kind='display guide at fixed consensus scores', parameter_scores=grid.tolist())

    display = {}
    for variant in VARIANTS:
        key = variant['id']; scores = [None]*len(points)
        for unit in COHORTS:
            co = c['cohorts'][unit]
            for m, value in zip(co['meta'], co['paths'][key]['scores']):
                scores[lookup[m['id']]] = value
        pp = np.array([scores[i] for i in para]); psd = float(np.std(pp, ddof=1))
        pool_ranks = dict(zip([m['id'] for m in c['cohorts']['bifolia']['meta']], c['cohorts']['bifolia']['paths'][key]['ranks']))
        rows = []
        for old in base['sheets']:
            row = dict(old); pool = lookup[row['bifolium']]
            ii = [i for i in para if points[i]['sheet'] == row['bifolium']]
            row.update(rank=pool_ranks[row['bifolium']], path_position=scores[pool], path_page_sd_normalized=float(np.std([scores[i] for i in ii], ddof=1)/psd) if len(ii)>1 else 0)
            rows.append(row)
        rows.sort(key=lambda r:r['rank'])
        total = float(np.sum((pp-pp.mean())**2)); within = 0.
        for sheet in sorted({points[i]['sheet'] for i in para}):
            vals = np.array([scores[i] for i in para if points[i]['sheet'] == sheet])
            within += float(np.sum((vals-vals.mean())**2))
        variation = [dict(r) for r in base['variation'][:8]] + [dict(label='V2 consensus score', fraction=within/total),
            dict(label='Held-out v2 score', fraction=None)]
        curves = {mode:dict(proposed=guide(scores, model['coordinates'], pools),
            pages=guide(scores, model['coordinates'], para),
            paragraphs_and_pools=guide(scores, model['coordinates'], para+pools)) for mode,model in base['modes'].items()}
        display[key] = dict(positions=scores, sheets=rows, variation=variation, curves=curves)
    display['legacy_proposed'] = dict(positions=base['modes']['combined']['curves']['proposed']['positions'],
        sheets=base['sheets'], variation=base['variation'], curves={})
    return dict(collection=c, default=DEFAULT, options=VARIANTS+[
        dict(id='legacy_proposed', name='Legacy · original proposed path (v1)', conditions=None)],
        display=display, guide_method=dict(method='Gaussian local linear regression at fixed consensus score',
            bandwidth=.06, points=160, dimensions=8, reference_geometry='unchanged ZL shared MDS',
            interpretation='Visual summary, not an ElPiGraph fit or nearest-curve projection. It cannot change consensus scores.'))


def package():
    c, selected = collection()
    files = {'path-v2/collection.json': encoded(c), 'path-v2/selected-fits.csv': csv_bytes(selected),
        'path-v2/README.md': (APP/'path-v2-methods.md').read_bytes(), 'path-v2/reproduce.py': Path(__file__).read_bytes(),
        'path-v2/.gitattributes': b'* -text\n'}
    for unit in COHORTS:
        files['path-v2/'+unit+'.csv'] = csv_bytes(wide_rows(c, unit))
    return files


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--site-root', type=Path); p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(); c, rows = build(args.site_root); args.out.mkdir(parents=True, exist_ok=True)
    (args.out/'collection.json').write_bytes(encoded(c))
    (args.out/'selected-fits.csv').write_bytes(csv_bytes(rows))
    for unit in COHORTS:
        (args.out/(unit+'.csv')).write_bytes(csv_bytes(wide_rows(c, unit)))
    print(json.dumps({'status': 'complete', 'default': c['default'], 'cohorts': {k:len(v['meta']) for k,v in c['cohorts'].items()},
        'variant_group_orders': {k:{p:x['resulting_group_order'] for p,x in v['paths'].items()} for k,v in c['cohorts'].items()}}))
