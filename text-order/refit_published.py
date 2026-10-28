"""Raw ZL/RF text -> fresh documented ElPiGraph ensembles -> published v2 scores.

The explicit published design includes labels and starting orders. This is
reproduction of that conditional analysis, not label-free order discovery.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import numpy as np
from threadpoolctl import threadpool_limits
import runtime as R


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(design, zl, rf):
    """No frozen distances, coordinates, nodes, positions or fitted results."""
    prep = R.preprocessing()
    core = R.engine().core
    pages, pools = design['pages'], design['bifolia']
    blocks = np.array([m['bifolium'] for m in pages])
    sheets = np.array(sorted({m['bifolium'] for m in pools}))
    sheet_draws = np.random.RandomState(2026101602).random_sample((100, len(sheets))) <= .9
    point_draws = np.random.RandomState(2026101601).random_sample((100, len(pages))) <= .9
    contexts, census = {}, {}
    for rep, path, units in [('ZL', zl, 'zattera'), ('RF', rf, 'raw')]:
        counts, info = prep.count_source(path, pages, units, verify=False)
        census[rep] = info
        alphabet = design['zl_alphabet'] if rep == 'ZL' else sorted(
            {'^', '$'} | {g for mode in counts.values() for c in mode.values() for pair in c for g in pair})
        distances = {'pages': [], 'bifolia': []}
        for mode in ('word', 'line'):
            a = prep.dense(counts[mode], pages, alphabet)
            b = np.array([a[blocks == m['bifolium']].sum(0) for m in pools])
            distances['pages'].append(prep.distances(a))
            distances['bifolia'].append(prep.distances(b))
        for unit, meta in [('pages', pages), ('bifolia', pools)]:
            d = sum(distances[unit]) / 2
            x, _ = core.COMMON.embed(d)
            bs = np.array([m['bifolium'] for m in meta])
            w = core.COMMON.sheet_weights(bs)
            labels = np.array([m['group'] for m in meta])
            centers = {g: w[labels == g] @ x[labels == g] / w[labels == g].sum() for g in core.ORDER}
            masks = point_draws if unit == 'pages' else sheet_draws[:, [list(sheets).index(b) for b in bs]]
            for penalty in (.01, .05):
                key = f"{'paragraph207' if unit == 'pages' else 'pool50'}_point_l{int(penalty * 1000):03d}"
                starts = []
                for s in design['starts'][key]:
                    initial = np.array([centers[g] for g in s['anchor_order']])
                    starts.append(dict(s, initial=initial))
                contexts[rep + '/' + key] = dict(key=key, cohort=key.split('_')[0], sampler='point',
                    penalty=penalty, meta=meta, x=x, weights=w, blocks=bs, labels=labels,
                    starts=starts, masks=masks, distance=d)
    return contexts, census


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design', type=Path, default=R.HERE / 'published-design.json')
    p.add_argument('--zl', type=Path, default=R.ROOT / 'ZL_ivtff_raw.txt')
    p.add_argument('--rf', type=Path, default=R.RF / 'RF1b-er.txt')
    p.add_argument('--out', type=Path, required=True, help='New directory; never overwrites a release')
    p.add_argument('--replicates', type=int, default=100, help='1..100; fewer is a smoke test only')
    p.add_argument('--start-limit', type=int, default=8, help='1..8; fewer is a smoke test only')
    args = p.parse_args()
    if not 1 <= args.replicates <= 100 or not 1 <= args.start_limit <= 8:
        p.error('replicates must be 1..100 and start-limit must be 1..8')
    if args.out.exists():
        p.error('output already exists; use a new directory')
    design = json.loads(args.design.read_text(encoding='utf-8'))
    for rep, path in [('ZL', args.zl), ('RF', args.rf)]:
        if sha(path) != design['raw_sha256'][rep]:
            p.error(f'{rep} raw input differs from the published transcription; use fit.py for other text')
    engine = R.engine()
    with threadpool_limits(1):
        contexts, census = prepare(design, args.zl, args.rf)
        engine._STATE = contexts
        args.out.mkdir(parents=True)
        full = args.replicates == 100 and args.start_limit == 8
        began = time.time()
        folders = {'ZL': 'golf_published_stability_20261016', 'RF': 'golf_rf1b_stability_20261016'}
        outputs = {rep: {} for rep in folders}
        versions = {name: importlib.metadata.version(name) for name in
            ['numpy', 'scipy', 'elpigraph-python', 'numba', 'scikit-learn', 'threadpoolctl']}
        audit = dict(status='running', full_published_design=full, versions=versions,
            design_sha256=sha(args.design), source_sha256=sha(__file__), census=census,
            engine_sha256=sha(engine.__file__), core_sha256=sha(engine.core.__file__),
            planned_fits=len(contexts) * args.replicates * args.start_limit,
            note='Fresh growth with published group anchors, orientation and downstream filters.')
        engine.dump(args.out / 'refit-audit.json', audit)
        completed = 0
        for compound, c in contexts.items():
            rep, key = compound.split('/')
            destination = args.out / 'chronology-stability/replication/solve/botanical_enrichment/team_golf' / folders[rep]
            destination.mkdir(parents=True, exist_ok=True)
            runs = []
            for sample in range(args.replicates):
                for start in range(args.start_limit):
                    _, run = engine.fit_one((compound, sample, start))
                    runs.append(run)
                    completed += 1
                print(f'{compound}: {sample + 1}/{args.replicates}; {completed}/{audit["planned_fits"]} fresh fits', flush=True)
            name = 'ensemble_' + key + '.json'
            engine.dump(destination / name, dict(status='complete' if full else 'smoke-test',
                design=dict(meta=c['meta'], masks=c['masks'], penalty=c['penalty']), runs=runs))
            outputs[rep][name] = sha(destination / name)
        for rep, folder in folders.items():
            destination = args.out / 'chronology-stability/replication/solve/botanical_enrichment/team_golf' / folder
            engine.dump(destination / 'ensemble_index.json', dict(status='complete' if full else 'smoke-test',
                result_sha256=outputs[rep]))
        if full:
            script = R.SITE / 'path-v2/reproduce.py'
            if not script.exists():
                script = R.SITE / 'build_path_v2.py'
            summary = R.module('refit_v2_summary', script)
            collection, selected = summary.build(args.out)
            target = args.out / 'path-v2'
            target.mkdir()
            (target / 'collection.json').write_bytes(summary.encoded(collection))
            (target / 'selected-fits.csv').write_bytes(summary.csv_bytes(selected))
            for unit in summary.COHORTS:
                (target / (unit + '.csv')).write_bytes(summary.csv_bytes(summary.wide_rows(collection, unit)))
        audit.update(status='complete' if full else 'smoke-test', completed_fits=completed,
            seconds=time.time() - began, result_sha256=outputs)
        engine.dump(args.out / 'refit-audit.json', audit)
        print('Published-score refit complete.' if full else 'Smoke test complete; no published consensus claimed.')


if __name__ == '__main__':
    main()
