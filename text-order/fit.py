"""Text -> transition JS -> positive classical MDS -> fresh ElPiGraph curves.

No Golf labels, saved fits or group-order filters enter fitting or selection.
An undirected curve and its reverse are equally valid. See README.md.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from collections import Counter
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import numpy as np
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
import runtime as R

TOL = 1e-7


def read_csv(path):
    """One row per text line; repeated IDs preserve independent line boundaries."""
    lines, meta = {}, {}
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        if not {'id', 'text'} <= set(reader.fieldnames or []):
            raise ValueError('CSV requires id,text; optional block,group columns')
        for row in reader:
            uid = row['id'].strip()
            if not uid or row['text'] is None:
                raise ValueError('Each row needs a nonempty id and a text cell')
            m = dict(id=uid, bifolium=(row.get('block') or uid).strip(), group=(row.get('group') or '').strip())
            if uid in meta and meta[uid] != m:
                raise ValueError('Conflicting metadata for ' + uid)
            meta[uid] = m
            lines.setdefault(uid, []).append(row['text'].split())
    # Input row order must not supply a proposed document order.
    return [meta[k] for k in sorted(meta)], lines


def count_csv(meta, lines, mode):
    from voynich_toolkit.glyph_units import units
    def glyphs(word):
        if mode == 'characters':
            return tuple(word)
        return tuple('other' if u.rare else u.text for u in units(word, mode))
    counts = {m: {r['id']: Counter() for r in meta} for m in ('word', 'line')}
    words, excluded = Counter(), Counter()
    for row in meta:
        uid = row['id']
        for tokens in lines[uid]:
            if any('^' in t or '$' in t for t in tokens):
                raise ValueError('^ and $ are reserved boundary markers')
            good = ['?' not in t for t in tokens]
            words[uid] += sum(good)
            excluded[uid] += len(tokens) - sum(good)
            for token, ok in zip(tokens, good):
                if ok:
                    seq = glyphs(token)
                    counts['word'][uid].update(zip(('^',) + seq, seq + ('$',)))
            start = None
            for i in range(len(tokens) + 1):
                ok = i < len(tokens) and good[i]
                if ok and start is None:
                    start = i
                if not ok and start is not None:
                    seq = glyphs(''.join(tokens[start:i]))
                    if seq:
                        counts['line'][uid].update(zip(seq, seq[1:]))
                        if start == 0:
                            counts['line'][uid][('^', seq[0])] += 1
                        if i == len(tokens):
                            counts['line'][uid][(seq[-1], '$')] += 1
                    start = None
    if any(words[m['id']] == 0 for m in meta):
        raise ValueError('Every observation must contain readable words')
    return counts, dict(readable_words=dict(words), excluded_words=dict(excluded))


def inputs(path, fmt, mode, metadata=None, pooled=False):
    prep = R.preprocessing()
    if fmt == 'csv':
        meta, lines = read_csv(path)
        counts, census = count_csv(meta, lines, mode)
    else:
        if mode == 'characters':
            raise ValueError('IVTFF requires --units raw or zattera')
        from voynich_toolkit.source_lattice import SourceLattice
        doc = SourceLattice(path, verify=False)
        ids = sorted({line.locus.page for line in doc.stream('all_spaces', 'preferred_first') if line.locus.category == 'para'})
        if metadata:
            # Only physical membership and reporting labels are used here.
            rows = json.loads(metadata.read_text(encoding='utf-8'))['pages']
            lookup = {m['id']: m for m in rows}
            if set(ids) != set(lookup):
                raise ValueError('IVTFF paragraph IDs must match metadata exactly')
            meta = [dict(id=uid, bifolium=lookup[uid]['bifolium'], group=lookup[uid].get('quire', '')) for uid in ids]
        else:
            if pooled:
                raise ValueError('Pooling IVTFF requires --metadata with physical membership')
            meta = [dict(id=uid, bifolium=uid, group='') for uid in ids]
        counts, census = prep.count_source(path, meta, mode, verify=False)
    alphabet = sorted({'^', '$'} | {g for m in counts.values() for c in m.values() for pair in c for g in pair})
    arrays = {mode: prep.dense(c, meta, alphabet) for mode, c in counts.items()}
    if pooled:
        blocks = np.array([m['bifolium'] for m in meta])
        new = []
        for block in sorted(set(blocks)):
            labels = sorted({m['group'] for m in meta if m['bifolium'] == block and m['group']})
            new.append(dict(id=block, bifolium=block, group=';'.join(labels)))
        arrays = {mode: np.array([a[blocks == m['id']].sum(0) for m in new]) for mode, a in arrays.items()}
        meta = new
    if len(meta) < 4:
        raise ValueError('Need at least four observations after pooling')
    d = sum(prep.distances(a) for a in arrays.values()) / 2
    if not np.any(d > 1e-12):
        raise ValueError('No usable transition differences; cannot infer an order')
    return meta, d, census, alphabet


def starts(x, weights, count, seed):
    """Declared unlabeled initialization design, not an ElPiGraph default."""
    mean = weights @ x
    _, _, axes = np.linalg.svd(np.sqrt(weights[:, None]) * (x - mean), full_matrices=False)
    for axis in axes:
        if axis[np.argmax(abs(axis))] < 0:
            axis *= -1
    rng = np.random.default_rng(seed)
    result = []
    for i in range(count):
        axis = axes[0] if i < 2 else rng.normal(size=x.shape[1])
        axis = axis / np.linalg.norm(axis)
        extent = np.sqrt(weights @ ((x - mean) @ axis) ** 2)
        initial = mean + np.array([-1., 0., 1.])[:, None] * extent * axis
        if i == 1 and len(axes) > 1:
            initial[1] += np.sqrt(weights @ ((x - mean) @ axes[1]) ** 2) * axes[1]
        result.append(initial)
    return result, (x - mean) @ axes[0]


def candidate(core, x, w, initial, penalty, nodes):
    grown = core.grow_curve(x, w, initial, lam=penalty, n_nodes=nodes)
    route = core.chain_order(grown['edges'], len(grown['nodes']))
    return core.refine_curve(x, w, grown['nodes'][route], lam=penalty)


def best_fit(core, x, w, initial, penalty, nodes):
    attempts, best = [], None
    for i, start in enumerate(initial):
        model = candidate(core, x, w, start, penalty, nodes)
        attempts.append(dict(start=i, objective=model['objective'], converged=model['converged']))
        if model['converged'] and (best is None or model['objective'] < best['objective']):
            best = dict(model, start=i)
    if best is None:
        raise ValueError('No converged candidate; no order silently substituted')
    return best, attempts


def align(t, reference):
    direct = np.sum((t - reference) ** 2)
    reverse = np.sum((1 - t - reference) ** 2)
    return (1 - t if reverse < direct else t), abs(direct - reverse)


def consensus(positions, masks):
    n = len(positions[0])
    joint = np.zeros((n, n), dtype=int)
    after = np.zeros((n, n), dtype=float)
    for t, keep in zip(positions, masks):
        both = keep[:, None] & keep[None, :]
        delta = t[:, None] - t[None, :]
        joint += both
        after += both * ((delta > TOL) + .5 * (abs(delta) <= TOL))
    if np.any(joint == 0):
        raise ValueError('Some pairs have no jointly retained fits; increase --samples')
    return after / joint, joint


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--format', choices=['csv', 'ivtff'], default='csv')
    p.add_argument('--units', choices=['characters', 'raw', 'zattera'], default='characters')
    p.add_argument('--metadata', type=Path, help='IVTFF physical membership; only pages id/bifolium/quire read')
    p.add_argument('--pool', action='store_true', help='Sum transition counts within physical blocks')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--samples', type=int, default=100)
    p.add_argument('--starts', type=int, default=8)
    p.add_argument('--nodes', type=int, default=28)
    p.add_argument('--penalties', type=float, nargs='+', default=[.01, .05])
    p.add_argument('--inclusion', type=float, default=.9)
    p.add_argument('--seed', type=int, default=20261028)
    args = p.parse_args()
    if args.out.exists():
        p.error('output already exists; use a new directory')
    if args.samples < 1 or args.starts < 1 or args.nodes < 3 or not 0 < args.inclusion <= 1:
        p.error('positive samples/starts, nodes >= 3 and 0 < inclusion <= 1 required')
    if any(not np.isfinite(v) or v <= 0 for v in args.penalties) or len(set(args.penalties)) != len(args.penalties):
        p.error('penalties must be distinct finite positive numbers')
    engine = R.engine()
    core = engine.core
    began = time.time()
    with threadpool_limits(1):
        meta, d, census, alphabet = inputs(args.input, args.format, args.units, args.metadata, args.pool)
        x, eigenvalues = core.COMMON.embed(d)
        weights = core.COMMON.sheet_weights([m['bifolium'] for m in meta])
        initial, pc1 = starts(x, weights, args.starts, args.seed)
        # A single reference resolves arbitrary reversal before averaging runs.
        # Its penalty is declared, not selected by comparing unlike objectives.
        reference, reference_attempts = best_fit(core, x, weights, initial, args.penalties[0], args.nodes)
        ref = core.project_curve(x, reference)['positions']
        if np.dot(ref - ref.mean(), pc1 - pc1.mean()) < 0:
            ref = 1 - ref
        if np.ptp(ref) <= TOL:
            raise ValueError('Reference curve has no resolved positions')
        masks = np.random.RandomState(args.seed).random_sample((args.samples, len(meta))) <= args.inclusion
        if np.any(masks.sum(1) < 3):
            raise ValueError('A sample has fewer than three observations; increase inclusion or sample size')
        fractions, supports, fits, selected = [], [], [], []
        for penalty in args.penalties:
            positions = []
            for sample, mask in enumerate(masks):
                model, attempts = best_fit(core, x[mask], weights[mask] / weights[mask].sum(), initial, penalty, args.nodes)
                t, margin = align(core.project_curve(x, model)['positions'], ref)
                positions.append(t)
                record = dict(penalty=penalty, sample=sample, start=model['start'], objective=model['objective'],
                    converged=model['converged'], orientation_margin=margin, positions=t.tolist(),
                    retained=np.flatnonzero(mask).tolist(), attempts=attempts)
                fits.append(record)
                selected.append({k: record[k] for k in ('penalty', 'sample', 'start', 'objective', 'orientation_margin')})
                print(f'lambda={penalty}: {sample + 1}/{args.samples} samples', flush=True)
            fraction, joint = consensus(positions, masks)
            fractions.append(fraction)
            supports.append(dict(penalty=penalty, minimum_joint_samples=int(joint.min())))
        pairwise = np.mean(fractions, axis=0)
        scores = (pairwise.sum(1) - .5) / (len(meta) - 1)
        rows = [dict(m, score=float(scores[i]), reverse_score=float(1 - scores[i]),
            rank=float(rankdata(scores, method='average')[i]),
            reverse_rank=float(rankdata(-scores, method='average')[i])) for i, m in enumerate(meta)]
        labels = np.array([m['group'] for m in meta])
        groups = {g: dict(n=int(sum(labels == g)), median_score=float(np.median(scores[labels == g])))
            for g in sorted(set(labels)) if g}
        # Reporting only: this executes after every curve has been fitted/selected.
        comparisons = {a + '|' + b: float(pairwise[np.ix_(labels == a, labels == b)].mean())
            for a in groups for b in groups if a < b}
        audit = dict(status='complete', method='ElPiGraph 0.3.4 principal-curve growth plus documented fixed-graph refinement',
            config={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
            metadata_sha256=None if not args.metadata else hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
            source_sha256={Path(f).name: hashlib.sha256(Path(f).read_bytes()).hexdigest()
                for f in [__file__, R.__file__, core.__file__, core.COMMON.__file__, core.MODELS.__file__]},
            versions={n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'elpigraph-python', 'numba', 'scikit-learn']},
            observations=len(meta), positive_dimensions=x.shape[1],
            negative_eigenvalue_mass=float(-eigenvalues[eigenvalues < 0].sum() / abs(eigenvalues).sum()),
            census=census, alphabet=alphabet, group_summary=groups, group_pairwise_after=comparisons,
            direction='Arbitrary reference convention. Reversing all scores is equally valid; no historical direction inferred.',
            selection='Minimum converged objective per penalty/sample, no labels or order filters; settings equally weighted.',
            reference_penalty=args.penalties[0], reference_attempts=reference_attempts,
            reference_positions=ref.tolist(), support=supports, fits=fits, seconds=time.time() - began)
        args.out.mkdir(parents=True)
        with (args.out / 'scores.csv').open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda r: (r['rank'], r['id'])))
        with (args.out / 'pairwise.csv').open('w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['id'] + [m['id'] for m in meta])
            writer.writerows([[m['id']] + pairwise[i].tolist() for i, m in enumerate(meta)])
        engine.dump(args.out / 'run.json', audit)
        print(f'{len(meta)} observations; fresh ElPiGraph results in {args.out}. Direction is arbitrary.')


if __name__ == '__main__':
    main()
