"""Actual ElPiGraph, nested sheet CV, and whole-manuscript topology controls.

Temporary dependencies live in _deps (remove after all team runs are done).
RAM budget <1 GB; persistent metrics/display <5 MB, no large arrays saved.
"""
from pathlib import Path
import contextlib
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import sys
import time
import warnings

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / '_deps'))
import elpigraph
from elpigraph.src.distutils import ComputePrimitiveGraphElasticEnergy

SEED = 20261015
CONFIGS = [dict(nodes=n, **{'lambda': lam}, mu=.1) for n in [8, 12] for lam in [.01, .05]]
FAMILIES = ['curve', 'tree', 'two_paths', 'plane']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


COMMON = module('whole_topology_common_for_models', HERE / 'common.py')


def project(x, model):
    """Per-point squared segment distance, selected edge and affine fraction."""
    x = np.asarray(x)
    if model['family'] == 'plane':
        delta = x - model['mean']
        prediction = model['mean'] + (delta @ model['axes'].T) @ model['axes']
        return np.sum((x - prediction) ** 2, axis=1), None, None
    nodes, edges = model['nodes'], model['edges']
    a = nodes[edges[:, 0]]; v = nodes[edges[:, 1]] - a
    l2 = np.sum(v * v, axis=1)
    dot = x @ v.T - np.sum(a * v, axis=1)
    frac = np.clip(dot / np.maximum(l2, 1e-15), 0, 1)
    d2 = np.maximum(np.sum(x*x, axis=1)[:, None] + np.sum(a*a, axis=1)[None, :]
                    - 2*x@a.T - 2*frac*dot + frac*frac*l2, 0)
    nearest = np.argmin(d2, axis=1)
    return d2[np.arange(len(x)), nearest], nearest, frac[np.arange(len(x)), nearest]


def project_errors(x, model):
    return project(x, model)[0]


def _initializations(x, weights):
    mean = weights @ x
    _, _, axes = np.linalg.svd(np.sqrt(weights[:, None]) * (x - mean), full_matrices=False)
    for row in axes[:2]:
        if row[np.argmax(abs(row))] < 0:
            row *= -1
    scores = (x - mean) @ axes[:2].T
    extent = np.sqrt(np.sum(weights * scores[:, 0] ** 2))
    line = mean + np.array([-1., 0., 1.])[:, None] * extent * axes[0]
    bent = line.copy()
    if axes.shape[0] > 1:
        bent[1] += np.sqrt(np.sum(weights * scores[:, 1] ** 2)) * axes[1]
    return [('weighted_PCA1', line), ('unlabeled_PCA_bend', bent)]


def _single_graph(x, blocks, family, config, point_weights=None):
    weights = COMMON.sheet_weights(blocks) if point_weights is None else np.asarray(point_weights) / np.sum(point_weights)
    fn = elpigraph.computeElasticPrincipalTree if family == 'tree' else elpigraph.computeElasticPrincipalCurve
    attempts = []
    for label, initial in _initializations(x, weights):
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured), warnings.catch_warnings(record=True) as caught:
            output = fn(np.asarray(x, dtype=float), NumNodes=int(config['nodes']),
                Lambda=float(config['lambda']), Mu=float(config['mu']),
                InitNodePositions=initial, InitEdges=np.array([[0, 1], [1, 2]]),
                Do_PCA=False, CenterData=False, PointWeights=weights[:, None],
                FinalEnergy='Base', alpha=0, beta=0, MaxNumberOfIterations=100,
                eps=1e-5, n_cores=1, nReps=1, ProbPoint=1, EmbPointProb=1,
                SampleIC=False, AvoidResampling=True, ComputeMSEP=False,
                verbose=False, ShowTimer=False, DisplayWarnings=False)[0]
        nodes = output['NodePositions']; edges = output['Edges'][0].astype(int)
        nearest = cdist(x, nodes, metric='sqeuclidean').min(axis=1)[:, None]
        energy = ComputePrimitiveGraphElasticEnergy(nodes, output['ElasticMatrix'], nearest, weights[:, None])
        attempts.append(dict(family=family, nodes=nodes, edges=edges,
            training_objective=float(energy[0]), node_mse=float(energy[1]),
            edge_energy=float(energy[2]), star_energy=float(energy[3]), initialization=label,
            warning_count=len(caught), printed_warning_lines=len(captured.getvalue().splitlines())))
    chosen = min(attempts, key=lambda m: m['training_objective'])
    chosen['starts'] = [{k: v for k, v in m.items() if k not in ['nodes', 'edges']} for m in attempts]
    chosen['config'] = dict(config)
    return chosen


def fit_graph(x, blocks, family, config=None, seed=SEED):
    """Public API used unchanged by calibration; no labels/group metadata."""
    x, blocks = np.asarray(x), np.asarray(blocks)
    weights = COMMON.sheet_weights(blocks)
    if family == 'plane':
        mean = weights @ x
        _, _, axes = np.linalg.svd(np.sqrt(weights[:, None]) * (x - mean), full_matrices=False)
        return dict(family='plane', mean=mean, axes=axes[:2], config=None)
    if family in ['curve', 'tree']:
        return _single_graph(x, blocks, family, config)
    assert family == 'two_paths'
    # This control has learned clusters but fixed disconnected-path topology.
    # It is not advertised as unconstrained forest topology discovery.
    clustering = KMeans(n_clusters=2, n_init=10, random_state=seed).fit(x, sample_weight=weights)
    nodes, edges, components = [], [], []
    for group in range(2):
        selected = clustering.labels_ == group
        local = dict(config); local['nodes'] = int(config['nodes']) // 2
        mass = float(weights[selected].sum())
        # Preserve each page's original sheet weight even when one sheet spans
        # both clusters. Scale elastic terms so the global objective is the
        # same weighted node-MSE plus sum of edge/star penalties.
        local['lambda'] /= mass; local['mu'] /= mass
        fit = _single_graph(x[selected], blocks[selected], 'curve', local, weights[selected])
        edges.append(fit['edges'] + sum(len(v) for v in nodes)); nodes.append(fit['nodes'])
        components.append(dict(n_points=int(selected.sum()), n_sheets=len(set(blocks[selected])),
                               point_weight_mass=mass,
                               global_objective_contribution=mass * fit['training_objective'], starts=fit['starts']))
    return dict(family=family, nodes=np.vstack(nodes), edges=np.vstack(edges), config=dict(config), components=components)


def graph_junctions(model):
    if model['family'] == 'plane':
        return []
    degree = np.bincount(model['edges'].ravel(), minlength=len(model['nodes']))
    return np.flatnonzero(degree >= 3).tolist()


def arm_memberships(model, x):
    """Point partitions after deleting each junction, nearest-node assignment."""
    if model['family'] == 'plane':
        return []
    adjacency = [set() for _ in model['nodes']]
    for a, b in model['edges']:
        adjacency[a].add(b); adjacency[b].add(a)
    nearest = cdist(x, model['nodes'], metric='sqeuclidean').argmin(axis=1)
    records = []
    for center in graph_junctions(model):
        components = []
        seen = {center}
        for neighbor in sorted(adjacency[center]):
            if neighbor in seen:
                continue
            component = set(); pending = [neighbor]
            while pending:
                node = pending.pop()
                if node in seen:
                    continue
                seen.add(node); component.add(node); pending.extend(adjacency[node] - seen)
            components.append(component)
        records.append(dict(node=center, degree=len(adjacency[center]),
                            arms=[np.flatnonzero(np.isin(nearest, list(c))) for c in components],
                            junction_points=np.flatnonzero(nearest == center)))
    return records


def _candidate_configs(family):
    return [None] if family == 'plane' else CONFIGS


def nested_cv(d, blocks, families=FAMILIES, seed=SEED, outer_splits=10, inner_splits=3, configs=None):
    """Full training-only geometry, inner selection, and outer evaluation."""
    blocks = np.asarray(blocks)
    grid = CONFIGS if configs is None else configs
    outer = COMMON.make_folds(blocks, outer_splits, seed)
    errors = {family: np.full(len(d), np.nan) for family in families}
    folds = []
    for fold, (train, test) in enumerate(outer):
        started = time.time()
        inner = []
        for itrain, itest in COMMON.make_folds(blocks[train], inner_splits, seed + fold + 1):
            a, b = train[itrain], train[itest]
            geometry = COMMON.fold_geometry(d[np.ix_(a, a)], d[np.ix_(b, a)])
            inner.append((a, b, geometry))
        selected, scores = {}, {}
        for family in families:
            candidate_scores = []
            candidates = [None] if family == 'plane' else grid
            for config in candidates:
                held = np.full(len(d), np.nan)
                for a, b, geometry in inner:
                    fit = fit_graph(geometry['x_train'], blocks[a], family, config, seed + fold)
                    held[b] = project_errors(geometry['x_test'], fit) + geometry['test_offset']
                candidate_scores.append(COMMON.summarize_errors(held[train], blocks[train])['equal_bifolium_mean'])
            chosen = int(np.argmin(candidate_scores))
            selected[family] = candidates[chosen]
            scores[family] = candidate_scores
        geometry = COMMON.fold_geometry(d[np.ix_(train, train)], d[np.ix_(test, train)])
        models = {}
        for family in families:
            fit = fit_graph(geometry['x_train'], blocks[train], family, selected[family], seed + fold)
            errors[family][test] = project_errors(geometry['x_test'], fit) + geometry['test_offset']
            record = dict(selected_config=selected[family], inner_equal_sheet_scores=scores[family],
                outer_equal_sheet_score=COMMON.summarize_errors(errors[family][test], blocks[test])['equal_bifolium_mean'],
                junction_count=len(graph_junctions(fit)),
                branch_arms=[dict(node=r['node'], degree=r['degree'],
                                 arm_training_point_indices=[train[v].tolist() for v in r['arms']],
                                 junction_training_point_indices=train[r['junction_points']].tolist())
                             for r in arm_memberships(fit, geometry['x_train'])])
            models[family] = record
        folds.append(dict(fold=fold, train_indices=train.tolist(), test_indices=test.tolist(),
                          heldout_bifolia=sorted(set(blocks[test])), geometry=geometry['info'],
                          heldout_offset_range=[float(geometry['test_offset'].min()), float(geometry['test_offset'].max())],
                          models=models, seconds=time.time()-started))
        print('outer', fold + 1, 'of', len(outer), 'n', len(d), 'seconds', round(time.time()-started, 2), flush=True)
    summaries = {family: COMMON.summarize_errors(v, blocks) for family, v in errors.items()}
    names = sorted(set(blocks)); rng = np.random.default_rng(seed)
    sample = rng.integers(0, len(names), size=(5000, len(names)))
    contrasts = {}
    for other in families:
        if other == 'curve':
            continue
        a = np.array([summaries['curve']['per_bifolium'][b] for b in names])
        b = np.array([summaries[other]['per_bifolium'][q] for q in names])
        delta = (a - b)[sample].mean(axis=1)
        gain = 1 - b[sample].mean(axis=1) / a[sample].mean(axis=1)
        contrasts[other + '_versus_curve'] = dict(equal_sheet_error_reduction=float((a-b).mean()),
            relative_reduction=float(1 - b.mean()/a.mean()),
            fixed_fold_sheet_bootstrap_95_error_reduction=np.quantile(delta, [.025, .975]).tolist(),
            fixed_fold_sheet_bootstrap_95_relative_reduction=np.quantile(gain, [.025, .975]).tolist(),
            sheets_better=int((b < a).sum()), total_sheets=len(a))
    return dict(summaries=summaries, contrasts=contrasts, folds=folds,
                point_errors={k:v.tolist() for k,v in errors.items()})


def select_full(d, blocks, family, seed=SEED):
    scores = []
    for config in _candidate_configs(family):
        errors = np.empty(len(d))
        for fold, (train, test) in enumerate(COMMON.make_folds(blocks, 3, seed + 500)):
            geometry = COMMON.fold_geometry(d[np.ix_(train, train)], d[np.ix_(test, train)])
            fit = fit_graph(geometry['x_train'], blocks[train], family, config, seed + fold)
            errors[test] = project_errors(geometry['x_test'], fit) + geometry['test_offset']
        scores.append(COMMON.summarize_errors(errors, blocks)['equal_bifolium_mean'])
    config = _candidate_configs(family)[int(np.argmin(scores))]
    x, values = COMMON.embed(d)
    return fit_graph(x, blocks, family, config, seed), x, scores


def _branch_display(model, x, meta, folds):
    result = []
    for branch in arm_memberships(model, x):
        arms = []
        for members in branch['arms']:
            counts = {}
            for i in members:
                group = meta[i]['group']; counts[group] = counts.get(group, 0) + 1
            arms.append(dict(point_indices=members.tolist(), pages=len(members),
                             bifolia=sorted(set(meta[i]['bifolium'] for i in members)), group_counts=counts))
        matches = []
        for fold in folds:
            eligible = set(fold['train_indices']); best = 0.
            for candidate in fold['models']['tree']['branch_arms']:
                if len(candidate['arm_training_point_indices']) != len(arms):
                    continue
                matrix = np.zeros((len(arms), len(arms)))
                for i, arm in enumerate(arms):
                    a = set(arm['point_indices']) & eligible
                    for j, other in enumerate(candidate['arm_training_point_indices']):
                        b = set(other)
                        matrix[i, j] = len(a & b) / max(len(a | b), 1)
                rows, cols = linear_sum_assignment(-matrix)
                best = max(best, float(matrix[rows, cols].mean()))
            matches.append(best)
        result.append(dict(node=branch['node'], degree=branch['degree'], arms=arms,
            junction_point_indices=branch['junction_points'].tolist(),
            outer_fold_best_same_degree_arm_jaccard=matches,
            outer_folds_with_match_at_least_0_7=sum(v >= .7 for v in matches)))
    return result


def main():
    started = time.time(); data = COMMON.load_data()
    results = {}
    for mode in ['combined', 'word', 'line']:
        families = FAMILIES if mode == 'combined' else ['curve', 'tree', 'plane']
        print('START paragraph', mode, flush=True)
        results['paragraph_' + mode] = nested_cv(data['ds'][mode], data['blocks'], families)
    pool_blocks = np.array([m['bifolium'] for m in data['pool_meta']])
    print('START pooled combined', flush=True)
    results['pooled_combined'] = nested_cv(data['pool_ds']['combined'], pool_blocks, ['curve', 'tree', 'plane'])
    display_models = {}
    for family in FAMILIES:
        fit, x, scores = select_full(data['ds']['combined'], data['blocks'], family)
        entry = dict(family=family, config=fit['config'], selection_3fold_scores=scores,
                     in_sample_errors=project_errors(x, fit).tolist())
        if family == 'plane':
            entry.update(mean_first4=fit['mean'][:4].tolist(), axes_first4=fit['axes'][:, :4].tolist())
        else:
            entry.update(nodes_first4=fit['nodes'][:, :4].tolist(), edges=fit['edges'].tolist(),
                         junctions=graph_junctions(fit), starts=fit.get('starts'),
                         branch_support=_branch_display(fit, x, data['meta'], results['paragraph_combined']['folds']) if family == 'tree' else [])
        display_models[family] = entry
    dependencies = {p: importlib.metadata.version(p) for p in ['elpigraph-python', 'numba', 'llvmlite', 'igraph', 'networkx', 'shapely', 'numpy', 'scipy', 'scikit-learn']}
    sources = list(data['sources']) + [HERE / 'models.py', HERE / 'common.py']
    package_sources = [HERE / '_deps/elpigraph' / p for p in ['_topologies.py', '_BaseElPiWrapper.py', 'src/BaseElPi.py', 'src/grammar_operations.py', 'src/core.py', 'src/distutils.py']]
    result = dict(question='Does a label-free learned principal tree predict unseen whole-bifolium profiles better than a curve?',
        implementation='Unmodified elpigraph-python0.3.4, official curve/tree grammar wrappers',
        package_url='https://github.com/j-bac/elpigraph-python', dependencies=dependencies,
        design=dict(seed=SEED, outer_folds=10, inner_folds=3, configurations=CONFIGS,
                    training_weights='Equal physical bifolium, uniform pages within bifolium',
                    starts=['weighted_PCA1', 'unlabeled_PCA_bend'], max_EM_iterations=100, eps=1e-5,
                    final_energy='Base; alpha=beta=0; exact weighted node-MSE + edge + star energy',
                    selection_score='Equal-bifolium heldout squared segment error plus Gower orthogonal signed-radius offset'),
        analyses=results, full_models=display_models, elapsed_seconds=time.time()-started,
        limits=[
            'All folds hold out entire physical bifolia; both inner and outer geometry fit training profiles only.',
            'ElPiGraph fits weighted distances to nodes plus elastic terms; validation uses distances to segments for all graph families.',
            'The tree starts from an unlabeled path and learns topology by node/edge graph grammars. Golf groups, section, hand, folio order and published curve positions never select graph connections or initial nodes.',
            'Two starts and a finite complexity/penalty grid do not guarantee a global optimum or exhaust topology possibilities.',
            'The disconnected control is two learned weighted-k-means clusters with independent paths, matched in total node count; it is not unrestricted forest learning.',
            'Positive training MDS discards negative eigenvectors. The test correction evaluates signed source-distance reconstruction in its positive affine span and can yield negative point scores for indefinite dissimilarities.',
            'The correction is common to candidates within a fold and prevents silently discarding unseen test variance.',
            'Fixed-fold paired block bootstrap describes variation across these 50 sheets, not uncertainty from refitting folds or an independent replication.',
            'Reported branch-arm matching across outer refits is descriptive stability, not a posterior probability or calibrated branch confidence.',
            'Word and line representations reuse the same manuscript and are sensitivity checks, not independent samples.',
            'No chronological arrow, missing-text reconstruction, or historical branching mechanism is inferred from graph geometry alone.'
        ], source_sha256=COMMON.hash_sources(sources),
        package_source_sha256={str(p.relative_to(HERE)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in package_sources})
    (HERE / 'models.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    display = dict(meta=data['meta'], coordinates_first4=x[:, :4].tolist(), models=display_models)
    (HERE / 'model_display.json').write_text(json.dumps(display, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    assert (HERE / 'models.json').stat().st_size + (HERE / 'model_display.json').stat().st_size < 5_000_000
    print('COMPLETE',time.time()-started,flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
