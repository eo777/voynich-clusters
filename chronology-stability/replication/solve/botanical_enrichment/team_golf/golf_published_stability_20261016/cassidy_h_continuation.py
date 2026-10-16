"""Cassidy Eq.8 predictor with explicit nearest-node assignment events.

The paper's smooth IFT calculation is used only inside assignment regions.
Crossing a region boundary invokes a corrector and a new derivative, not an
assertion that the globally nonsmooth principal-curve objective is smooth.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'golf_influence_20261016'
sys.path.insert(0, str(OLD))
import coordinate_slices as grid

LAM = .01
MU = .1
NODES = 28
NUDGE = 1e-6
LAP, EDGE, STAR = grid.matrices(NODES, LAM, MU)
COUNTERS = dict(corrector_calls=0, corrector_iterations=0, derivative_solves=0,
                smooth_steps=0, assignment_events=0, paths=0)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(name, value):
    (HERE/name).write_text(json.dumps(value, separators=(',', ':'), allow_nan=False)+'\n', encoding='utf-8')


def assignments(x, y):
    return cdist(x, y, 'sqeuclidean').argmin(axis=1)


def gradient(x, y, assignment):
    mass = np.bincount(assignment, minlength=NODES)/len(x)
    sums = np.zeros_like(y)
    np.add.at(sums, assignment, x/len(x))
    return 2*((np.diag(mass)+LAP)@y-sums)


def correct(x, y):
    COUNTERS['corrector_calls'] += 1
    ans, energy, iterations, increase = grid.refine(x, y, LAM, tol=1e-11)
    COUNTERS['corrector_iterations'] += iterations
    assert np.max(np.abs(gradient(x, ans, assignments(x, ans)))) < 1e-10
    return ans


def derivative(x, y, hi):
    """Solve H_theta_theta J = -H_theta_phi, Cassidy Eq.8/section2.3.

    The 49 identical coordinate blocks need only one 28-dimensional solve.
    J for coordinate d is coeff in column d and zero in other columns.
    """
    assignment = assignments(x, y)
    mass = np.bincount(assignment, minlength=NODES)/len(x)
    hessian = 2*(np.diag(mass)+LAP)
    mixed = np.zeros(NODES)
    mixed[assignment[hi]] = -2/len(x)
    coeff = np.linalg.solve(hessian, -mixed)
    eigen = float(np.linalg.eigvalsh(hessian)[0])
    assert eigen > 0
    COUNTERS['derivative_solves'] += 1
    return coeff, assignment, eigen


def first_assignment_event(x, y, dx, dy, assignment):
    """Earliest distance-equality crossing on the exact affine prediction.

    For every observation/competing node, g(t)=d_competing²-d_assigned²
    is a quadratic. Only positive-to-negative crossings invalidate the region.
    An equality without crossing (tangency) is not a required event here.
    """
    r = x[:, None, :]-y[None, :, :]
    v = dx[:, None, :]-dy[None, :, :]
    ar = r[np.arange(len(x)), assignment]
    av = v[np.arange(len(x)), assignment]
    q0 = np.sum(r*r, axis=2)-np.sum(ar*ar, axis=1)[:, None]
    q1 = 2*(np.sum(r*v, axis=2)-np.sum(ar*av, axis=1)[:, None])
    q2 = np.sum(v*v, axis=2)-np.sum(av*av, axis=1)[:, None]
    q0[np.arange(len(x)), assignment] = np.inf
    roots = np.full(q0.shape, np.inf)
    linear = (np.abs(q2) < 1e-14) & (q1 < -1e-14) & np.isfinite(q0)
    roots[linear] = -q0[linear]/q1[linear]
    nonlinear = (np.abs(q2) >= 1e-14) & np.isfinite(q0)
    disc = q1*q1-4*q2*np.where(np.isfinite(q0), q0, 0)
    valid = nonlinear & (disc >= 0)
    for sign in [-1, 1]:
        candidate = np.full(q0.shape, np.inf)
        candidate[valid] = (-q1[valid]+sign*np.sqrt(disc[valid]))/(2*q2[valid])
        downward = valid & (candidate >= -1e-10) & (candidate <= 1+1e-10)
        downward &= q1+2*q2*np.where(np.isfinite(candidate), candidate, 0) < -1e-12
        roots = np.minimum(roots, np.where(downward, np.maximum(candidate, 0), np.inf))
    roots[(roots < -1e-10) | (roots > 1+1e-10)] = np.inf
    index = np.unravel_index(np.argmin(roots), roots.shape)
    return float(roots[index]), int(index[0]), int(index[1])


def empty_stats():
    return dict(steps=0, events=0, assignment_changes=0, max_correction_rms=0.,
                max_smooth_correction_rms=0., min_hessian_eigenvalue=1.,
                max_event_distance_residual=0., event_log=[])


def add_stats(a, b, retain=False):
    for key in ['steps', 'events', 'assignment_changes']:
        a[key] += b[key]
    for key in ['max_correction_rms', 'max_smooth_correction_rms', 'max_event_distance_residual']:
        a[key] = max(a[key], b[key])
    a['min_hessian_eigenvalue'] = min(a['min_hessian_eigenvalue'], b['min_hessian_eigenvalue'])
    if retain:
        a['event_log'].extend(b['event_log'])


def continue_segment(base_x, y, hi, plane, start, target, nudge=NUDGE, retain=False):
    x = base_x.copy()
    current = np.asarray(start, dtype=float).copy()
    target = np.asarray(target, dtype=float)
    x[hi, plane] = current
    stats = empty_stats()
    for _ in range(100):
        delta = target-current
        if np.linalg.norm(delta) < 1e-13:
            return y, stats
        coeff, assignment, mineig = derivative(x, y, hi)
        dx = np.zeros_like(x); dx[hi, plane] = delta
        dy = coeff[:, None]*dx[hi][None, :]
        fraction, observation, competitor = first_assignment_event(x, y, dx, dy, assignment)
        is_event = fraction < 1-1e-10
        step = min(1., fraction+nudge) if is_event else 1.
        predicted = y+step*dy
        updated_x = x+step*dx
        fitted = correct(updated_x, predicted)
        changes = int(np.sum(assignments(updated_x, fitted) != assignment))
        error = float(np.sqrt(np.mean((fitted-predicted)**2)))
        stats['steps'] += 1
        stats['assignment_changes'] += changes
        stats['min_hessian_eigenvalue'] = min(stats['min_hessian_eigenvalue'], mineig)
        stats['max_correction_rms'] = max(stats['max_correction_rms'], error)
        if is_event:
            COUNTERS['assignment_events'] += 1
            stats['events'] += 1
            event_x = x+fraction*dx
            event_y = y+fraction*dy
            diff = abs(np.sum((event_x[observation]-event_y[assignment[observation]])**2)
                       -np.sum((event_x[observation]-event_y[competitor])**2))
            stats['max_event_distance_residual'] = max(stats['max_event_distance_residual'], float(diff))
            if retain:
                stats['event_log'].append(dict(coordinates=(current+fraction*delta).tolist(),
                    observation=int(observation), previous_node=int(assignment[observation]),
                    competing_node=int(competitor), equality_residual=float(diff),
                    correction_rms=error, assignment_changes=changes))
        else:
            COUNTERS['smooth_steps'] += 1
            stats['max_smooth_correction_rms'] = max(stats['max_smooth_correction_rms'], error)
            assert changes == 0, 'Undetected nearest-node assignment crossing'
            assert error < 1e-9, 'Affine predictor should be exact in a smooth region'
        y, x = fitted, updated_x
        current = x[hi, plane].copy()
    raise AssertionError('More than 100 events on one segment')


def path(base_x, y, hi, plane, waypoints, nudge=NUDGE, retain=False):
    COUNTERS['paths'] += 1
    stats = empty_stats()
    for start, target in zip(waypoints[:-1], waypoints[1:]):
        y, sub = continue_segment(base_x, y, hi, plane, start, target, nudge, retain)
        add_stats(stats, sub, retain)
    return y, stats


def summary(x, y, meta):
    objective, _, _, _, _ = grid.components(x, y, EDGE, STAR, LAM)
    position, _ = grid.project(x, y)
    return float(objective), grid.signs(position, meta)


def finite_difference_checks(x, seeds, hi):
    rows = []
    for seed in seeds:
        y = seed['nodes']
        coeff, assigned, eigen = derivative(x, y, hi)
        for dimension in [0, 1, 2]:
            expected = np.zeros_like(y); expected[:, dimension] = coeff
            for epsilon in [1e-5, 1e-6]:
                plus = x.copy(); minus = x.copy()
                plus[hi, dimension] += epsilon; minus[hi, dimension] -= epsilon
                yp, ym = correct(plus, y), correct(minus, y)
                assert np.array_equal(assignments(plus, yp), assigned)
                assert np.array_equal(assignments(minus, ym), assigned)
                numeric = (yp-ym)/(2*epsilon)
                relative = float(np.linalg.norm(numeric-expected)/np.linalg.norm(expected))
                assert relative < 1e-7
                mixed_numeric = (gradient(plus, y, assigned)-gradient(minus, y, assigned))/(2*epsilon)
                mixed_expected = np.zeros_like(y); mixed_expected[assigned[hi], dimension] = -2/len(x)
                mixed_error = float(np.max(np.abs(mixed_numeric-mixed_expected)))
                assert mixed_error < 1e-9
                # Central gradient differences independently verify a Hessian action.
                hessian_numeric = (gradient(x, y+epsilon*expected, assigned)
                                   -gradient(x, y-epsilon*expected, assigned))/(2*epsilon)
                hessian_error = float(np.max(np.abs(hessian_numeric+mixed_expected)))
                assert hessian_error < 1e-9
                rows.append(dict(seed=seed['name'], dimension=dimension, epsilon=epsilon,
                    derivative_relative_error=relative, mixed_derivative_max_error=mixed_error,
                    hessian_action_max_error=hessian_error, minimum_hessian_eigenvalue=eigen))
    return rows


def preference(objectives, signs, pair=1):
    minima = []
    for direction in [-1, 1]:
        eligible = [e for e, s in zip(objectives, signs) if s[0] == -1 and s[pair] == direction]
        minima.append(min(eligible) if eligible else None)
    gap = 100*(minima[1]-minima[0])/min(minima) if all(v is not None for v in minima) else None
    return dict(before_objective=minima[0], after_objective=minima[1], gap_percent=gap)


def main():
    started = time.time()
    data = grid.COMMON.load_data()
    x, eigenvalues = grid.COMMON.embed(data['pool_ds']['combined'])
    old = json.loads((OLD/'expanded_influence_grid.json').read_text())
    basins_path = HERE.parent/'golf_stability_20261016/basins.json'
    basins = json.loads(basins_path.read_text())
    base = next(s for s in basins['settings'] if s['cohort']=='pool50' and s['lambda']==LAM)
    enriched = next(s for s in old['seed_sets'] if s['penalty']==LAM)
    seeds = [dict(name=r['start'], nodes=np.asarray(r['nodes'])) for r in base['baseline_runs']]
    seeds += [dict(name=r['name'], nodes=np.asarray(r['nodes'])) for r in enriched['added_full_models']]
    assert len(seeds) == 10 and [s['name'] for s in seeds] == enriched['initializations']
    hi = old['h_index']; meta = old['meta']
    paths = [Path(__file__), HERE/'CASSIDY_DESIGN.md', OLD/'coordinate_slices.py',
             OLD/'expanded_influence_grid.json', basins_path]
    hashes = {str(p.relative_to(HERE.parent)): sha(p) for p in paths}
    out = dict(status='running', source_sha256=hashes, citation=dict(
        title='A Continuation Technique for Maximum Likelihood Estimators in Biological Models',
        author='Tyler Cassidy', year=2023, doi='10.1007/s11538-023-01200-0',
        implementation='Equations7-8 analytic implicit derivative, section2.3 linear solve; event-aware correction is the declared piecewise-smooth adaptation.'),
        design=dict(penalty=LAM, mu=MU, nodes=NODES, dimensions=x.shape[1],
            initializations=[s['name'] for s in seeds], event_crossing_fraction=NUDGE,
            paths=['direct', 'x_then_y', 'y_then_x'], reverse_check='direct endpoint back to actual H',
            historical_probability=False, scope='Fixed hidden geometry; final1a-first filtered; not arbitrary missing-count completion.'),
        meta=meta, h_index=hi, reference_coords3=x[:, :3].tolist(), palette=old['palette'],
        positive_variance_shares=old['positive_variance_shares'], settings=[], path_refinement_checks=[])
    aggregate = empty_stats()
    with threadpool_limits(limits=1):
        out['finite_difference_checks'] = finite_difference_checks(x, seeds, hi)
        for oldrow in [s for s in old['settings'] if s['penalty']==LAM]:
            plane = oldrow['plane']; start = x[hi, plane].copy()
            row = {k: oldrow[k] for k in ['plane', 'xs', 'ys', 'baseline_cell_index']}
            row.update(cells=[], selected_event_traces=[])
            for iy, yy in enumerate(row['ys']):
                for ix, xx in enumerate(row['xs']):
                    target = np.array([xx, yy]); changed = x.copy(); changed[hi, plane] = target
                    cell = dict(paths=[], maximum_return_node_rms=0., return_order_changes=0,
                                maximum_path_node_rms=0., path_order_changes=0)
                    pathmodels = []
                    for route in ['direct', 'x_then_y', 'y_then_x']:
                        waypoints = [start, target] if route=='direct' else [start,
                            np.array([xx, start[1]]) if route=='x_then_y' else np.array([start[0], yy]), target]
                        energies=[]; relations=[]; eventcounts=[]; corrections=[]; models=[]
                        for si, seed in enumerate(seeds):
                            retain = ix in [0, 12, 24] and iy in [0, 12, 24] and si in [0, 9] and route=='direct'
                            fitted, stats = path(x, seed['nodes'], hi, plane, waypoints, retain=retain)
                            add_stats(aggregate, stats)
                            energy, relation = summary(changed, fitted, meta)
                            energies.append(energy); relations.append(relation)
                            eventcounts.append(stats['events']); corrections.append(stats['max_correction_rms']); models.append(fitted)
                            if retain:
                                row['selected_event_traces'].append(dict(cell_index=iy*25+ix, seed=seed['name'], **stats))
                                half, _ = path(x, seed['nodes'], hi, plane, waypoints, nudge=NUDGE/2)
                                half_energy, half_relation = summary(changed, half, meta)
                                out['path_refinement_checks'].append(dict(plane=plane, cell_index=iy*25+ix, seed=seed['name'],
                                    nudge_half_endpoint_rms=float(np.sqrt(np.mean((fitted-half)**2))),
                                    objective_difference=float(half_energy-energy), same_relation=half_relation==relation))
                            if route=='direct':
                                returned, reverse_stats = path(x, fitted, hi, plane, [target, start])
                                add_stats(aggregate, reverse_stats)
                                rms = float(np.sqrt(np.mean((returned-seed['nodes'])**2)))
                                cell['maximum_return_node_rms'] = max(cell['maximum_return_node_rms'], rms)
                                cell['return_order_changes'] += summary(x, returned, meta)[1] != summary(x, seed['nodes'], meta)[1]
                            else:
                                rms = float(np.sqrt(np.mean((fitted-pathmodels[0][si])**2)))
                                cell['maximum_path_node_rms'] = max(cell['maximum_path_node_rms'], rms)
                                cell['path_order_changes'] += relation != cell['paths'][0]['signs'][si]
                        cell['paths'].append(dict(name=route, objectives=energies, signs=relations,
                            assignment_events=eventcounts, maximum_correction_rms=corrections,
                            herbal_preference=preference(energies, relations), terminal_preference=preference(energies, relations, 2)))
                        pathmodels.append(models)
                    allenergies = [e for p in cell['paths'] for e in p['objectives']]
                    allsigns = [s for p in cell['paths'] for s in p['signs']]
                    cell['combined_herbal_preference'] = preference(allenergies, allsigns)
                    cell['combined_terminal_preference'] = preference(allenergies, allsigns, 2)
                    original = oldrow['cells'][iy*25+ix]
                    cell['original_herbal_preference'] = preference(original['objectives'], original['signs'])
                    row['cells'].append(cell)
                print('CASSIDY', plane, 'row', iy+1, 'events', COUNTERS['assignment_events'],
                      'elapsed', round(time.time()-started, 1), flush=True)
            out['settings'].append(row)
            save('cassidy_h_checkpoint.json', out)
    aggregate.pop('event_log')
    out.update(status='complete', counts=COUNTERS, aggregate=aggregate, elapsed_seconds=time.time()-started)
    assert all(sha(HERE.parent/p)==h for p,h in hashes.items())
    save('cassidy_h_results.json', out)
    print('COMPLETE', COUNTERS, round(out['elapsed_seconds'], 1), flush=True)


if __name__ == '__main__':
    main()
