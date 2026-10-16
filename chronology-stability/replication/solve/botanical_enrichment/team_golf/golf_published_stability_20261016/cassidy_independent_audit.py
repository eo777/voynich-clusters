"""Bounded independent audit of Cassidy continuation; no growth or model archive."""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cassidy_h_continuation as impl


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def nearest(x, y):
    return ((x[:, None, :]-y[None, :, :])**2).sum(axis=2).argmin(axis=1)


def elastic(n):
    edge = np.zeros((n-1, n))
    for j in range(n-1):
        edge[j, j:j+2] = [-1, 1]
    bend = np.zeros((n-2, n))
    for j in range(n-2):
        bend[j, j:j+3] = [-.5, 1, -.5]
    return .01*edge.T@edge+.1*bend.T@bend


def solve_branch(x, y):
    lap = elastic(len(y))
    for iteration in range(1000):
        assigned = nearest(x, y)
        masses = np.array([(assigned == j).mean() for j in range(len(y))])
        sums = np.array([x[assigned == j].sum(axis=0)/len(x) for j in range(len(y))])
        fitted = np.linalg.solve(lap+np.diag(masses), sums)
        if np.array_equal(assigned, nearest(x, fitted)):
            return fitted
        y = fitted
    raise AssertionError('Independent corrector did not converge')


def coefficient(x, y, hi):
    assigned = nearest(x, y)
    matrix = elastic(len(y))+np.diag([(assigned == j).mean() for j in range(len(y))])
    rhs = np.zeros(len(y)); rhs[assigned[hi]] = 1/len(x)
    return np.linalg.solve(matrix, rhs), assigned


def event(x, y, dx, dy, assigned):
    # Scalar polynomial roots provide a structurally separate check of the
    # production vectorized quadratic formula, including both candidate roots.
    winner = (float('inf'), -1, -1)
    for i in range(len(x)):
        a = assigned[i]
        r0 = x[i]-y[a]; v0 = dx[i]-dy[a]
        for j in range(len(y)):
            if j == a:
                continue
            r = x[i]-y[j]; v = dx[i]-dy[j]
            c = float(r@r-r0@r0)
            b = float(2*(r@v-r0@v0))
            aa = float(v@v-v0@v0)
            coefficients = [aa, b, c] if abs(aa) >= 1e-14 else [b, c]
            if len(coefficients) == 2 and abs(b) < 1e-14:
                continue
            for root in np.roots(coefficients):
                if abs(root.imag) > 1e-10:
                    continue
                t = float(root.real)
                if -1e-10 <= t <= 1+1e-10 and b+2*aa*t < -1e-12:
                    candidate = (max(t, 0.), i, j)
                    if candidate < winner:
                        winner = candidate
    return winner


def follow(base, y, hi, plane, waypoints, nudge=1e-6):
    events = 0
    for start, target in zip(waypoints[:-1], waypoints[1:]):
        x = base.copy(); x[hi, plane] = start
        for iteration in range(100):
            delta = np.asarray(target)-x[hi, plane]
            if np.linalg.norm(delta) < 1e-13:
                break
            coeff, assigned = coefficient(x, y, hi)
            dx = np.zeros_like(x); dx[hi, plane] = delta
            dy = coeff[:, None]*dx[hi]
            fraction, _, _ = event(x, y, dx, dy, assigned)
            boundary = fraction < 1-1e-10
            step = min(1., fraction+nudge) if boundary else 1.
            prediction = y+step*dy
            x += step*dx
            y = solve_branch(x, prediction)
            if not boundary:
                assert np.array_equal(assigned, nearest(x, y))
                assert np.max(abs(y-prediction)) < 1e-9
            events += boundary
        else:
            raise AssertionError('Too many independent events')
    return y, int(events)


def readout(x, y, meta):
    reconstruction = ((x[:, None, :]-y[None, :, :])**2).sum(axis=2).min(axis=1).mean()
    objective = reconstruction+np.sum(y*(elastic(len(y))@y))
    lengths = np.linalg.norm(np.diff(y, axis=0), axis=1)
    arc = np.r_[0., np.cumsum(lengths)]
    positions = []
    for point in x:
        options = []
        for j in range(len(y)-1):
            v = y[j+1]-y[j]
            f = np.clip((point-y[j])@v/max(v@v, 1e-15), 0, 1)
            options.append((np.sum((point-y[j]-f*v)**2), j, f))
        _, j, f = min(options)
        positions.append((arc[j]+f*lengths[j])/arc[-1])
    groups = np.array([m['group'] for m in meta]); positions = np.asarray(positions)
    contrasts = [float(np.median(positions[groups == a])-np.median(positions[groups == b]))
                 for a, b in [('1a', '1b'), ('3a', '3b'), ('4a', '4c')]]
    signs = [(-1 if v < -1e-8 else 1 if v > 1e-8 else 0) for v in contrasts]
    return float(objective), signs


def main():
    started = time.time()
    result_path = HERE/'cassidy_h_results.json'
    result = json.loads(result_path.read_text())
    assert result['status'] == 'complete'
    for relative, expected in result['source_sha256'].items():
        assert digest(HERE.parent/relative) == expected
    old = json.loads((impl.OLD/'expanded_influence_grid.json').read_text())
    basins = json.loads((HERE.parent/'golf_stability_20261016/basins.json').read_text())
    native = next(s for s in basins['settings'] if s['cohort'] == 'pool50' and s['lambda'] == .01)
    extra = next(s for s in old['seed_sets'] if s['penalty'] == .01)
    seeds = [dict(name=s['start'], nodes=np.asarray(s['nodes'])) for s in native['baseline_runs']]
    seeds += [dict(name=s['name'], nodes=np.asarray(s['nodes'])) for s in extra['added_full_models']]
    data = impl.grid.COMMON.load_data()
    x, _ = impl.grid.COMMON.embed(data['pool_ds']['combined'])
    hi = result['h_index']; meta = result['meta']
    out = dict(status='running', derivative_checks=[], event_checks=[], path_checks=[],
               source_sha256={'cassidy_h_results.json': digest(result_path),
                              'cassidy_h_continuation.py': digest(HERE/'cassidy_h_continuation.py'),
                              'cassidy_independent_audit.py': digest(__file__)})
    with threadpool_limits(limits=1):
        for seed in seeds:
            y = seed['nodes']; coeff, assigned = coefficient(x, y, hi)
            production, _, _ = impl.derivative(x, y, hi)
            assert np.max(abs(coeff-production)) < 1e-13
            for dimension in [0, 1, 2]:
                epsilon = 1e-6
                plus = x.copy(); minus = x.copy()
                plus[hi, dimension] += epsilon; minus[hi, dimension] -= epsilon
                yp, ym = solve_branch(plus, y), solve_branch(minus, y)
                assert np.array_equal(nearest(plus, yp), assigned)
                assert np.array_equal(nearest(minus, ym), assigned)
                expected = np.zeros_like(y); expected[:, dimension] = coeff
                error = float(np.linalg.norm((yp-ym)/(2*epsilon)-expected)/np.linalg.norm(expected))
                assert error < 1e-7
                out['derivative_checks'].append(dict(seed=seed['name'], dimension=dimension, relative_error=error))
        for row in result['settings']:
            plane = row['plane']; start = x[hi, plane].copy()
            # Fixed baseline/edge/corner coverage plus largest-event and
            # path-preference-disagreement cells selected from completed data.
            selected = {row['baseline_cell_index'], 0, 24, 600, 624}
            selected.add(max(range(625), key=lambda k: sum(row['cells'][k]['paths'][0]['assignment_events'])))
            for k, cell in enumerate(row['cells']):
                signs = [np.sign(p['herbal_preference']['gap_percent']) if p['herbal_preference']['gap_percent'] is not None else None for p in cell['paths']]
                if len(set(signs)) > 1:
                    selected.add(k); break
            for index in sorted(selected):
                target = np.array([row['xs'][index%25], row['ys'][index//25]])
                changed = x.copy(); changed[hi, plane] = target
                for seed_index in [0, 9]:
                    seed = seeds[seed_index]; y = seed['nodes']
                    coeff, assigned = coefficient(x, y, hi)
                    dx = np.zeros_like(x); dx[hi, plane] = target-start
                    dy = coeff[:, None]*dx[hi]
                    independent = event(x, y, dx, dy, assigned)
                    production = impl.first_assignment_event(x, y, dx, dy, assigned)
                    finite = np.isfinite(independent[0])
                    assert finite == np.isfinite(production[0])
                    root_error = abs(independent[0]-production[0]) if finite else 0.
                    assert root_error < 1e-8
                    out['event_checks'].append(dict(plane=plane, cell=index, seed=seed['name'],
                        has_event=bool(finite), root_error=float(root_error)))
                    for route_index, route in enumerate(['direct', 'x_then_y', 'y_then_x']):
                        middle = np.array([target[0], start[1]]) if route_index == 1 else np.array([start[0], target[1]])
                        waypoints = [start, target] if route_index == 0 else [start, middle, target]
                        fitted, events = follow(x, y, hi, plane, waypoints)
                        objective, signs = readout(changed, fitted, meta)
                        reference = row['cells'][index]['paths'][route_index]
                        error = abs(objective-reference['objectives'][seed_index])
                        assert error < 1e-11 and signs == reference['signs'][seed_index]
                        assert events == reference['assignment_events'][seed_index]
                        returned_error = None
                        if route_index == 0:
                            returned, _ = follow(x, fitted, hi, plane, [target, start])
                            _, returned_signs = readout(x, returned, meta)
                            _, original_signs = readout(x, y, meta)
                            assert returned_signs == original_signs
                            returned_error = float(np.sqrt(np.mean((returned-y)**2)))
                        out['path_checks'].append(dict(plane=plane, cell=index, seed=seed['name'], route=route,
                            objective_error=float(error), signs=signs, assignment_events=events,
                            reverse_node_rms=returned_error))
    assert digest(result_path) == out['source_sha256']['cassidy_h_results.json']
    out.update(status='passed', elapsed_seconds=time.time()-started,
        scope='30 independent finite differences; scalar-polynomial event roots; independent fixed-chain corrector and segment readout; selected direct/axis/reverse paths only; no new growth.')
    (HERE/'cassidy_independent_audit.json').write_text(json.dumps(out, separators=(',', ':'), allow_nan=False)+'\n')
    print(json.dumps(dict(status=out['status'], derivatives=len(out['derivative_checks']),
        events=len(out['event_checks']), paths=len(out['path_checks']), elapsed=out['elapsed_seconds'])))


if __name__ == '__main__':
    main()
