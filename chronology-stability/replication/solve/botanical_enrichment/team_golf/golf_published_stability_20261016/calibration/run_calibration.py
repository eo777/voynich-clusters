"""Conditional manuscript calibration of the verified, source-defined EGRET score.

See DESIGN.md before interpreting these simulation references as evidence.
No synthetic feature matrices are retained. Original R outputs are immutable.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import time

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
SEED = 20261016
CONFIG = dict(seed=SEED, bandwidths=[.08, .12, .20], models=['sheet_residual', 'dirichlet_multinomial'],
              fit_reps=400, calibration_reps=999, evaluation_reps=500,
              power_calibration_reps=399, power_evaluation_reps=300,
              gap_centers=[.35, .65, .85], gap_widths=[.10, .20, .30])
spec = importlib.util.spec_from_file_location('egret_audit', STUDY/'scegret_independent_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def write(name, value):
    (HERE/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rng(*key):
    return np.random.default_rng(np.random.SeedSequence([SEED]+list(key)))


def wilson(k, n):
    if not n:
        return [None, None]
    z = 1.959963984540054
    p = k/n
    den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(center-half), float(center+half)]


def rate(values):
    n = len(values)
    k = int(np.count_nonzero(values))
    return dict(events=k, trials=n, rate=k/n if n else None, interval95=wilson(k, n))


def tail(samples, observed):
    return float((1+np.count_nonzero(samples >= observed))/(len(samples)+1))


class Scorer:
    def __init__(self, positions, steps, h_indices=(), scales=('per_million', 'raw')):
        self.t = np.asarray(positions)
        self.views = []
        self.labels = []
        for step in steps:
            assignment, _, occupancy = audit.bins(self.t, step)
            kept = np.flatnonzero(occupancy > 3)+1
            weights = np.column_stack([(assignment == b)/occupancy[b-1] for b in kept]) if len(kept) else np.empty((len(self.t), 0))
            for scale in scales:
                self.views.append((step, scale, kept, weights, assignment, occupancy))
                for a, b in zip(kept[:-1], kept[1:]):
                    hb = set(assignment[list(h_indices)]) if len(h_indices) else set()
                    self.labels.append(dict(steps=step, scale=scale, left=int(a), right=int(b),
                        left_n=int(occupancy[a-1]), right_n=int(occupancy[b-1]),
                        touches_H=bool(a in hb or b in hb)))

    def score(self, counts):
        detected = np.count_nonzero(counts, axis=0) > len(counts)*.1
        if not np.any(detected):
            raise ValueError('No detected features')
        raw = counts[:, detected]
        cpm = raw/counts.sum(axis=1)[:, None]*1e6
        out = []
        for _, scale, kept, weights, _, _ in self.views:
            matrix = cpm if scale == 'per_million' else raw
            if len(kept) < 2:
                continue
            means = np.log1p(weights.T @ matrix)
            denominator = np.log1p(matrix.mean(axis=0))
            out.extend(np.mean(np.abs(np.diff(means, axis=0))/denominator, axis=1))
        return np.asarray(out)


class Generator:
    def __init__(self, counts, positions, sheets, bandwidth):
        self.x = counts
        self.n = counts.sum(axis=1).astype(int)
        self.t = np.asarray(positions)
        self.p = counts/self.n[:, None]
        self.sheets = np.asarray(sheets)
        self.groups = [np.flatnonzero(self.sheets == s) for s in sorted(set(sheets))]
        self.groups = [g[np.argsort(self.t[g], kind='stable')] for g in self.groups]
        self.h = bandwidth
        w = np.exp(-.5*((self.t[:, None]-self.t[None, :])/bandwidth)**2)
        self.mu = w @ self.p/w.sum(axis=1)[:, None]
        loo = w*(self.sheets[:, None] != self.sheets[None, :])
        prediction = loo @ self.p/loo.sum(axis=1)[:, None]
        residual_sq = np.sum((self.p-prediction)**2, axis=1)
        base = 1-np.sum(prediction**2, axis=1)
        # E ||p-mu||^2 = (1-||mu||^2) [1/N + rho(1-1/N)].
        self.rho = float(np.clip(np.sum(residual_sq-base/self.n)/np.sum(base*(1-1/self.n)), 1e-7, .95))
        self.alpha = 1/self.rho-1
        self.residual = np.sqrt(self.p)-np.sqrt(prediction)
        self.residual -= self.residual.mean(axis=0)
        self.root_mu = np.sqrt(self.mu)
        self.donor_maps = []
        for target in self.groups:
            self.donor_maps.append([donor[np.minimum(((np.arange(len(target))+.5)/len(target)*len(donor)).astype(int), len(donor)-1)] for donor in self.groups])
        self.observed_diag = self.diagnostics(counts)

    def diagnostics(self, counts):
        p = counts/counts.sum(axis=1)[:, None]
        return dict(detected_features=int(np.sum(np.count_nonzero(counts, axis=0) > len(counts)*.1)),
                    zero_fraction=float(np.mean(counts == 0)),
                    mean_squared_residual=float(np.mean(np.sum((p-self.mu)**2, axis=1))))

    def draw(self, model, random):
        if model == 'dirichlet_multinomial':
            pp = random.gamma(np.maximum(self.mu*self.alpha, 1e-12), 1.)
            pp /= pp.sum(axis=1)[:, None]
            return np.array([random.multinomial(int(n), p) for n, p in zip(self.n, pp)])
        donors = np.empty(len(self.n), dtype=int)
        for target, maps, j in zip(self.groups, self.donor_maps, random.integers(0, len(self.groups), len(self.groups))):
            donors[target] = maps[j]
        scale = np.sqrt((self.rho+1/self.n)/(self.rho+1/self.n[donors]))
        root = np.maximum(self.root_mu+self.residual[donors]*scale[:, None], 0)
        pp = root*root
        pp /= pp.sum(axis=1)[:, None]
        expected = pp*self.n[:, None]
        result = np.floor(expected).astype(int)
        short = self.n-result.sum(axis=1)
        fractions = expected-result
        for i, k in enumerate(short):
            if k:
                chosen = np.argpartition(fractions[i], -k)[-k:]
                result[i, chosen] += 1
        assert np.array_equal(result.sum(axis=1), self.n)
        return result


def load():
    data = json.loads((STUDY/'scegret_inputs.json').read_text())
    out = []
    for ci, c in enumerate(data['cohorts']):
        t = np.asarray(next(r['positions'] for r in c['candidates'] if r['id'] == 'original_golf_spline'))
        t = (t-t.min())/(t.max()-t.min())
        h_ids = [i for i, m in enumerate(c['meta']) if m['id'] in ['f58r', 'f58v', 'H:58-65']]
        for channel, mode in enumerate(['word', 'line']):
            ids, features, raw = audit.read_matrix(STUDY/c['matrices'][mode+'_raw'])
            x = raw.T
            active = x.sum(axis=0) > 0
            out.append(dict(ci=ci, channel=channel, cohort=c['cohort'], mode=mode, ids=ids,
                features=np.asarray(features)[active].tolist(), x=x[:, active], t=t,
                meta=c['meta'], sheets=[m['bifolium'] for m in c['meta']], h_ids=h_ids,
                steps=data['design']['steps'][c['cohort']], source=STUDY/c['matrices'][mode+'_raw']))
    return out


def verify(cases):
    old = json.loads((STUDY/'scegret_independent_audit.json').read_text())
    assert old['status'] == 'passed'
    max_error = 0.
    checks = 0
    for c in cases:
        for mask in [np.ones(len(c['t']), bool), ~((c['t'] > .53) & (c['t'] < .75))]:
            t, x = c['t'][mask], c['x'][mask]
            scorer = Scorer(t, c['steps'])
            reference = []
            for steps in c['steps']:
                for scale in ['per_million', 'raw']:
                    matrix = (x/x.sum(axis=1)[:, None]*1e6 if scale == 'per_million' else x).T
                    e = audit.expected(matrix, c['features'], np.asarray(c['ids'])[mask].tolist(), t, steps)
                    if e['sufficient']:
                        reference.extend(e['adjacent'])
            error = float(np.max(np.abs(np.asarray(reference)-scorer.score(x))))
            max_error = max(max_error, error)
            checks += 1
    assert max_error < 1e-12
    out = dict(status='passed', fast_scorer_comparisons=checks, max_absolute_error=max_error,
        original_R_audit_sha256=sha(STUDY/'scegret_independent_audit.json'),
        original_R_results_sha256=sha(STUDY/'scegret_results.json'))
    write('scorer_verification.json', out)


def calibrate(c, bandwidth, model, key):
    gen = Generator(c['x'], c['t'], c['sheets'], bandwidth)
    scorer = Scorer(c['t'], c['steps'], c['h_ids'])
    observed = scorer.score(c['x'])
    all_sets, diagnostic = [], []
    for stage, n in enumerate([CONFIG['fit_reps'], CONFIG['calibration_reps'], CONFIG['evaluation_reps']]):
        random = rng(*key, stage)
        values = []
        for i in range(n):
            x = gen.draw(model, random)
            values.append(scorer.score(x))
            if stage == 1 and i < 100:
                diagnostic.append(gen.diagnostics(x))
        all_sets.append(np.asarray(values))
    fit, calibration, evaluation = all_sets
    means, sd = fit.mean(axis=0), fit.std(axis=0, ddof=1)
    assert np.all(sd > 0)
    calibration_z = (calibration-means)/sd
    max_z = calibration_z.max(axis=1)
    evaluated_max_z = ((evaluation-means)/sd).max(axis=1)
    observed_z = (observed-means)/sd
    threshold = float(np.quantile(max_z, .95, method='higher'))
    boundaries = []
    for j, label in enumerate(scorer.labels):
        scan_tail = tail(max_z, observed_z[j])
        boundaries.append(dict(label, observed=float(observed[j]), null_mean=float(means[j]), null_sd=float(sd[j]),
            z=float(observed_z[j]), local_tail=tail(calibration[:, j], observed[j]),
            scan_tail=scan_tail, whole_study_tail=min(1., 4*scan_tail),
            reference95=np.quantile(calibration[:, j], [.025, .975]).tolist()))
    observed_max = float(observed_z.max())
    diag = {name: dict(observed=gen.observed_diag[name], simulated_median=float(np.median([r[name] for r in diagnostic])),
        simulated95=np.quantile([r[name] for r in diagnostic], [.025, .975]).tolist()) for name in gen.observed_diag}
    result = dict(status='complete', cohort=c['cohort'], channel=c['mode'], bandwidth=bandwidth, model=model,
        alpha=gen.alpha, excess_variance_rho=gen.rho, diagnostics=diag,
        threshold95_max_standardized_score=threshold, heldout_false_alarm=rate(evaluated_max_z > threshold),
        observed_max_z=observed_max, overall_scan_tail=tail(max_z, observed_max),
        overall_whole_study_tail=min(1., 4*tail(max_z, observed_max)), boundaries=boundaries,
        calibration_max_z=max_z.tolist(), heldout_max_z=evaluated_max_z.tolist())
    return result


def power(c, model, key):
    gen = Generator(c['x'], c['t'], c['sheets'], .12)
    t, n, steps = c['t'], len(c['t']), c['steps'][1]
    endpoint = np.flatnonzero((t == t.min()) | (t == t.max()))
    interior = np.setdiff1d(np.arange(n), endpoint)
    out = []
    for center_i, center in enumerate(CONFIG['gap_centers']):
        for width_i, width in enumerate(CONFIG['gap_widths']):
            gap = (t > max(0, center-width/2)) & (t < min(1, center+width/2))
            gap[endpoint] = False
            removed = int(gap.sum())
            gap_scorer = Scorer(t[~gap], [steps], scales=('per_million',))
            def measure(x, mask):
                sc = Scorer(t[mask], [steps], scales=('per_million',))
                scores = sc.score(x[mask])
                return float(scores.max()) if len(scores) else None
            def random_mask(random):
                m = np.ones(n, bool)
                m[random.choice(interior, removed, replace=False)] = False
                return m
            random = rng(*key, 40, center_i, width_i, 0)
            null_scores = []
            for _ in range(CONFIG['power_calibration_reps']):
                x = gen.draw(model, random)
                score = measure(x, random_mask(random))
                if score is not None:
                    null_scores.append(score)
            threshold = float(np.quantile(null_scores, .95, method='higher'))
            random = rng(*key, 40, center_i, width_i, 1)
            gap_scores, random_scores = [], []
            for _ in range(CONFIG['power_evaluation_reps']):
                x = gen.draw(model, random)
                ss = gap_scorer.score(x[~gap])
                gap_scores.append(float(ss.max()) if len(ss) else None)
                random_scores.append(measure(x, random_mask(random)))
            gv = np.array([s for s in gap_scores if s is not None])
            rv = np.array([s for s in random_scores if s is not None])
            out.append(dict(center=center, width=width, removed=removed, retained=n-removed,
                steps=steps, retained_gap_bin_ids=gap_scorer.views[0][2].tolist(),
                threshold95_max_adjacent=threshold, gap_detection=rate(gv > threshold),
                random_deletion_false_alarm=rate(rv > threshold),
                gap_score_median=float(np.median(gv)), random_score_median=float(np.median(rv)),
                calibration_evaluable=len(null_scores), gap_evaluable=len(gv), random_evaluable=len(rv)))
    return dict(cohort=c['cohort'], channel=c['mode'], model=model, bandwidth=.12, conditions=out)


def attribution(cases):
    result = []
    for c in cases:
        if c['cohort'] != 'paragraph207':
            continue
        h58 = np.array([s in ['f58r', 'f58v'] for s in c['ids']])
        astro = np.array([s in ['f67v2', 'f69r', 'f70r2'] for s in c['ids']])
        hsheet = np.array([s == 'H:58-65' for s in c['sheets']])
        for name, omit in [('full', np.zeros(len(h58), bool)), ('omit_f58rv', h58),
                           ('omit_other_three_H_bin_profiles', astro), ('omit_H_sheet', hsheet)]:
            mask = ~omit
            x, ids, t = c['x'][mask], np.array(c['ids'])[mask], c['t'][mask]
            e = audit.expected((x/x.sum(axis=1)[:, None]*1e6).T, c['features'], ids.tolist(), t, 12)
            pairs = list(zip(e['retained_bins'][:-1], e['retained_bins'][1:]))
            top = np.argsort(e['contributions'].sum(axis=1))[-12:][::-1]
            result.append(dict(channel=c['mode'], omission=name, omitted=np.array(c['ids'])[omit].tolist(),
                occupancy=e['occupancy'].tolist(), retained_bins=e['retained_bins'].tolist(),
                adjacent=[dict(left=int(a), right=int(b), score=float(v)) for (a,b),v in zip(pairs,e['adjacent'])],
                bin8_members=e['bin_observations'].get('8', []), bin8_eligible=8 in e['retained_bins'],
                top_overall_features=[e['features'][i] for i in top]))
    write('attribution.json', dict(status='complete', runs=result,
        warning='Omissions that make a bin ineligible do not isolate causal responsibility.'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['all','verify','null','power'], default='all')
    parser.add_argument('--case', type=int)
    args = parser.parse_args()
    cases = load()
    write('configuration.json', dict(CONFIG, design_sha256=sha(HERE/'DESIGN.md'),
        code_sha256=sha(Path(__file__)), input_sha256={c['source'].name:sha(c['source']) for c in cases},
        original_positions_sha256=sha(STUDY/'scegret_inputs.json')))
    verify(cases)
    attribution(cases)
    if args.stage == 'verify':
        print('Verification and attribution complete.', flush=True)
        return
    start = time.monotonic()
    with threadpool_limits(limits=1):
        for case_i, c in enumerate(cases):
            if args.case is not None and case_i != args.case:
                continue
            for model_i, model in enumerate(CONFIG['models']):
                if args.stage in ['all','null']:
                    for bw_i, bandwidth in enumerate(CONFIG['bandwidths']):
                        filename = f"null_{c['cohort']}_{c['mode']}_{model}_h{bandwidth:.2f}.json"
                        if (HERE/filename).exists():
                            cached = json.loads((HERE/filename).read_text())
                            assert cached['code_sha256'] == sha(Path(__file__))
                            assert cached['design_sha256'] == sha(HERE/'DESIGN.md')
                            print('Cached '+filename, flush=True)
                            continue
                        out = calibrate(c, bandwidth, model, [case_i,model_i,bw_i])
                        out['code_sha256'] = sha(Path(__file__))
                        out['design_sha256'] = sha(HERE/'DESIGN.md')
                        write(filename, out)
                        print(json.dumps(dict(file=filename, seconds=round(time.monotonic()-start),
                            false_alarm=out['heldout_false_alarm']['rate'], scan_tail=out['overall_scan_tail'])), flush=True)
                if args.stage in ['all','power']:
                    filename = f"power_{c['cohort']}_{c['mode']}_{model}.json"
                    if (HERE/filename).exists():
                        cached = json.loads((HERE/filename).read_text())
                        assert cached['code_sha256'] == sha(Path(__file__))
                        assert cached['design_sha256'] == sha(HERE/'DESIGN.md')
                        print('Cached '+filename, flush=True)
                        continue
                    out = power(c, model, [case_i,model_i,1])
                    out['code_sha256'] = sha(Path(__file__))
                    out['design_sha256'] = sha(HERE/'DESIGN.md')
                    write(filename, out)
                    print(json.dumps(dict(file=filename, seconds=round(time.monotonic()-start))), flush=True)


if __name__ == '__main__':
    main()
