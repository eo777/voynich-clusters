"""Power of the exact multiplicity-adjusted test, using fixed-count latent gaps."""
import json
import time
from pathlib import Path
import numpy as np
import run_calibration as cal

HERE = Path(__file__).resolve().parent


def main():
    start = time.monotonic()
    with cal.threadpool_limits(limits=1):
        for case_i, c in enumerate(cal.load()):
            scorer = cal.Scorer(c['t'], c['steps'], c['h_ids'])
            for model_i, model in enumerate(cal.CONFIG['models']):
                source = HERE/f"null_{c['cohort']}_{c['mode']}_{model}_h0.12.json"
                ref = json.loads(source.read_text())
                assert ref['code_sha256'] == cal.sha(HERE/'run_calibration.py')
                means = np.array([r['null_mean'] for r in ref['boundaries']])
                sd = np.array([r['null_sd'] for r in ref['boundaries']])
                max_ref = np.array(ref['calibration_max_z'])
                def detected(max_z):
                    return cal.tail(max_ref, max_z)*4 <= .05
                results = []
                for ci, center in enumerate(cal.CONFIG['gap_centers']):
                    for wi, width in enumerate(cal.CONFIG['gap_widths']):
                        gen = cal.Generator(c['x'], c['t'], c['sheets'], .12)
                        u = (1-width)*c['t']+width*(c['t'] >= center)
                        w = np.exp(-.5*((u[:, None]-c['t'][None, :])/.12)**2)
                        gen.mu = w @ gen.p/w.sum(axis=1)[:, None]
                        gen.root_mu = np.sqrt(gen.mu)
                        random = cal.rng(case_i, model_i, 91, ci, wi)
                        max_scores = []
                        for _ in range(300):
                            scores = scorer.score(gen.draw(model, random))
                            max_scores.append(float(np.max((scores-means)/sd)))
                        results.append(dict(center=center, width=width,
                            skipped_latent_interval=[(1-width)*center, (1-width)*center+width],
                            detection=cal.rate([detected(z) for z in max_scores]),
                            median_max_z=float(np.median(max_scores))))
                out = dict(status='complete', cohort=c['cohort'], channel=c['mode'], model=model, bandwidth=.12,
                    threshold_rule='4 * scan_simulation_tail <= 0.05',
                    continuous_heldout_false_alarm=cal.rate([detected(z) for z in ref['heldout_max_z']]),
                    conditions=results, source_reference_sha256=cal.sha(source),
                    code_sha256=cal.sha(Path(__file__)), design_sha256=cal.sha(HERE/'MATCHED_DESIGN.md'))
                name = f"matched_{c['cohort']}_{c['mode']}_{model}.json"
                cal.write(name, out)
                print(json.dumps(dict(file=name, seconds=round(time.monotonic()-start))), flush=True)


if __name__ == '__main__':
    main()
