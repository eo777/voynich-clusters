"""Read-only reconstruction from the packaged inputs; NumPy and SciPy."""
from pathlib import Path
import csv
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE/'replication'
TEAM = ROOT/'solve/botanical_enrichment/team_golf'
STUDY = TEAM/'golf_published_stability_20261016'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def main():
    common = module('package_common', TEAM/'whole_manuscript_topology_20261015/common.py')
    driver = module('package_counts', TEAM/'b_horseshoe_20261015/drivers.py')
    data = common.load_data()
    alphabet = [r['unit'] for r in common.rows(common.DATA/'alphabet.csv')]
    maximum = 0.; matrix_checks = 0
    for cohort, meta, distances in [('paragraph207',data['meta'],data['ds']),
                                    ('pool50',data['pool_meta'],data['pool_ds'])]:
        recovered = []
        for mode in ('word','line'):
            counts = driver.read_counts(mode,[m['profile_id'] for m in meta],alphabet)
            js = driver.REPRO.weighted_js(counts)[0]
            error = float(np.max(np.abs(js-distances[mode])))
            assert error < 1e-12, (cohort,mode,error)
            maximum = max(maximum,error); recovered.append(js)
            flat = counts.reshape(len(meta),-1)
            for scale in ('raw','per_million'):
                with (STUDY/f'scegret_{cohort}_{mode}_{scale}.csv').open(newline='',encoding='utf-8') as f:
                    rows = list(csv.reader(f))
                assert rows[0][1:] == [m['id'] for m in meta]
                observed = np.array([r[1:] for r in rows[1:]],dtype=float).T
                expected = flat if scale=='raw' else flat/flat.sum(1)[:,None]*1e6
                np.testing.assert_allclose(observed,expected,rtol=1e-10,atol=1e-10)
                matrix_checks += 1
        np.testing.assert_allclose(sum(recovered)/2,distances['combined'],rtol=0,atol=1e-12)
    print(json.dumps(dict(status='passed',pages=len(data['meta']),pools=len(data['pool_meta']),
        maximum_distance_error=maximum,feature_matrices=matrix_checks)))


if __name__ == '__main__': main()
