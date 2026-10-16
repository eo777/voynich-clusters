"""Matched RF ensemble using the frozen original growth/refinement/readout code."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import copy
import hashlib
import importlib.metadata
import json
import sys
import time
HERE=Path(__file__).resolve().parent
TEAM=HERE.parent
ROOT=TEAM.parents[2]
OLD=TEAM/'golf_published_stability_20261016'
sys.path.insert(0,str(HERE/'_deps'))
sys.path.insert(1,str(OLD))
import numpy as np
from threadpoolctl import threadpool_limits
import ensemble as engine
core=engine.core
sha=engine.sha
dump=engine.dump

def prepare():
    info=json.loads((HERE/'inputs.json').read_text())
    assert sha(HERE/'inputs.npz')==info['data_sha256']
    for p,d in info['source_sha256'].items():assert sha(ROOT/p)==d,p
    a=np.load(HERE/'inputs.npz');contexts={};paths=[]
    baseline=json.loads((OLD/'ensemble_index.json').read_text())
    for name in baseline['setting_files']:
        path=OLD/name
        assert sha(path)==baseline['result_sha256'][name],name
        d=json.loads(path.read_text())['design'];prefix='page' if d['cohort']=='paragraph207' else 'pool'
        x=a[prefix+'_coords'];dist=a[prefix+'_distance'];meta=info['meta' if prefix=='page' else 'pool_meta']
        assert [m['id'] for m in meta]==[m['id'] for m in d['meta']]
        assert [m['group'] for m in meta]==[m['group'] for m in d['meta']]
        w=np.array(d['original_weights']);labels=np.array([m['group'] for m in meta])
        centers={g:w[labels==g]@x[labels==g]/w[labels==g].sum() for g in core.ORDER}
        starts=[]
        for old in d['starts']:
            s=copy.deepcopy(old);initial=np.array([centers[g] for g in s['anchor_order']])
            s['zl_initial_sha256']=s.pop('initial_sha256')
            s['zl_selection_objective']=s.pop('prior_objective')
            s['initial']=initial;s['initial_sha256']=hashlib.sha256(initial.tobytes()).hexdigest()
            starts.append(s)
        key=name.removeprefix('ensemble_').removesuffix('.json')
        contexts[key]=dict(key=key,cohort=d['cohort'],sampler=d['sampler'],penalty=d['penalty'],meta=meta,
            x=x,weights=w,blocks=np.array(d['blocks']),labels=labels,starts=starts,masks=np.array(d['masks'],dtype=bool),
            mask_seed=d['mask_seed'],mask_unit=d['mask_unit'],
            geometry_sha256=hashlib.sha256(x.tobytes()).hexdigest(),distance_sha256=hashlib.sha256(dist.tobytes()).hexdigest())
        paths.append(path)
    paths += [Path(__file__),HERE/'DESIGN.md',HERE/'inputs.json',HERE/'inputs.npz',OLD/'ensemble.py',
        OLD/'ENSEMBLE_DESIGN.md',OLD/'ensemble_index.json',Path(core.__file__),
        TEAM/'whole_manuscript_topology_20261015/common.py',TEAM/'whole_manuscript_topology_20261015/models.py',
        HERE/'_deps/elpigraph/_BaseElPiWrapper.py',HERE/'_deps/elpigraph/_topologies.py']
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths};sources.update(info['source_sha256'])
    return contexts,sources

def init_worker():
    threadpool_limits(limits=1)
    engine._STATE=prepare()[0]

def fit_one(task):return engine.fit_one(task)

def pilot():
    # Execute the unchanged original native-mask verification on RF coordinates.
    engine.prepare=prepare;engine.HERE=HERE
    with threadpool_limits(limits=1):
        engine.native_pilot()
        contexts,sources=prepare()
        original=core.load_data();c=copy.deepcopy(contexts['paragraph207_point_l010'])
        c['x']=core.COMMON.embed(original['ds']['combined'])[0]
        x=c['x'];w=c['weights'];labels=c['labels']
        centers={g:w[labels==g]@x[labels==g]/w[labels==g].sum() for g in core.ORDER}
        for s in c['starts']:s['initial']=np.array([centers[g] for g in s['anchor_order']])
        engine._STATE={'zl_replay':c}
        _,got=fit_one(('zl_replay',0,0))
        saved=json.loads((OLD/'ensemble_paragraph207_point_l010.json').read_text())['runs'][0]
        errors={k:float(abs(np.asarray(got[k])-np.asarray(saved[k])).max()) for k in ['positions','objective','node_mse','edge_energy','star_energy']}
        assert max(errors.values())<1e-10,errors
        engine._STATE=contexts
        checks=[]
        for key in ['paragraph207_sheet_l050','pool50_point_l010']:
            _,r=fit_one((key,0,0))
            checks.append(dict(setting=key,converged=r['converged'],objective=r['objective'],seconds=r['seconds'],
                positions=r['positions'],full_node_sha256=r['full_node_sha256']))
            assert r['converged']
    versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn','threadpoolctl','elpigraph-python','numba','llvmlite','networkx','igraph','shapely','texttable']}
    dump(HERE/'execution_verification.json',dict(status='passed',source_sha256=sources,
        original_zl_replay_max_errors=errors,rf_checks=checks,additional_growth_calls=3,versions=versions))
    print(json.dumps(dict(zl_replay_errors=errors,rf_seconds=[c['seconds'] for c in checks],status='passed')),flush=True)

def run(workers):
    assert json.loads((HERE/'ensemble_pilot.json').read_text())['status']=='passed'
    verified=json.loads((HERE/'execution_verification.json').read_text());assert verified['status']=='passed'
    contexts,sources=prepare();assert sources==verified['source_sha256']
    output={};complete=set()
    for key,c in contexts.items():
        path=HERE/f'ensemble_{key}.json'
        if path.exists():
            obj=json.loads(path.read_text());assert obj['source_sha256']==sources
        else:
            design=engine.metadata(c)
            design.update(representation='RF1b-er / basic EVA characters',matched_zl_setting=f'ensemble_{key}.json')
            obj=dict(status='running',source_sha256=sources,design=design,runs=[])
        output[key]=obj
        complete.update((key,r['replicate'],r['start_index']) for r in obj['runs'])
    tasks=[(key,rep,s) for rep in range(100) for key in contexts for s in range(8) if (key,rep,s) not in complete]
    index=dict(status='running',planned_fits=4800,completed_fits=len(complete),source_sha256=sources,
        setting_files=[f'ensemble_{key}.json' for key in contexts],workers=workers,blas_threads_per_worker=1,
        pilot_growth_calls=2,additional_verification_growth_calls=3,
        representation='RF1b-er / basic EVA characters',matched_masks_and_anchor_permutations=True,
        interpretation='Equal-design retained-fit sensitivity; not chronology probabilities. Changes transcription and glyph units.')
    began=time.time();prior_elapsed=0.
    if (HERE/'ensemble_index.json').exists():prior_elapsed=json.loads((HERE/'ensemble_index.json').read_text()).get('elapsed_seconds',0.)
    def checkpoint():
        for key,obj in output.items():
            obj['runs'].sort(key=lambda r:(r['replicate'],r['start_index']))
            obj['status']='complete' if len(obj['runs'])==800 else 'running'
            dump(HERE/f'ensemble_{key}.json',obj)
        index.update(elapsed_seconds=prior_elapsed+time.time()-began,
            converged_fits=sum(r['converged'] for o in output.values() for r in o['runs']))
        dump(HERE/'ensemble_index.json',index)
    checkpoint()
    with ProcessPoolExecutor(max_workers=workers,initializer=init_worker) as pool:
        futures=[pool.submit(fit_one,t) for t in tasks]
        for future in as_completed(futures):
            key,r=future.result();output[key]['runs'].append(r);index['completed_fits']+=1
            if index['completed_fits']%100==0:
                checkpoint()
                print(f"{index['completed_fits']}/4800 fits; {index['elapsed_seconds']:.1f}s; converged {index['converged_fits']}",flush=True)
    assert index['completed_fits']==4800
    index['status']='complete';checkpoint()
    index['result_sha256']={name:sha(HERE/name) for name in index['setting_files']}
    index['retained_result_bytes']=sum((HERE/name).stat().st_size for name in index['setting_files'])
    assert index['retained_result_bytes']<100_000_000
    dump(HERE/'ensemble_index.json',index)
    print(json.dumps({k:v for k,v in index.items() if k not in ['source_sha256','result_sha256']},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');p.add_argument('--workers',type=int,default=6)
    args=p.parse_args()
    if args.pilot:pilot()
    else:run(args.workers)
