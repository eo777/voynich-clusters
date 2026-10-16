"""Published ElPiGraph resampled ensembles with explicit auditable weights/masks."""
from pathlib import Path
import os
for _name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_name] = '1'
import argparse
import contextlib
import hashlib
import io
import itertools
import json
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent
sys.path.insert(0, str(HERE/'_deps'))
sys.path.insert(1, str(TEAM/'golf_stability_20261016'))
import numpy as np
import core
from threadpoolctl import threadpool_limits

GROUPS = core.ORDER
PAIRS = [('1a','1b'),('3a','3b'),('4a','4c')]
TOL = 1e-7
SEEDS = {'paragraph_point':2026101601,'sheet':2026101602,'pilot':2026101603}
_STATE = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    payload = json.dumps(core.jsonable(data), separators=(',',':'), allow_nan=False)
    assert len(payload.encode('utf-8')) < 12_000_000, (path, len(payload))
    temporary = Path(str(path)+'.tmp')
    temporary.write_text(payload,encoding='utf-8')
    temporary.replace(path)


def order(t, labels, included=None):
    use = np.ones(len(t),dtype=bool) if included is None else included
    medians = {g:float(np.median(t[(labels==g)&use])) if np.any((labels==g)&use) else None for g in GROUPS}
    pairs = {}
    for a,b in PAIRS:
        delta = None if medians[a] is None or medians[b] is None else medians[a]-medians[b]
        pairs[a+'|'+b] = {'delta':delta,'sign':None if delta is None else -1 if delta < -TOL else 1 if delta > TOL else 0}
    return {'group_medians':medians,'contrasts':pairs}


def prepare():
    data = core.load_data()
    prior = json.loads((TEAM/'golf_influence_20261016/b_permutations.json').read_text())
    contexts = {}
    pool_blocks = np.array([m['bifolium'] for m in data['pool_meta']])
    sheets = np.array(sorted(pool_blocks))
    sheet_draws = np.random.RandomState(SEEDS['sheet']).random_sample((100,len(sheets))) <= .9
    for name,meta,d,blocks in [
        ('paragraph207',data['meta'],data['ds']['combined'],np.asarray(data['blocks'])),
        ('pool50',data['pool_meta'],data['pool_ds']['combined'],pool_blocks)]:
        x, eig = core.COMMON.embed(d)
        weights = core.COMMON.sheet_weights(blocks)
        labels = np.array([m['group'] for m in meta])
        centers = {g:weights[labels==g]@x[labels==g]/weights[labels==g].sum() for g in GROUPS}
        sheet_masks = sheet_draws[:,np.array([np.flatnonzero(sheets==b)[0] for b in blocks])]
        if name=='paragraph207':
            point_masks = np.random.RandomState(SEEDS['paragraph_point']).random_sample((100,len(meta))) <= .9
        else:
            point_masks = sheet_masks.copy()
        for lam in [.01,.05]:
            s = next(s for s in prior['settings'] if s['cohort']==name and s['penalty']==lam)
            starts = []
            for hs,ts in itertools.product([-1,1],repeat=2):
                candidates = [r for r in s['runs'] if r['early_relation']==-1 and r['b_relations']['3a|3b']==hs and r['b_relations']['4a|4c']==ts]
                representative = min(candidates,key=lambda r:r['objective'])
                for early in [0,1]:
                    prefix = ['1a','1b','1.5a','2a'] if early==0 else ['1b','1a','1.5a','2a']
                    anchor_order = prefix+representative['initial_B_order']
                    initial = np.array([centers[g] for g in anchor_order])
                    starts.append(dict(id=f'e{early}_h{hs:+d}_t{ts:+d}', prior_candidate_id=representative['id'],
                        anchor_order=anchor_order, prior_final_signs=[hs,ts],
                        prior_objective=representative['objective'], initial=initial,
                        initial_sha256=hashlib.sha256(initial.tobytes()).hexdigest()))
            assert len({tuple(s['anchor_order']) for s in starts})==8
            for sampler in (['point','sheet'] if name=='paragraph207' else ['point']):
                key=f'{name}_{sampler}_l{int(lam*1000):03d}'
                contexts[key]=dict(key=key,cohort=name,sampler=sampler,penalty=lam,meta=meta,
                    x=x,weights=weights,blocks=blocks,labels=labels,starts=starts,
                    masks=point_masks if sampler=='point' else sheet_masks,
                    mask_seed=SEEDS['paragraph_point'] if name=='paragraph207' and sampler=='point' else SEEDS['sheet'],
                    mask_unit='row' if sampler=='point' else 'physical_bifolium',
                    geometry_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
                    distance_sha256=hashlib.sha256(np.asarray(d).tobytes()).hexdigest())
    paths=[Path(__file__),HERE/'ENSEMBLE_DESIGN.md',TEAM/'golf_stability_20261016/core.py',
        TEAM/'golf_influence_20261016/b_permutations.json',
        HERE/'_deps/elpigraph/_BaseElPiWrapper.py',HERE/'_deps/elpigraph/_topologies.py',*data['sources']]
    paths += [core.COMMON.DATA/'profiles.csv']
    paths += [core.COMMON.DATA/f'pairwise/shared-{mode}-js.csv' for mode in ['word','line','combined']]
    hashes={str(p.relative_to(TEAM.parents[2])):sha(p) for p in dict.fromkeys(paths)}
    return contexts,hashes


def metadata(c):
    return dict(cohort=c['cohort'],sampler=c['sampler'],penalty=c['penalty'],mu=.1,nodes=28,
        meta=c['meta'],blocks=c['blocks'],original_weights=c['weights'],positive_dimensions=c['x'].shape[1],
        reference_coords3=np.round(c['x'][:,:3],10),mask_seed=c['mask_seed'],mask_unit=c['mask_unit'],
        masks=c['masks'].astype(int),geometry_sha256=c['geometry_sha256'],distance_sha256=c['distance_sha256'],
        starts=[{k:v for k,v in s.items() if k!='initial'} for s in c['starts']],
        tie_tolerance=TOL,position_readout='closest-segment normalized arc length; all original rows projected',
        primary_frequency_scope='row retained; pair retained jointly',mask_probability=.9,
        weight_rule='subset original equal-sheet weights and renormalize',expected_runs=800)


def init_worker():
    global _STATE
    _STATE = prepare()[0]
    threadpool_limits(limits=1)


def fit_one(task):
    key,replicate,start_index = task
    c = _STATE[key]
    mask = c['masks'][replicate]
    w = core.normalized(c['weights'][mask])
    assert len(w)==int(mask.sum()) and np.isclose(w.sum(),1)
    start=c['starts'][start_index]
    began=time.time()
    grown=core.grow_curve(c['x'][mask],w,start['initial'],lam=c['penalty'])
    grown=core.orient_curve(grown,c['x'],c['meta'])
    model=core.refine_curve(c['x'][mask],w,grown['nodes'],grown['edges'],lam=c['penalty'])
    projection=core.project_curve(c['x'],model)
    t=projection['positions']
    full=core.energy(c['x'],c['weights'],model['nodes'],model['edges'],c['penalty'],.1)
    assert abs(model['objective']-(model['node_mse']+model['edge_energy']+model['star_energy']))<1e-12
    result=dict(replicate=replicate,start_index=start_index,start_id=start['id'],retained=int(mask.sum()),
        retained_weight_before_normalization=float(c['weights'][mask].sum()),
        normalized_weight_sha256=hashlib.sha256(w.tobytes()).hexdigest(),
        positions=t,nodes3=np.round(model['nodes'][:,:3],10),
        full_node_sha256=hashlib.sha256(model['nodes'].tobytes()).hexdigest(),
        group_readout_all=order(t,c['labels']),group_readout_retained=order(t,c['labels'],mask),
        objective=model['objective'],node_mse=model['node_mse'],edge_energy=model['edge_energy'],star_energy=model['star_energy'],
        full_reference_objective=full['objective'],growth_objective=grown['objective'],
        training_segment_mse=float(w@projection['errors'][mask]),
        full_reference_segment_mse=float(c['weights']@projection['errors']),path_length=projection['arc_length'],
        converged=model['converged'],iterations=model['iterations'],
        assignment_changes=int(sum(model['assignment_changes'])),
        maximum_energy_increase=max(0.,float(np.diff(model['energy_trace']).max())),
        warning_count=grown.get('warning_count',0),seconds=time.time()-began)
    return key,result


def native_pilot():
    contexts,hashes=prepare()
    c=contexts['pool50_point_l010']
    initial=c['starts'][0]['initial']
    seed=SEEDS['pilot']
    mask=np.random.RandomState(seed).random_sample(len(c['x']))<=.9
    np.random.seed(seed)
    began=time.time()
    with contextlib.redirect_stdout(io.StringIO()),warnings.catch_warnings(record=True):
        native=core.elpigraph.computeElasticPrincipalCurve(c['x'],NumNodes=28,Lambda=.01,Mu=.1,
            InitNodePositions=initial.copy(),InitEdges=np.array([[i,i+1] for i in range(8)]),
            Do_PCA=False,CenterData=False,PointWeights=None,FinalEnergy='Base',alpha=0,beta=0,
            MaxNumberOfIterations=100,eps=1e-5,n_cores=1,nReps=1,ProbPoint=.9,EmbPointProb=1,
            SampleIC=False,AvoidResampling=True,ComputeMSEP=False,verbose=False,ShowTimer=False,DisplayWarnings=False)[0]
    native_seconds=time.time()-began
    began=time.time()
    explicit=core.grow_curve(c['x'][mask],np.ones(mask.sum()),initial,lam=.01)
    explicit_seconds=time.time()-began
    native_model=dict(nodes=native['NodePositions'],edges=native['Edges'][0].astype(int))
    a=core.orient_curve(native_model,c['x'],c['meta'])
    b=core.orient_curve(explicit,c['x'],c['meta'])
    node_error=float(abs(a['nodes']-b['nodes']).max())
    native_energy=core.energy(c['x'][mask],np.ones(mask.sum())/mask.sum(),a['nodes'],a['edges'],.01,.1)
    energy_error=float(abs(native_energy['objective']-explicit['objective']))
    assert node_error<1e-10 and energy_error<1e-12,(node_error,energy_error)
    audit=dict(status='passed',source_sha256=hashes,seed=seed,retained=int(mask.sum()),mask=mask.astype(int),
        growth_calls=2,native_seconds=native_seconds,explicit_seconds=explicit_seconds,
        maximum_node_difference=node_error,objective_difference=energy_error,
        native_probability=.9,native_replicates=1,no_average_graph=True,
        weight_wrapper_finding='X is selected at wrapper line 357; PointWeights forwarded unchanged at line 393. Explicit resampling subsets and renormalizes weights.',
        storage_estimate_MB=60,estimated_wall_minutes_two_workers=[20,40])
    dump(HERE/'ensemble_pilot.json',audit)
    print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','mask']},indent=2),flush=True)


def run(workers=2):
    assert json.loads((HERE/'ensemble_pilot.json').read_text())['status']=='passed'
    contexts,hashes=prepare()
    output={key:dict(status='running',source_sha256=hashes,design=metadata(c),runs=[]) for key,c in contexts.items()}
    index=dict(status='running',source_sha256=hashes,planned_fits=4800,completed_fits=0,
        pilot_growth_calls=2,workers=workers,blas_threads_per_worker=1,
        published_method='ElPiGraph 2020 section 2.10 / Fig1C random-resampled graph ensembles',
        setting_files=[f'ensemble_{key}.json' for key in contexts],
        interpretation='Sampling and algorithm sensitivity under the declared design; not historical probability.',
        primary_scope='retained-row and jointly retained pair frequencies; all-reference projections retained separately')
    for key,obj in output.items():dump(HERE/f'ensemble_{key}.json',obj)
    dump(HERE/'ensemble_index.json',index)
    tasks=[(key,replicate,start) for replicate in range(100) for key in contexts for start in range(8)]
    began=time.time()
    with ProcessPoolExecutor(max_workers=workers,initializer=init_worker) as pool:
        futures=[pool.submit(fit_one,t) for t in tasks]
        for future in as_completed(futures):
            key,record=future.result()
            output[key]['runs'].append(record)
            index['completed_fits']+=1
            if len(output[key]['runs'])%40==0 or index['completed_fits']==4800:
                obj=output[key]
                obj['runs'].sort(key=lambda r:(r['replicate'],r['start_index']))
                if len(obj['runs'])==800:obj['status']='complete'
                dump(HERE/f'ensemble_{key}.json',obj)
            if index['completed_fits']%100==0:
                index['elapsed_seconds']=time.time()-began
                index['converged_fits']=sum(r['converged'] for obj in output.values() for r in obj['runs'])
                dump(HERE/'ensemble_index.json',index)
                print(f"{index['completed_fits']}/4800 fits; {index['elapsed_seconds']:.1f}s; converged {index['converged_fits']}",flush=True)
    for key,obj in output.items():
        assert len(obj['runs'])==800
        obj['status']='complete';obj['runs'].sort(key=lambda r:(r['replicate'],r['start_index']))
        dump(HERE/f'ensemble_{key}.json',obj)
    index['status']='complete';index['elapsed_seconds']=time.time()-began
    index['result_sha256']={f'ensemble_{key}.json':sha(HERE/f'ensemble_{key}.json') for key in contexts}
    index['retained_result_bytes']=sum((HERE/f'ensemble_{key}.json').stat().st_size for key in contexts)
    assert index['retained_result_bytes']<100_000_000
    dump(HERE/'ensemble_index.json',index)
    print(json.dumps(index,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        if args.pilot:native_pilot()
        else:run(args.workers)
