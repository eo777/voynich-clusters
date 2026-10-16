"""Two post-production weighted package replays with independent calculations."""
from pathlib import Path
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import sys
import time
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name] = '1'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE/'_deps'))
sys.path.insert(1, str(HERE.parent/'golf_stability_20261016'))
import numpy as np
import elpigraph
import core
from threadpoolctl import threadpool_limits


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def laplacian(n, lam):
    matrix = np.zeros((n,n))
    for j in range(n-1):
        v = np.zeros(n); v[j:j+2] = [-1,1]
        matrix += lam*np.outer(v,v)
    for j in range(n-2):
        v = np.zeros(n); v[j:j+3] = [-.5,1,-.5]
        matrix += .1*np.outer(v,v)
    return matrix


def energy(x,w,y,lam):
    differences = ((x[:,None,:]-y[None,:,:])**2).sum(axis=2)
    node = float(w@differences.min(axis=1))
    edge = float(lam*np.sum(np.diff(y,axis=0)**2))
    star = float(.1*np.sum((y[1:-1]-.5*(y[:-2]+y[2:]))**2))
    return np.array([node+edge+star,node,edge,star]),differences.argmin(axis=1)


def orient(y,edges,x,labels):
    adjacent=[[] for _ in range(len(y))]
    for a,b in edges:
        adjacent[a].append(b);adjacent[b].append(a)
    ends=[j for j,a in enumerate(adjacent) if len(a)==1]
    assert len(ends)==2 and max(map(len,adjacent))<=2
    order=[min(ends)]
    while len(order)<len(y):
        order.append(next(j for j in adjacent[order[-1]] if j not in order))
    y=y[order]
    early=x[np.isin(labels,['1a','1b'])].mean(axis=0)
    late=x[np.isin(labels,['4a','4c'])].mean(axis=0)
    direct=np.sum((y[0]-early)**2)+np.sum((y[-1]-late)**2)
    reverse=np.sum((y[-1]-early)**2)+np.sum((y[0]-late)**2)
    return y[::-1].copy() if reverse<direct else y.copy()


def polish(x,w,y,lam):
    lap=laplacian(len(y),lam)
    _,assignment=energy(x,w,y,lam)
    for iteration in range(1000):
        masses=np.array([w[assignment==j].sum() for j in range(len(y))])
        sums=np.array([(w[assignment==j,None]*x[assignment==j]).sum(axis=0) for j in range(len(y))])
        updated=np.linalg.solve(lap+np.diag(masses),sums)
        _,new_assignment=energy(x,w,updated,lam)
        if np.array_equal(assignment,new_assignment) and np.linalg.norm(updated-y,axis=1).max()<=1e-9:
            return updated,iteration+1
        y,assignment=updated,new_assignment
    raise AssertionError('Independent polish failed')


def project(x,y):
    lengths=np.linalg.norm(np.diff(y,axis=0),axis=1)
    arc=np.r_[0.,np.cumsum(lengths)]
    positions=[];errors=[]
    for point in x:
        options=[]
        for j in range(len(y)-1):
            delta=y[j+1]-y[j]
            fraction=float(np.clip((point-y[j])@delta/max(delta@delta,1e-15),0,1))
            options.append((float(np.sum((point-y[j]-fraction*delta)**2)),j,fraction))
        error,j,fraction=min(options)
        positions.append((arc[j]+fraction*lengths[j])/arc[-1]);errors.append(error)
    return np.asarray(positions),np.asarray(errors),float(arc[-1])


def main():
    began=time.time()
    index_path=HERE/'ensemble_index.json';index=json.loads(index_path.read_text())
    assert index['status']=='complete' and index['completed_fits']==4800,'Wait until production completes.'
    for path,digest in index['source_sha256'].items():
        assert sha(ROOT/path)==digest,path
    common_path=HERE.parent/'whole_manuscript_topology_20261015/common.py'
    spec=importlib.util.spec_from_file_location('weighted_audit_common',common_path)
    common=importlib.util.module_from_spec(spec);spec.loader.exec_module(common)
    data=common.load_data();x,_=common.embed(data['ds']['combined'])
    meta=data['meta'];labels=np.array([m['group'] for m in meta]);blocks=np.array(data['blocks'])
    unique=sorted(set(blocks))
    weights=np.array([1./np.sum(blocks==b) for b in blocks]);weights/=weights.sum()
    out=dict(status='running',new_growth_calls=2,replays=[],source_sha256={
        'ensemble_index.json':sha(index_path),'ensemble_independent_audit.py':sha(__file__)},
        scope='Two direct unmodified-package weighted growth replays; independent weights, anchors, chain orientation, polish, objective components and segment projection; no full models retained.')
    with threadpool_limits(limits=1):
        for key,replicate,start_index in [('paragraph207_point_l010',0,0),('paragraph207_sheet_l050',7,7)]:
            path=HERE/f'ensemble_{key}.json';saved=json.loads(path.read_text());design=saved['design']
            assert saved['status']=='complete' and sha(path)==index['result_sha256'][path.name]
            out['source_sha256'][path.name]=sha(path)
            reference=next(r for r in saved['runs'] if r['replicate']==replicate and r['start_index']==start_index)
            mask=np.asarray(design['masks'][replicate],dtype=bool)
            if design['sampler']=='point':
                independent_mask=np.random.RandomState(design['mask_seed']).random_sample((100,len(x)))[replicate]<=.9
            else:
                sheet_draws=np.random.RandomState(design['mask_seed']).random_sample((100,len(unique)))<=.9
                independent_mask=sheet_draws[replicate,[unique.index(b) for b in blocks]]
            assert np.array_equal(mask,independent_mask)
            assert np.max(abs(weights-np.asarray(design['original_weights'])))<1e-15
            w=weights[mask]/weights[mask].sum();lam=design['penalty']
            assert hashlib.sha256(w.tobytes()).hexdigest()==reference['normalized_weight_sha256']
            # Production's public grow/refine functions normalize their input
            # once more; reproduce that arithmetic for an exact package replay.
            package_weights=w/w.sum()
            initial=np.array([weights[labels==g]@x[labels==g]/weights[labels==g].sum()
                              for g in design['starts'][start_index]['anchor_order']])
            assert hashlib.sha256(initial.tobytes()).hexdigest()==design['starts'][start_index]['initial_sha256']
            with contextlib.redirect_stdout(io.StringIO()):
                native=elpigraph.computeElasticPrincipalCurve(x[mask],NumNodes=28,Lambda=lam,Mu=.1,
                    InitNodePositions=initial.copy(),InitEdges=np.array([[j,j+1] for j in range(len(initial)-1)]),
                    Do_PCA=False,CenterData=False,PointWeights=package_weights[:,None],FinalEnergy='Base',alpha=0,beta=0,
                    MaxNumberOfIterations=100,eps=1e-5,n_cores=1,nReps=1,ProbPoint=1,EmbPointProb=1,
                    SampleIC=False,AvoidResampling=True,ComputeMSEP=False,verbose=False,ShowTimer=False,
                    DisplayWarnings=False)[0]
            y=orient(native['NodePositions'],native['Edges'][0].astype(int),x,labels)
            growth_error=abs(energy(x[mask],package_weights,y,lam)[0][0]-reference['growth_objective'])
            # Reproduce the production polish for exact full-node hash recovery,
            # then independently solve and compare every hidden coordinate.
            reproduced=core.refine_curve(x[mask],w,y,lam=lam)
            reproduced_hash=hashlib.sha256(reproduced['nodes'].tobytes()).hexdigest()
            assert reproduced_hash==reference['full_node_sha256']
            y,iterations=polish(x[mask],package_weights,y,lam)
            independent_node_error=float(np.max(abs(y-reproduced['nodes'])))
            values,_=energy(x[mask],package_weights,y,lam)
            target=np.array([reference[k] for k in ['objective','node_mse','edge_energy','star_energy']])
            component_error=float(np.max(abs(values-target)))
            positions,errors,length=project(x,y)
            position_error=float(np.max(abs(positions-np.asarray(reference['positions']))))
            node3_error=float(np.max(abs(y[:,:3]-np.asarray(reference['nodes3']))))
            full_error=abs(energy(x,weights,y,lam)[0][0]-reference['full_reference_objective'])
            assert max(growth_error,component_error,full_error)<1e-11
            assert independent_node_error<1e-11
            assert position_error<1e-10 and node3_error<5.1e-11
            assert abs(length-reference['path_length'])<1e-10
            assert abs(w@errors[mask]-reference['training_segment_mse'])<1e-11
            out['replays'].append(dict(setting=key,replicate=replicate,start_index=start_index,
                retained=int(mask.sum()),growth_objective_error=float(growth_error),
                objective_component_max_error=component_error,position_max_error=position_error,
                saved_rounded_nodes3_max_error=node3_error,full_reference_objective_error=float(full_error),
                reproduced_full_node_sha256=reproduced_hash,full_node_hash_matches=True,
                independent_full_node_max_error=independent_node_error,
                independent_polish_iterations=iterations))
    out.update(status='passed',elapsed_seconds=time.time()-began)
    (HERE/'ensemble_independent_audit.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n')
    print(json.dumps(out))


if __name__=='__main__':
    main()
