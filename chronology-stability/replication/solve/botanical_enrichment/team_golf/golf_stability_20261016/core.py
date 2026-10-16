"""Shared exact data and elastic-curve operations for stability experiments.

Growth uses unmodified ElPiGraph. Fixed-graph refinement exposes the same
alternating quadratic solve with assignment and energy diagnostics, checked
against ElPiGraph's primitive embedding. No smoothing of the JS eligibility rule.
"""
from pathlib import Path
import contextlib
import hashlib
import importlib.util
import io
import itertools
import json
import sys
import warnings

import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent
sys.path.insert(0, str(HERE / '_deps'))
import elpigraph
from elpigraph.src.core import MakeUniformElasticMatrix, ComputeSpringLaplacianMatrix
from elpigraph.src.distutils import ComputePrimitiveGraphElasticEnergy


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


COMMON = module('stability_frozen_data', TEAM/'whole_manuscript_topology_20261015/common.py')
MODELS = module('stability_frozen_projection', TEAM/'whole_manuscript_topology_20261015/models.py')
SEED = 20261016
TOL = 1e-10
ORDER = ['1a','1b','1.5a','2a','3a','3b','3c','4a','4c']
CONTRASTS = [('1a','1b'),('3a','3b'),('4a','4c')]


def load_data():
    return COMMON.load_data()


def load_count_data():
    data = load_data()
    loader = module('stability_raw_counts', TEAM/'b_horseshoe_20261015/drivers.py')
    alphabet = [r['unit'] for r in loader.REPRO.rows('alphabet.csv')]
    meta = data['pool_meta']
    pages = data['meta']
    target_ids = ['f65v','f58r','f58v']
    selected = [next(m for m in pages if m['id']==uid) for uid in target_ids]
    ids = [m['profile_id'] for m in meta+selected]
    arrays = {mode:loader.read_counts(mode,ids,alphabet) for mode in ['word','line']}
    pools = {mode:x[:len(meta)].copy() for mode,x in arrays.items()}
    partial = {mode:x[len(meta)].copy() for mode,x in arrays.items()}
    added = {mode:x[len(meta)+1:].sum(0) for mode,x in arrays.items()}
    hi = next(i for i,m in enumerate(meta) if m['bifolium']=='H:58-65')
    for mode in pools:
        np.testing.assert_array_equal(pools[mode][hi], partial[mode]+added[mode])
    got = distance_from_counts(pools)
    discrepancy = float(abs(got-data['pool_ds']['combined']).max())
    assert discrepancy < 1e-12
    paths = [loader.DATA/'alphabet.csv',loader.DATA/'profiles.csv',loader.DATA/'reproduce.py']
    paths += [loader.DATA/f'transitions-{mode}-counts.csv' for mode in pools]
    return dict(meta=meta,pool_counts=pools,h_partial=partial,h_added=added,h_index=hi,
        alphabet=alphabet,source_hashes=COMMON.hash_sources(paths),distance_recovery_max_error=discrepancy)


def distance_from_counts(counts):
    repro = module('stability_exact_js', COMMON.DATA/'reproduce.py')
    return sum(repro.weighted_js(np.asarray(x,dtype=float))[0] for x in counts.values())/len(counts)


def align_embedding(x, reference):
    width = max(x.shape[1],reference.shape[1])
    a = np.pad(x,((0,0),(0,width-x.shape[1])))
    b = np.pad(reference,((0,0),(0,width-reference.shape[1])))
    rotation,_ = orthogonal_procrustes(a,b)
    return a@rotation,rotation


def normalized(weights):
    w = np.asarray(weights,dtype=float).reshape(-1)
    assert np.all(w>=0) and w.sum()>0
    return w/w.sum()


def make_starts(x,meta,weights):
    w=normalized(weights); groups=np.array([m['group'] for m in meta])
    centers={}
    for g in ORDER:
        mask=groups==g
        assert (w[mask]>0).any(), 'No positive-weight anchor: '+g
        centers[g]=w[mask]@x[mask]/w[mask].sum()
    starts={}
    for e,h,t in itertools.product([0,1],repeat=3):
        order=ORDER.copy()
        for flag,i in [(e,0),(h,4),(t,7)]:
            if flag:order[i:i+2]=order[i:i+2][::-1]
        starts[f'e{e}_h{h}_t{t}']=np.array([centers[g] for g in order])
    mean=w@x
    _,_,axes=np.linalg.svd(np.sqrt(w[:,None])*(x-mean),full_matrices=False)
    axis=axes[0]
    if axis[np.argmax(abs(axis))]<0:axis=-axis
    extent=np.sqrt(w@((x-mean)@axis)**2)
    starts['PCA']=mean+np.array([-1.,0.,1.])[:,None]*extent*axis
    return starts


def energy(x,w,nodes,edges,lam,mu):
    matrix=MakeUniformElasticMatrix(np.asarray(edges,dtype=int),lam,mu)
    d=cdist(x,nodes,'sqeuclidean')
    assignment=d.argmin(1)
    parts=ComputePrimitiveGraphElasticEnergy(nodes,matrix,d.min(1)[:,None],w[:,None])
    return dict(objective=float(parts[0]),training_objective=float(parts[0]),
        node_mse=float(parts[1]),edge_energy=float(parts[2]),star_energy=float(parts[3]),
        assignments=assignment)


def grow_curve(x,weights,initial,lam,mu=.1,n_nodes=28):
    w=normalized(weights); initial=np.asarray(initial,dtype=float)
    with contextlib.redirect_stdout(io.StringIO()),warnings.catch_warnings(record=True) as caught:
        output=elpigraph.computeElasticPrincipalCurve(np.asarray(x,dtype=float),NumNodes=n_nodes,
            Lambda=lam,Mu=mu,InitNodePositions=initial.copy(),
            InitEdges=np.array([[i,i+1] for i in range(len(initial)-1)]),
            Do_PCA=False,CenterData=False,PointWeights=w[:,None],FinalEnergy='Base',
            alpha=0,beta=0,MaxNumberOfIterations=100,eps=1e-5,n_cores=1,nReps=1,
            ProbPoint=1,EmbPointProb=1,SampleIC=False,AvoidResampling=True,
            ComputeMSEP=False,verbose=False,ShowTimer=False,DisplayWarnings=False)[0]
    nodes=output['NodePositions'];edges=output['Edges'][0].astype(int)
    return dict(nodes=nodes,edges=edges,lam=lam,mu=mu,warning_count=len(caught),
        **energy(x,w,nodes,edges,lam,mu))


def refine_curve(x,weights,nodes,edges=None,lam=.01,mu=.1,tol=1e-9,maxiter=1000):
    x=np.asarray(x,dtype=float);w=normalized(weights);nodes=np.asarray(nodes,dtype=float).copy()
    if edges is None:edges=np.array([[i,i+1] for i in range(len(nodes)-1)])
    edges=np.asarray(edges,dtype=int)
    elastic=MakeUniformElasticMatrix(edges,lam,mu)
    lap=ComputeSpringLaplacianMatrix(elastic)
    history=[];changes=[];converged=False;movement=float('inf')
    start=energy(x,w,nodes,edges,lam,mu)
    history.append(start['objective']);assignment=start['assignments']
    for iteration in range(maxiter):
        mass=np.bincount(assignment,weights=w,minlength=len(nodes))
        sums=np.zeros_like(nodes)
        np.add.at(sums,assignment,w[:,None]*x)
        updated=np.linalg.solve(np.diag(mass)+lap,sums)
        current=energy(x,w,updated,edges,lam,mu)
        movement=float(np.max(np.linalg.norm(updated-nodes,axis=1)))
        changed=int(np.sum((current['assignments']!=assignment)&(w>0)))
        changes.append(changed);history.append(current['objective'])
        assert history[-1]<=history[-2]+1e-11, (history[-2],history[-1])
        nodes=updated;assignment=current['assignments']
        if changed==0 and movement<=tol:
            converged=True
            break
    mass=np.bincount(assignment,weights=w,minlength=len(nodes))
    curvature=np.linalg.eigvalsh(2*(np.diag(mass)+lap))
    return dict(nodes=nodes,edges=edges,lam=lam,mu=mu,converged=converged,iterations=iteration+1,
        final_node_max_move=movement,energy_trace=history,assignment_changes=changes,
        minimum_fixed_assignment_hessian_eigenvalue=float(curvature[0]),
        maximum_fixed_assignment_hessian_eigenvalue=float(curvature[-1]),**current)


def chain_order(edges,n):
    neighbors=[[] for _ in range(n)]
    for a,b in edges:neighbors[a].append(int(b));neighbors[b].append(int(a))
    assert max(map(len,neighbors))<=2
    leaves=[i for i,v in enumerate(neighbors) if len(v)==1]
    assert len(leaves)==2
    route=[min(leaves)]
    while len(route)<n:
        route.append(next(v for v in neighbors[route[-1]] if len(route)==1 or v!=route[-2]))
    return np.array(route)


def orient_curve(model,x,meta):
    out=dict(model);nodes=np.asarray(model['nodes'])[chain_order(model['edges'],len(model['nodes']))]
    groups=np.array([m['group'] for m in meta])
    early=x[np.isin(groups,['1a','1b'])].mean(0)
    late=x[np.isin(groups,['4a','4c'])].mean(0)
    direct=np.sum((nodes[0]-early)**2)+np.sum((nodes[-1]-late)**2)
    reverse=np.sum((nodes[-1]-early)**2)+np.sum((nodes[0]-late)**2)
    if reverse<direct:nodes=nodes[::-1]
    out['nodes']=nodes.copy();out['edges']=np.array([[i,i+1] for i in range(len(nodes)-1)])
    # Assignments, if present, use the old node numbering and are recomputed by refinement.
    out.pop('assignments',None)
    return out


def project_curve(x,model):
    nodes=np.asarray(model['nodes']);edges=np.asarray(model['edges'])
    assert np.array_equal(edges,np.array([[i,i+1] for i in range(len(nodes)-1)])), 'Canonical ordered chain required'
    errors,edge,fraction=MODELS.project(np.asarray(x),dict(family='curve',nodes=nodes,edges=edges))
    lengths=np.linalg.norm(np.diff(nodes,axis=0),axis=1)
    arc=np.r_[0,np.cumsum(lengths)]
    positions=(arc[edge]+fraction*lengths[edge])/max(arc[-1],1e-15)
    return dict(positions=positions,errors=errors,edge=edge,fraction=fraction,arc_length=float(arc[-1]))


def summarize_curve(model,x,meta,weights=None):
    result=project_curve(x,model);positions=result['positions']
    groups=np.array([m['group'] for m in meta])
    medians={g:float(np.median(positions[groups==g])) for g in sorted(set(groups))}
    contrasts={};signature=[]
    for a,b in CONTRASTS:
        if a not in medians or b not in medians:continue
        value=medians[a]-medians[b]
        sign=0 if abs(value)<=TOL else int(np.sign(value))
        contrasts[a+'_vs_'+b]=dict(difference=value,sign=sign,relation='tied' if sign==0 else a+'_before_'+b if sign<0 else b+'_before_'+a)
        signature.append(str(sign))
    w=np.ones(len(x))/len(x) if weights is None else normalized(weights)
    result.update(group_medians=medians,contrasts=contrasts,signature=','.join(signature),
        segment_error=float(w@result['errors']),objective=model['objective'])
    return result


def jsonable(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,Path):return str(value)
    if isinstance(value,dict):return {str(k):jsonable(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [jsonable(v) for v in value]
    return value


def save(name,data,max_bytes=12_000_000):
    payload=json.dumps(jsonable(data),allow_nan=False,separators=(',',':'))+'\n'
    assert len(payload.encode('utf-8'))<max_bytes
    (HERE/name).write_text(payload,encoding='utf-8')


def source_hashes(paths=()):
    return COMMON.hash_sources([Path(p).resolve() for p in [Path(__file__),TEAM/'whole_manuscript_topology_20261015/common.py',
        TEAM/'whole_manuscript_topology_20261015/models.py',*paths]])
