"""Conditional H-location maps using saved ElPiGraph chains and explicit refinement.

This changes geometric coordinates, not counts or recovered missing text.
All grid seeds are the same saved baselines, independent of scan order.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import time

import numpy as np
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STABILITY = HERE.parent/'golf_stability_20261016'
COMMON_PATH = HERE.parent/'whole_manuscript_topology_20261015/common.py'
spec = importlib.util.spec_from_file_location('influence_frozen_common', COMMON_PATH)
COMMON = importlib.util.module_from_spec(spec)
spec.loader.exec_module(COMMON)
PAIRS = [('1a','1b'),('3a','3b'),('4a','4c')]
TOL = 1e-8


def matrices(n, lam, mu=.1):
    edge = np.diff(np.eye(n), axis=0)
    star = np.zeros((n-2,n))
    for i in range(1,n-1):
        star[i-1,i] = 1
        star[i-1,[i-1,i+1]] = -.5
    return lam*edge.T@edge+mu*star.T@star, edge, star


def components(x, nodes, edge, star, lam, mu=.1):
    distances = cdist(x, nodes, 'sqeuclidean')
    node = float(distances.min(axis=1).mean())
    edges = float(lam*np.sum((edge@nodes)**2))
    stars = float(mu*np.sum((star@nodes)**2))
    return node+edges+stars, node, edges, stars, distances.argmin(axis=1)


def project(x, nodes):
    a = nodes[:-1]
    v = np.diff(nodes, axis=0)
    l2 = np.sum(v*v,axis=1)
    dot = x@v.T-np.sum(a*v,axis=1)
    frac = np.clip(dot/np.maximum(l2,1e-15),0,1)
    errors = np.maximum(np.sum(x*x,axis=1)[:,None]+np.sum(a*a,axis=1)[None,:]
        -2*x@a.T-2*frac*dot+frac*frac*l2,0)
    nearest = errors.argmin(axis=1)
    lengths = np.sqrt(l2)
    arc = np.r_[0,np.cumsum(lengths)]
    positions = (arc[nearest]+frac[np.arange(len(x)),nearest]*lengths[nearest])/arc[-1]
    return positions, errors[np.arange(len(x)),nearest]


def refine(x, initial, lam, tol=1e-9):
    nodes = initial.copy()
    lap, edge, star = matrices(len(nodes),lam)
    energy,_,_,_,assignment = components(x,nodes,edge,star,lam)
    max_increase = 0.
    for iteration in range(1000):
        mass = np.bincount(assignment,minlength=len(nodes))/len(x)
        sums = np.zeros_like(nodes)
        np.add.at(sums,assignment,x/len(x))
        updated = np.linalg.solve(np.diag(mass)+lap,sums)
        new_energy,_,_,_,new_assignment = components(x,updated,edge,star,lam)
        max_increase = max(max_increase,new_energy-energy)
        assert new_energy<=energy+1e-11
        movement = np.linalg.norm(updated-nodes,axis=1).max()
        same = np.array_equal(assignment,new_assignment)
        nodes,assignment,energy = updated,new_assignment,new_energy
        if same and movement<=tol:
            return nodes,energy,iteration+1,max_increase
    raise AssertionError('Fixed-chain refinement did not converge')


def signs(positions, meta):
    groups = np.array([m['group'] for m in meta])
    contrasts = np.array([np.median(positions[groups==a])-np.median(positions[groups==b]) for a,b in PAIRS])
    return np.where(contrasts < -TOL,-1,np.where(contrasts > TOL,1,0)).tolist()


def save(name,value):
    (HERE/name).write_text(json.dumps(value,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8')


def main():
    started = time.time()
    frozen = COMMON.load_data()
    meta = frozen['pool_meta']
    x,eigenvalues = COMMON.embed(frozen['pool_ds']['combined'])
    hi = next(i for i,m in enumerate(meta) if m['bifolium']=='H:58-65')
    basins = json.loads((STABILITY/'basins.json').read_text())
    settings = [s for s in basins['settings'] if s['cohort']=='pool50']
    palette = json.loads((HERE.parent/'nonlinear_sequence_20261010/plot_data.json').read_text())['palette']
    palette['']='#8b8f97'
    sources = [Path(__file__).resolve(),HERE/'DESIGN.md',COMMON_PATH,STABILITY/'basins.json']
    output = dict(status='running', design=dict(grid_size=25,planes=[[0,1],[0,2]],
        initializations=[r['start'] for r in settings[0]['baseline_runs']],
        contrast_pairs=PAIRS, tie_tolerance=TOL, scope='Hypothetical coordinate slices with H hidden coordinates and other49points fixed; not inferred text profiles or universal2Dregions.',
        fitting='All positive dimensions; fixed 28-node topology; independent refits from saved ElPiGraph baselines.',
        winner_scopes=['all','final_1a_first']),
        meta=meta, h_index=hi, reference_coords3=np.round(x[:,:3],10).tolist(),
        full_dimension=x.shape[1],positive_variance_shares=(eigenvalues[:3]/eigenvalues[eigenvalues>1e-10].sum()).tolist(),
        palette=palette, source_sha256={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        baseline_checks=[],settings=[])
    totalfits=0;maxiter=0;maxincrease=0.
    with threadpool_limits(limits=1):
        for setting in settings:
            lam=setting['lambda']; bases=setting['baseline_runs']
            initial=[np.asarray(r['nodes']) for r in bases]
            for base,nodes in zip(bases,initial):
                lap,edge,star=matrices(len(nodes),lam)
                energy,node,ee,se,assignment=components(x,nodes,edge,star,lam)
                errors=[abs(energy-base['objective']),abs(node-base['node_mse']),
                        abs(ee-base['edge_energy']),abs(se-base['star_energy'])]
                positions,_=project(x,nodes)
                position_error=float(abs(positions-np.asarray(base['positions'])).max())
                recovered,_,_,_=refine(x,nodes,lam)
                recovery_error=float(abs(recovered-nodes).max())
                assert max(errors)<1e-12 and position_error<1e-12 and recovery_error<1e-10
                output['baseline_checks'].append(dict(penalty=lam,start=base['start'],
                    component_max_error=max(errors),position_max_error=position_error,recovery_max_error=recovery_error))
            for plane in [[0,1],[0,2]]:
                axes=[]
                for axis in plane:
                    lo,hi_axis=float(x[:,axis].min()),float(x[:,axis].max())
                    width=hi_axis-lo
                    grid=np.linspace(lo-.1*width,hi_axis+.1*width,25)
                    grid[np.argmin(abs(grid-x[hi,axis]))]=x[hi,axis]
                    assert (np.diff(grid)>0).all()
                    axes.append(grid)
                row=dict(penalty=lam,plane=plane,xs=axes[0].tolist(),ys=axes[1].tolist(),cells=[])
                for iy,y in enumerate(axes[1]):
                    for ix,xx in enumerate(axes[0]):
                        perturbed=x.copy();perturbed[hi,plane]=[xx,y]
                        models=[];objectives=[];relations=[];iterations=[]
                        for nodes in initial:
                            fitted,energy,it,increase=refine(perturbed,nodes,lam)
                            position,_=project(perturbed,fitted)
                            models.append(fitted);objectives.append(energy)
                            relations.append(signs(position,meta));iterations.append(it)
                            totalfits+=1;maxiter=max(maxiter,it);maxincrease=max(maxincrease,increase)
                        best=int(np.argmin(objectives))
                        early=[i for i,s in enumerate(relations) if s[0]==-1]
                        bestearly=min(early,key=lambda i:objectives[i]) if early else None
                        row['cells'].append(dict(objectives=objectives,signs=relations,best_all=best,best_1a_first=bestearly,
                            winner_nodes3=np.round(models[best][:,:3],6).tolist(),
                            winner_1a_nodes3=np.round(models[bestearly][:,:3],6).tolist() if bestearly is not None else None))
                    print('MAP',lam,plane,'row',iy+1,'fits',totalfits,'elapsed',round(time.time()-started,1),flush=True)
                baselineix=int(np.argmin(abs(axes[0]-x[hi,plane[0]])))
                baselineiy=int(np.argmin(abs(axes[1]-x[hi,plane[1]])))
                baselinecell=row['cells'][baselineiy*25+baselineix]
                assert max(abs(a-b['objective']) for a,b in zip(baselinecell['objectives'],bases))<1e-12
                row['baseline_cell_index']=baselineiy*25+baselineix
                output['settings'].append(row)
                output['progress_fit_count']=totalfits
                save('influence_grid_checkpoint.json',output)
    assert totalfits==22500
    output.update(status='complete',fit_count=totalfits,all_converged=True,
        maximum_iterations=maxiter,maximum_energy_increase=maxincrease,elapsed_seconds=time.time()-started)
    save('influence_grid.json',output)
    print('COMPLETE',totalfits,'fits',round(time.time()-started,1),'seconds',flush=True)


if __name__=='__main__':
    main()
