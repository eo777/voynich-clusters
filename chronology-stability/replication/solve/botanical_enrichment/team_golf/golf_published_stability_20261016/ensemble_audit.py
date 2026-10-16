"""Check ensemble accounting, masks, weights, energies and source integrity; no fits."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    index=json.loads((HERE/'ensemble_index.json').read_text())
    assert index['status']=='complete' and index['completed_fits']==4800
    checked_sources={}
    for name,digest in index['source_sha256'].items():
        actual=sha(ROOT/name);assert actual==digest,name
        checked_sources[name]=actual
    sources={'ensemble_audit.py':sha(Path(__file__)),'ensemble_index.json':sha(HERE/'ensemble_index.json')}
    settings=[];loaded={};maximum_component_error=0.;maximum_weight_error=0.
    for name in index['setting_files']:
        path=HERE/name;obj=json.loads(path.read_text());d=obj['design'];runs=obj['runs']
        sources[name]=sha(path);assert sources[name]==index['result_sha256'][name]
        assert obj['status']=='complete' and len(runs)==800
        assert {(r['replicate'],r['start_index']) for r in runs}=={(i,s) for i in range(100) for s in range(8)}
        n=len(d['meta']);masks=np.array(d['masks'],dtype=bool);w=np.array(d['original_weights']);blocks=np.array(d['blocks'])
        assert masks.shape==(100,n) and len(d['starts'])==8 and abs(w.sum()-1)<1e-12
        assert len({tuple(s['anchor_order']) for s in d['starts']})==8
        if d['cohort']=='paragraph207' and d['sampler']=='point':
            expected=np.random.RandomState(d['mask_seed']).random_sample((100,n))<=.9
        else:
            sheets=np.array(sorted(set(blocks)))
            draws=np.random.RandomState(d['mask_seed']).random_sample((100,len(sheets)))<=.9
            expected=draws[:,[np.flatnonzero(sheets==b)[0] for b in blocks]]
        assert np.array_equal(masks,expected)
        per_sheet=np.array([w[blocks==b].sum() for b in set(blocks)])
        assert abs(per_sheet-1/len(per_sheet)).max()<1e-12
        for r in runs:
            mask=masks[r['replicate']];sub=w[mask]/w[mask].sum()
            assert len(r['positions'])==n and np.shape(r['nodes3'])==(28,3)
            assert r['retained']==mask.sum() and r['converged']
            assert np.all(np.isfinite(r['positions'])) and min(r['positions'])>=-1e-12 and max(r['positions'])<=1+1e-12
            assert r['normalized_weight_sha256']==hashlib.sha256(sub.tobytes()).hexdigest()
            maximum_weight_error=max(maximum_weight_error,abs(r['retained_weight_before_normalization']-w[mask].sum()))
            error=abs(r['objective']-sum(r[k] for k in ['node_mse','edge_energy','star_energy']))
            maximum_component_error=max(maximum_component_error,error)
            assert error<1e-12 and r['maximum_energy_increase']<1e-11
            assert r['objective']<=r['growth_objective']+1e-11
            for kind,selected in [('all',np.ones(n,dtype=bool)),('retained',mask)]:
                readout=r['group_readout_'+kind];t=np.array(r['positions']);labels=np.array([m['group'] for m in d['meta']])
                for group,value in readout['group_medians'].items():
                    values=t[(labels==group)&selected]
                    assert value is None if len(values)==0 else abs(value-float(np.median(values)))<1e-12
                for pair,value in readout['contrasts'].items():
                    a,b=pair.split('|');aa=readout['group_medians'][a];bb=readout['group_medians'][b]
                    if aa is None or bb is None:assert value['delta'] is None and value['sign'] is None
                    else:
                        delta=aa-bb;sign=-1 if delta < -d['tie_tolerance'] else 1 if delta > d['tie_tolerance'] else 0
                        assert abs(value['delta']-delta)<1e-12 and value['sign']==sign
        loaded[(d['cohort'],d['sampler'],d['penalty'])]=d
        settings.append(dict(file=name,count=len(runs),converged=len(runs),retained_min=int(masks.sum(1).min()),retained_max=int(masks.sum(1).max())))
    for cohort,sampler in [('paragraph207','point'),('paragraph207','sheet'),('pool50','point')]:
        a=loaded[(cohort,sampler,.01)];b=loaded[(cohort,sampler,.05)]
        assert a['masks']==b['masks'] and a['geometry_sha256']==b['geometry_sha256']
    for lam in [.01,.05]:
        pool=loaded[('pool50','point',lam)];page=loaded[('paragraph207','sheet',lam)]
        mapping={b:i for i,b in enumerate(pool['blocks'])}
        pm=np.array(pool['masks']);qm=np.array(page['masks'])
        assert np.array_equal(pm[:,[mapping[b] for b in page['blocks']]],qm)
    result=dict(status='passed',source_sha256=sources,verified_frozen_sources=checked_sources,settings=settings,
        native_mask_law='independent Bernoulli inclusion p=.9',matched_pool_sheet_masks=True,paired_lambda_masks=True,
        production_fit_count=4800,pilot_growth_calls=2,no_growth_calls_in_this_audit=True,
        maximum_objective_component_error=maximum_component_error,maximum_weight_sum_error=maximum_weight_error,
        interpretation='Frequency checks verify the declared ensemble, not confidence coverage or historical ordering.')
    (HERE/'ensemble_audit.json').write_text(json.dumps(result,separators=(',',':')),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_sha256','verified_frozen_sources']},indent=2))


if __name__=='__main__':main()
