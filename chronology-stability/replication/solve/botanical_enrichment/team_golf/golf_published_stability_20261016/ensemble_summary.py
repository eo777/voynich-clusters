"""Retained-row frequency and jointly retained pair summaries; no model fits."""
from pathlib import Path
import collections
import hashlib
import itertools
import json
import numpy as np

HERE=Path(__file__).resolve().parent
GROUPS=['1a','1b','1.5a','2a','3a','3b','3c','4a','4c','5a','6','']
TOL=1e-7
BINS=np.linspace(0,1,101)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path,value):
    def default(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    payload=json.dumps(value,default=default,separators=(',',':'),allow_nan=False)
    assert len(payload.encode())<12_000_000
    path.write_text(payload,encoding='utf-8')


def describe_runs(obj,runs):
    design=obj['design'];n=len(design['meta']);masks=np.asarray(design['masks'],dtype=bool)
    blocks=np.array(design['blocks']);positions=np.array([r['positions'] for r in runs])
    labels=np.array([m['group'] for m in design['meta']])
    groups=[g for g in GROUPS if g in set(labels)]+sorted(set(labels)-set(GROUPS))
    retained=masks[[r['replicate'] for r in runs]] if runs else np.empty((0,n),dtype=bool)
    row=[]
    for i in range(n):
        values=positions[retained[:,i],i] if runs else np.array([])
        histogram=np.histogram(values,BINS)[0]
        assert histogram.sum()==len(values), 'Projected positions outside frequency bins'
        row.append(dict(id=design['meta'][i]['id'],retained_fits=len(values),
            retained_mask_draws=len({r['replicate'] for r,m in zip(runs,retained[:,i]) if m}),
            quantiles=np.quantile(values,[.05,.25,.5,.75,.95]) if len(values) else None,
            histogram=histogram))
    # Compact upper-triangle arrays preserve denominator, signs, ties separately.
    ii,jj=np.triu_indices(n,1)
    cross=blocks[ii]!=blocks[jj]
    before=np.zeros(len(ii),dtype=int);tie=before.copy();after=before.copy();denom=before.copy()
    for r,mask in zip(runs,retained):
        t=np.asarray(r['positions']);eligible=mask[ii]&mask[jj];d=t[ii]-t[jj]
        before+=(d < -TOL)&eligible;after+=(d > TOL)&eligible;tie+=(abs(d)<=TOL)&eligible;denom+=eligible
    assert np.array_equal(before+tie+after,denom)
    group_summaries={}
    for readout in ['all','retained']:
        field='group_readout_'+readout
        median_readouts=[]
        for ri,r in enumerate(runs):
            selected=np.ones(n,dtype=bool) if readout=='all' else retained[ri]
            median_readouts.append({g:float(np.median(positions[ri,(labels==g)&selected]))
                if np.any((labels==g)&selected) else None for g in groups})
        bygroup={}
        for g in groups:
            values=[m[g] for m in median_readouts if m[g] is not None]
            bygroup[g]=dict(count=len(values),histogram=np.histogram(values,BINS)[0],
                quantiles=np.quantile(values,[.05,.25,.5,.75,.95]) if values else None)
        pair_counts={}
        for a,b in itertools.combinations(groups,2):
            key=a+'|'+b;signs=[]
            for medians in median_readouts:
                aa=medians[a];bb=medians[b]
                signs.append(None if aa is None or bb is None else -1 if aa-bb < -TOL else 1 if aa-bb > TOL else 0)
            pair_counts[key]=dict(before=signs.count(-1),tie=signs.count(0),after=signs.count(1),missing=signs.count(None))
        joint=collections.Counter(tuple(r[field]['contrasts'][key]['sign'] for key in ['1a|1b','3a|3b','4a|4c']) for r in runs)
        group_summaries[readout]=dict(groups=bygroup,pair_counts=pair_counts,
            joint_signatures=[dict(signs=list(k),count=v) for k,v in sorted(joint.items(),key=lambda kv:str(kv[0]))])
    def selected_pairs(mask):
        return dict(i=ii[mask],j=jj[mask],before=before[mask],tie=tie[mask],after=after[mask],denominator=denom[mask])
    return dict(fits=len(runs),mask_draws=len({r['replicate'] for r in runs}),rows=row,
        cross_sheet_pair_counts=selected_pairs(cross),within_sheet_pair_counts=selected_pairs(~cross),
        group_readouts=group_summaries)


def main():
    index=json.loads((HERE/'ensemble_index.json').read_text())
    assert index['status']=='complete'
    output=dict(status='complete',histogram_edges=BINS,tie_tolerance=TOL,
        interpretation='Empirical sensitivity frequencies; 100 shared masks and eight declared starts per setting, not 800 independent samples and not historical probabilities.',
        scopes={'all_starts':'All eight starts, retained-row / jointly retained cross-sheet pairs.',
                'final_1a_first':'Subset whose all-reference group medians put 1a before 1b; this is an outcome filter, not a chronology probability.'},
        source_sha256={'ensemble_summary.py':sha(Path(__file__)),'ensemble_index.json':sha(HERE/'ensemble_index.json')},settings=[])
    for filename in index['setting_files']:
        path=HERE/filename;obj=json.loads(path.read_text());assert obj['status']=='complete' and len(obj['runs'])==800
        output['source_sha256'][filename]=sha(path)
        runs=obj['runs'];d=obj['design']
        scoped={'all_starts':runs,'final_1a_first':[r for r in runs if r['group_readout_all']['contrasts']['1a|1b']['sign']==-1]}
        detailed={scope:describe_runs(obj,rs) for scope,rs in scoped.items()}
        summary_file=filename.replace('ensemble_','ensemble_frequency_')
        dump(HERE/summary_file,dict(cohort=d['cohort'],sampler=d['sampler'],penalty=d['penalty'],
            meta=d['meta'],histogram_edges=BINS,tie_tolerance=TOL,source_sha256=sha(path),scopes=detailed))
        retained=np.asarray(d['masks'],dtype=bool)
        settingsummary=dict(cohort=d['cohort'],sampler=d['sampler'],penalty=d['penalty'],
            source_file=filename,frequency_file=summary_file,fit_count=len(runs),
            converged=sum(r['converged'] for r in runs),
            retained_per_mask_range=[int(retained.sum(1).min()),int(retained.sum(1).max())],
            retained_per_mask_mean=float(retained.sum(1).mean()),
            maximum_energy_increase=max(r['maximum_energy_increase'] for r in runs),
            max_objective_component_discrepancy=max(abs(r['objective']-r['node_mse']-r['edge_energy']-r['star_energy']) for r in runs),
            scopes={scope:dict(fits=value['fits'],group_readouts=value['group_readouts']) for scope,value in detailed.items()},
            per_start=[])
        for si,start in enumerate(d['starts']):
            rr=[r for r in runs if r['start_index']==si]
            settingsummary['per_start'].append(dict(start=start,fit_count=len(rr),
                final_joint_signatures=dict(collections.Counter(','.join(str(r['group_readout_all']['contrasts'][key]['sign']) for key in ['1a|1b','3a|3b','4a|4c']) for r in rr))))
        output['settings'].append(settingsummary)
    dump(HERE/'ensemble_summary.json',output)
    print(json.dumps([{k:s[k] for k in ['cohort','sampler','penalty','fit_count','converged','retained_per_mask_range']} for s in output['settings']],indent=2))


if __name__=='__main__':main()
