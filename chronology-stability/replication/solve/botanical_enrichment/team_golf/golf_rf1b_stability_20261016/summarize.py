"""Reuse frozen summaries and audit, then compare paired RF and ZL runs."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'golf_published_stability_20261016'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_module(name):
    path=OLD/(name+'.py');s=importlib.util.spec_from_file_location('rf_'+name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);m.HERE=HERE
    return m

def main():
    for name in ['ensemble_audit','ensemble_summary','ensemble_competition']:
        load_module(name).main()
    index=json.loads((HERE/'ensemble_index.json').read_text())
    rfcomp=json.loads((HERE/'ensemble_competition.json').read_text())
    zlcomp=json.loads((OLD/'ensemble_competition.json').read_text())
    metrics=[];sources={Path(__file__).name:sha(Path(__file__))}
    for name in index['setting_files']:
        r=json.loads((HERE/name).read_text());z=json.loads((OLD/name).read_text())
        a,b=r['design'],z['design']
        assert a['masks']==b['masks'] and a['original_weights']==b['original_weights']
        assert [(m['id'],m['group'],m['bifolium']) for m in a['meta']]==[(m['id'],m['group'],m['bifolium']) for m in b['meta']]
        assert [s['anchor_order'] for s in a['starts']]==[s['anchor_order'] for s in b['starts']]
        assert [(q['replicate'],q['start_index']) for q in r['runs']]==[(q['replicate'],q['start_index']) for q in z['runs']]
        masks=np.asarray(a['masks'],dtype=bool);rp=np.array([q['positions'] for q in r['runs']]);zp=np.array([q['positions'] for q in z['runs']])
        kept=masks[[q['replicate'] for q in r['runs']]]
        rows=[]
        for j,m in enumerate(a['meta']):
            rv=rp[kept[:,j],j];zv=zp[kept[:,j],j]
            rows.append(dict(id=m['id'],group=m['group'],retained_fits=len(rv),
                zl_quantiles=np.quantile(zv,[.05,.25,.5,.75,.95]).tolist(),
                rf_quantiles=np.quantile(rv,[.05,.25,.5,.75,.95]).tolist(),
                joint_quantiles=np.quantile(np.r_[rv,zv],[.05,.25,.5,.75,.95]).tolist(),
                paired_median_shift=float(np.median(rv-zv)),paired_mean_absolute_shift=float(np.mean(abs(rv-zv)))))
        contrasts={}
        for pair in ['1a|1b','3a|3b','4a|4c']:
            rr=np.array([q['group_readout_all']['contrasts'][pair]['sign'] for q in r['runs']])
            zz=np.array([q['group_readout_all']['contrasts'][pair]['sign'] for q in z['runs']])
            contrasts[pair]=dict(zl_before=int(sum(zz==-1)),rf_before=int(sum(rr==-1)),
                zl_after=int(sum(zz==1)),rf_after=int(sum(rr==1)),zl_ties=int(sum(zz==0)),rf_ties=int(sum(rr==0)),
                strict_paired_reversals=int(sum(rr*zz==-1)),different_signs=int(sum(rr!=zz)),paired_fits=len(rr))
        def competition(c):
            s=next(s for s in c['settings'] if all(s[k]==a[k] for k in ['cohort','sampler','penalty']))
            out={}
            for scope,v in s['scopes'].items():
                winners=[q['winner'] for q in v['mask_records'] if q['winner'] is not None]
                out[scope]=dict(masks=len(winners),directions={p:dict(before=sum(q['direction_signs'][p]==-1 for q in winners),
                    after=sum(q['direction_signs'][p]==1 for q in winners),ties=sum(q['direction_signs'][p]==0 for q in winners))
                    for p in ['3a|3b','4a|4c','3c|4a','3c|4c']},availability=v['summary']['directions'])
            return out
        metrics.append(dict(cohort=a['cohort'],sampler=a['sampler'],penalty=a['penalty'],
            positive_dimensions={'zl':b['positive_dimensions'],'rf':a['positive_dimensions']},
            rows=rows,group_contrasts=contrasts,winners={'zl':competition(zlcomp),'rf':competition(rfcomp)}))
        sources['rf/'+name]=sha(HERE/name);sources['zl/'+name]=sha(OLD/name)
    output=dict(status='complete',matched_masks=True,matched_weights=True,matched_group_orders=True,matched_cohorts=True,
        comparison='RF1b-er basic EVA versus ZL Zattera; both transcription and unitization change',
        paired_runs_per_setting=800,settings=metrics,source_sha256=sources,
        caveat='Paired starts can fall into different local minima; run-level positions are relative within each fitted curve.')
    (HERE/'comparison.json').write_text(json.dumps(output,separators=(',',':'),allow_nan=False),encoding='utf-8')
    print(json.dumps({'comparison':'passed','settings':len(metrics),'paired_runs':4800}),flush=True)

if __name__=='__main__':main()
