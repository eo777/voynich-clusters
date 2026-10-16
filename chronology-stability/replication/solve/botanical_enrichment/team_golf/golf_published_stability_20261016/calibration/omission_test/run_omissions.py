"""Paired omissions of observed data and the existing continuous-control streams."""
from pathlib import Path
import sys
import json
import time
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import run_calibration as cal


def write(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def cases(c):
    ids=np.array(c['ids']);groups=np.array([m['group'] for m in c['meta']])
    sections=np.array([m['section'] for m in c['meta']])
    f58=np.isin(ids,['f58r','f58v'])
    nearby=np.isin(ids,['f67v2','f69r','f70r2'])
    return [
        ('full','Full data',np.zeros(len(ids),bool)),
        ('omit_f58','Remove f58r/v',f58),
        ('omit_focal_five','Remove the five-profile focal bin',f58|nearby),
        ('omit_f58_5a','Remove f58 and Golf 5a',f58|(groups=='5a')),
        ('omit_f58_5a_3c_4c','Remove f58 and Golf 5a, 3c, 4c',f58|np.isin(groups,['5a','3c','4c'])),
        ('omit_sections_A_C_S','Remove source sections A, C, S',np.isin(sections,['A','C','S'])),
        ('retain_f58_omit_other_three','Retain f58; remove its three bin companions',nearby)]


def prepare(c):
    plans=[];max_error=0.
    focal=float(np.mean(c['t'][np.isin(c['ids'],['f58r','f58v'])]))
    for ident,title,omit in cases(c):
        keep=~omit;t=c['t'][keep];x=c['x'][keep];ids=np.array(c['ids'])[keep]
        sc=cal.Scorer(t,c['steps'])
        memberships={}
        for steps in c['steps']:
            assignment,_,occupancy=cal.audit.bins(t,steps)
            # Focal position remains internal, so adding it does not move cut endpoints.
            hbin=int(cal.audit.bins(np.r_[t,focal],steps)[0][-1])
            memberships[str(steps)]=dict(focal_bin=hbin,occupancy=occupancy.tolist(),
                retained_bins=(np.flatnonzero(occupancy>3)+1).tolist(),
                members={str(b):ids[assignment==b].tolist() for b in range(1,steps+1)})
        for b in sc.labels:
            hb=memberships[str(b['steps'])]['focal_bin']
            b['focal_join']=bool(b['left']<=hb<=b['right'])
        observed=sc.score(x)
        ref=[]
        for steps in c['steps']:
            for scale in ['per_million','raw']:
                mat=x/x.sum(axis=1)[:,None]*1e6 if scale=='per_million' else x
                e=cal.audit.expected(mat.T,c['features'],ids.tolist(),t,steps)
                if e['sufficient']:ref.extend(e['adjacent'])
        error=float(np.max(np.abs(observed-ref)));max_error=max(error,max_error)
        assert error<1e-12
        plans.append(dict(id=ident,title=title,keep=keep,scorer=sc,observed=observed,
            removed=np.array(c['ids'])[omit].tolist(),retained=int(keep.sum()),membership=memberships,
            position_range=[float(t.min()),float(t.max())],focal_position=focal))
    return plans,max_error


def main():
    start=time.monotonic();outputs=[]
    for case_i,c in enumerate(cal.load()[:2]):
        plans,error=prepare(c)
        for model_i,model in enumerate(cal.CONFIG['models']):
            generator=cal.Generator(c['x'],c['t'],c['sheets'],.12)
            arrays=[]
            with cal.threadpool_limits(limits=1):
                for stage,n in enumerate([400,999,500]):
                    random=cal.rng(case_i,model_i,1,stage)
                    rows=[[] for _ in plans]
                    for _ in range(n):
                        x=generator.draw(model,random)
                        for j,p in enumerate(plans):rows[j].append(p['scorer'].score(x[p['keep']]))
                    arrays.append([np.array(r) for r in rows])
            fit,calibration,evaluation=arrays
            means=[r.mean(axis=0) for r in fit];sd=[r.std(axis=0,ddof=1) for r in fit]
            assert all(np.all(v>0) for v in sd)
            oldpath=HERE.parent/f"null_paragraph207_{c['mode']}_{model}_h0.12.json"
            old=json.loads(oldpath.read_text())
            assert np.array_equal(means[0],np.array([b['null_mean'] for b in old['boundaries']]))
            assert np.array_equal(sd[0],np.array([b['null_sd'] for b in old['boundaries']]))
            focal_masks=[np.array([b['focal_join'] for b in p['scorer'].labels]) for p in plans]
            cz=[(v-m)/s for v,m,s in zip(calibration,means,sd)]
            ez=[(v-m)/s for v,m,s in zip(evaluation,means,sd)]
            full_max=np.concatenate(cz,axis=1).max(axis=1)
            focal_max=np.concatenate([z[:,mask] for z,mask in zip(cz,focal_masks)],axis=1).max(axis=1)
            efull=np.concatenate(ez,axis=1).max(axis=1)
            efocal=np.concatenate([z[:,mask] for z,mask in zip(ez,focal_masks)],axis=1).max(axis=1)
            results=[]
            for j,p in enumerate(plans):
                bounds=[]
                oz=(p['observed']-means[j])/sd[j]
                for k,b in enumerate(p['scorer'].labels):
                    bounds.append(dict(b,observed=float(p['observed'][k]),null_mean=float(means[j][k]),
                        null_sd=float(sd[j][k]),z=float(oz[k]),local_tail=cal.tail(calibration[j][:,k],p['observed'][k]),
                        targeted_adjusted_tail=min(1.,2*cal.tail(focal_max,oz[k])) if b['focal_join'] else None,
                        whole_path_adjusted_tail=min(1.,2*cal.tail(full_max,oz[k])),
                        reference95=np.quantile(calibration[j][:,k],[.025,.975]).tolist()))
                result={k:p[k] for k in ['id','title','removed','retained','membership','position_range','focal_position']}
                result['boundaries']=bounds
                results.append(result)
            out=dict(status='complete',channel=c['mode'],cohort=c['cohort'],model=model,bandwidth=.12,
                continuous_draws=1899,source_calibration_sha256=cal.sha(oldpath),
                code_sha256=cal.sha(Path(__file__)),design_sha256=cal.sha(HERE/'DESIGN.md'),
                scorer_max_error=error,full_data_reference_reproduced_exactly=True,cases=results,
                false_alarm_targeted=cal.rate([2*cal.tail(focal_max,z)<=.05 for z in efocal]),
                false_alarm_whole_path=cal.rate([2*cal.tail(full_max,z)<=.05 for z in efull]),
                reference_targeted_max=focal_max.tolist(),reference_whole_path_max=full_max.tolist())
            name=f"results_{c['mode']}_{model}.json";write(name,out);outputs.append(name)
            primary=[(r['id'],[(b['left'],b['right'],round(b['observed'],3),b['local_tail'],b['targeted_adjusted_tail'])
                for b in r['boundaries'] if b['focal_join'] and b['steps']==12 and b['scale']=='per_million']) for r in results]
            print(json.dumps(dict(file=name,seconds=round(time.monotonic()-start),H12=primary)),flush=True)
    write('audit.json',dict(status='passed',files={name:cal.sha(HERE/name) for name in outputs},
        continuous_draws=4*1899,observed_cases=7*2,models=2,
        source_R_score_verification='Every observed omission matches the independently R-validated calculation to <1e-12.',
        reproduction='All four full-data reference means and SDs exactly reproduce the frozen previous calibration.',
        conditioning='Same omissions applied to observed and continuous data; fixed original Golf order; native rebinning.'))


if __name__=='__main__':main()
