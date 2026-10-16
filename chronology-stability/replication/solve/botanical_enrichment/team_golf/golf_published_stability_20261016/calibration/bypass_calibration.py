"""Calibrate the supplementary two-sided diagnostic; it is not original EGRET."""
import json
import time
from pathlib import Path
import numpy as np
import run_calibration as cal

HERE = Path(__file__).resolve().parent


class BypassScorer(cal.Scorer):
    def __init__(self, positions, steps, h_indices=(), scales=('per_million','raw')):
        super().__init__(positions, steps, h_indices, scales)
        self.labels = []
        for step,scale,kept,weights,assignment,occupancy in self.views:
            hb=set(assignment[list(h_indices)]) if len(h_indices) else set()
            for a,b,c in zip(kept[:-2],kept[1:-1],kept[2:]):
                self.labels.append(dict(steps=step,scale=scale,left=int(a),middle=int(b),right=int(c),
                    left_n=int(occupancy[a-1]),middle_n=int(occupancy[b-1]),right_n=int(occupancy[c-1]),
                    touches_H=bool(b in hb)))

    def score(self, counts):
        detected=np.count_nonzero(counts,axis=0)>len(counts)*.1
        raw=counts[:,detected]
        cpm=raw/counts.sum(axis=1)[:,None]*1e6
        out=[]
        for _,scale,kept,weights,_,_ in self.views:
            if len(kept)<3:continue
            matrix=cpm if scale=='per_million' else raw
            means=np.log1p(weights.T@matrix)
            denominator=np.log1p(matrix.mean(axis=0))
            adjacent=np.mean(np.abs(np.diff(means,axis=0))/denominator,axis=1)
            bypass=np.mean(np.abs(means[2:]-means[:-2])/denominator,axis=1)
            out.extend(np.minimum(adjacent[:-1],adjacent[1:])-bypass)
        return np.array(out)


def main():
    start=time.monotonic()
    cases=cal.load()
    errors=[]
    for c in cases:
        reference=[]
        for steps in c['steps']:
            for scale in ['per_million','raw']:
                x=c['x'];matrix=(x/x.sum(axis=1)[:,None]*1e6 if scale=='per_million' else x).T
                e=cal.audit.expected(matrix,c['features'],c['ids'],c['t'],steps)
                lookup={tuple(pair):score for pair,score in zip(e['full_pairs'],e['full_enthalpy'])}
                kept=e['retained_bins']
                for a,b,d in zip(kept[:-2],kept[1:-1],kept[2:]):
                    reference.append(min(lookup[a,b],lookup[b,d])-lookup[a,d])
        error=float(np.max(np.abs(BypassScorer(c['t'],c['steps']).score(c['x'])-reference)))
        assert error<1e-12
        errors.append(error)
    cal.Scorer=BypassScorer
    with cal.threadpool_limits(limits=1):
        for case_i,c in enumerate(cases):
            for model_i,model in enumerate(cal.CONFIG['models']):
                for bw_i,bw in enumerate(cal.CONFIG['bandwidths']):
                    if model_i==1 and bw!=.12:continue
                    result=cal.calibrate(c,bw,model,[case_i,model_i,bw_i,129])
                    for b in result['boundaries']:
                        b['two_family_adjusted_tail']=min(1.,8*b['scan_tail'])
                    result['two_family_adjusted_tail']=min(1.,8*result['overall_scan_tail'])
                    result['statistic']='min(E_left_middle,E_middle_right)-E_left_right; supplementary, not published scEGRET criterion'
                    result['code_sha256']=cal.sha(Path(__file__))
                    result['generator_sha256']=cal.sha(HERE/'run_calibration.py')
                    result['design_sha256']=cal.sha(HERE/'BYPASS_DESIGN.md')
                    result['verification_max_error']=max(errors)
                    name=f"bypass_{c['cohort']}_{c['mode']}_{model}_h{bw:.2f}.json"
                    cal.write(name,result)
                    print(json.dumps(dict(file=name,seconds=round(time.monotonic()-start))),flush=True)


if __name__=='__main__':main()
