"""Within-mask, within-penalty objective competition among ensemble starts.

This complements deliberately balanced route frequencies. It does not compare
training energies across masks or penalties, nor infer chronology probabilities.
"""
from pathlib import Path
import collections
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
TOL=1e-7
ENERGY_TIE=1e-12
PRIMARY=['1a|1b','3a|3b','4a|4c']
DIRECTIONS=['3a|3b','4a|4c','3c|4a','3c|4c']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def sign_for(r,key):
    a,b=key.split('|');medians=r['group_readout_all']['group_medians'];delta=medians[a]-medians[b]
    return -1 if delta < -TOL else 1 if delta > TOL else 0


def compact(r,meta):
    labels=np.array([m['group'] for m in meta]);t=np.array(r['positions'])
    medians={g:float(np.median(t[labels==g])) for g in sorted(set(labels))}
    return dict(start_index=r['start_index'],start_id=r['start_id'],objective=r['objective'],
        all_reference_group_medians=medians,
        groups_sorted_by_median=sorted(medians,key=medians.get),
        primary_signs={key:sign_for(r,key) for key in PRIMARY},
        direction_signs={key:sign_for(r,key) for key in DIRECTIONS})


def group_signatures(winners):
    return [dict(signs=list(k),count=v) for k,v in sorted(collections.Counter(
        tuple(w['primary_signs'][key] for key in PRIMARY) for w in winners).items())]


def analyze_scope(runs,meta):
    records=[]
    for replicate in range(100):
        candidates=[r for r in runs if r['replicate']==replicate]
        record=dict(replicate=replicate,candidates=len(candidates),winner=None,objective_tied_best_starts=[],directions={})
        if not candidates:
            records.append(record);continue
        winner=min(candidates,key=lambda r:(r['objective'],r['start_index']))
        record['winner']=compact(winner,meta)
        record['objective_tied_best_starts']=[r['start_index'] for r in candidates if r['objective']<=winner['objective']+ENERGY_TIE]
        for key in DIRECTIONS:
            bysign={s:[r for r in candidates if sign_for(r,key)==s] for s in [-1,0,1]}
            before=min(bysign[-1],key=lambda r:r['objective']) if bysign[-1] else None
            after=min(bysign[1],key=lambda r:r['objective']) if bysign[1] else None
            gap=None if before is None or after is None else 100*(after['objective']-before['objective'])/min(before['objective'],after['objective'])
            difference=None if before is None or after is None else after['objective']-before['objective']
            record['directions'][key]=dict(before=compact(before,meta) if before else None,
                after=compact(after,meta) if after else None,tied_direction_candidates=len(bysign[0]),
                gap_percent_after_minus_before=gap,
                objective_preference=None if difference is None else 0 if abs(difference)<=ENERGY_TIE else 1 if difference>0 else -1)
        records.append(record)
    winners=[r['winner'] for r in records if r['winner'] is not None]
    direction_summary={}
    for key in DIRECTIONS:
        available=[r['directions'][key] for r in records if key in r['directions']]
        gaps=[r['gap_percent_after_minus_before'] for r in available if r['gap_percent_after_minus_before'] is not None]
        preferences=[r['objective_preference'] for r in available if r['objective_preference'] is not None]
        direction_summary[key]=dict(both_directions=len(gaps),
            only_before=sum(r['before'] is not None and r['after'] is None for r in available),
            only_after=sum(r['after'] is not None and r['before'] is None for r in available),
            neither=sum(r['before'] is None and r['after'] is None for r in available),
            before_lower=preferences.count(1),after_lower=preferences.count(-1),
            numerical_energy_ties=preferences.count(0),
            gap_percent_quantiles=np.quantile(gaps,[0,.05,.25,.5,.75,.95,1]).tolist() if gaps else None)
    quality=[]
    for start in range(8):
        selected=[r for r in runs if r['start_index']==start]
        gaps=[];best_count=0
        for r in selected:
            winner=records[r['replicate']]['winner'];best=winner['objective']
            gaps.append(100*(r['objective']-best)/best)
            best_count+=r['objective']<=best+ENERGY_TIE
        quality.append(dict(start_index=start,candidates=len(selected),
            tied_for_best=best_count,
            excess_percent_quantiles=np.quantile(gaps,[0,.05,.25,.5,.75,.95,1]).tolist() if gaps else None))
    return dict(mask_records=records,summary=dict(masks_with_candidates=len(winners),
        winner_primary_signatures=group_signatures(winners),directions=direction_summary,per_start_quality=quality))


def main():
    index=json.loads((HERE/'ensemble_index.json').read_text());assert index['status']=='complete'
    output=dict(status='complete',source_sha256={'ensemble_competition.py':sha(Path(__file__)),'ensemble_index.json':sha(HERE/'ensemble_index.json')},
        method='Best recovered of the declared eight starts, comparing training objectives only within identical mask and lambda.',
        gap_definition='100*(E_after-E_before)/min(E_after,E_before); positive favors first-listed group.',
        group_direction_readout='All-reference group medians of closest-segment normalized arc length.',
        directional_contrasts=DIRECTIONS,
        tie_tolerance=TOL,absolute_objective_tie_tolerance=ENERGY_TIE,
        scope_note='final_1a_first filters fitted outcomes; absent alternatives remain search-limited, not impossible.',settings=[])
    for filename in index['setting_files']:
        path=HERE/filename;obj=json.loads(path.read_text());assert obj['status']=='complete'
        output['source_sha256'][filename]=sha(path);d=obj['design'];runs=obj['runs']
        output['settings'].append(dict(cohort=d['cohort'],sampler=d['sampler'],penalty=d['penalty'],starts=d['starts'],
            scopes={'all_starts':analyze_scope(runs,d['meta']),
                'final_1a_first':analyze_scope([r for r in runs if r['group_readout_all']['contrasts']['1a|1b']['sign']==-1],d['meta'])}))
    payload=json.dumps(output,separators=(',',':'),allow_nan=False);assert len(payload.encode())<10_000_000
    (HERE/'ensemble_competition.json').write_text(payload,encoding='utf-8')
    print(json.dumps([dict(cohort=s['cohort'],sampler=s['sampler'],penalty=s['penalty'],
        final_1a_first=s['scopes']['final_1a_first']['summary']) for s in output['settings']],indent=2))


if __name__=='__main__':main()
