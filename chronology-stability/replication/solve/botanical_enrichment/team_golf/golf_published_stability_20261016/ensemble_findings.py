"""Compact report tables from the completed ensemble; no new fits/inference."""
from pathlib import Path
import collections
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
PAIRS=['3a|3b','4a|4c','3c|4a','3c|4c']
FOCAL=['3a|3b','4a|4c']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def winner_counts(scope,key):
    signs=[r['winner']['direction_signs'][key] if r['winner'] else None for r in scope['mask_records']]
    return dict(before=signs.count(-1),tie=signs.count(0),after=signs.count(1),missing=signs.count(None))


def main():
    frequency=json.loads((HERE/'ensemble_summary.json').read_text())
    competition=json.loads((HERE/'ensemble_competition.json').read_text())
    output=dict(status='complete',source_sha256={name:sha(HERE/name) for name in
        ['ensemble_findings.py','ensemble_summary.json','ensemble_competition.json']},settings=[],page_vs_sheet=[])
    for s in frequency['settings']:
        match=next(c for c in competition['settings'] if (c['cohort'],c['sampler'],c['penalty'])==(s['cohort'],s['sampler'],s['penalty']))
        rawpath=HERE/s['source_file'];raw=json.loads(rawpath.read_text());output['source_sha256'][rawpath.name]=sha(rawpath)
        table=dict(cohort=s['cohort'],sampler=s['sampler'],penalty=s['penalty'],candidate_fits=800,independent_masks=100,
            retained_group_arrival_counts={key:s['scopes']['all_starts']['group_readouts']['retained']['pair_counts'][key] for key in PAIRS},
            all_reference_group_arrival_counts={key:s['scopes']['all_starts']['group_readouts']['all']['pair_counts'][key] for key in PAIRS},
            winners={},directional_gaps={},per_start=[],initial_3c_relation_coverage={})
        for key in ['3c|4a','3c|4c']:
            a,b=key.split('|');values=[-1 if start['anchor_order'].index(a)<start['anchor_order'].index(b) else 1 for start in raw['design']['starts']]
            table['initial_3c_relation_coverage'][key]=dict(before=values.count(-1),after=values.count(1))
        for scope in ['all_starts','final_1a_first']:
            value=match['scopes'][scope]
            table['winners'][scope]={key:winner_counts(value,key) for key in PAIRS}
            table['directional_gaps'][scope]={}
            for key in PAIRS:
                direction=value['summary']['directions'][key];q=direction['gap_percent_quantiles']
                table['directional_gaps'][scope][key]={k:v for k,v in direction.items() if k!='gap_percent_quantiles'}
                table['directional_gaps'][scope][key]['gap_percent_p05_median_p95']=[q[i] for i in [1,3,5]] if q else None
        for i,start in enumerate(raw['design']['starts']):
            runs=[r for r in raw['runs'] if r['start_index']==i]
            target=tuple(start['prior_final_signs'])
            counts=collections.Counter(tuple(r['group_readout_all']['contrasts'][key]['sign'] for key in FOCAL) for r in runs)
            early=collections.Counter(r['group_readout_all']['contrasts']['1a|1b']['sign'] for r in runs)
            other={}
            for key in ['3c|4a','3c|4c']:
                a,b=key.split('|');delta=[r['group_readout_all']['group_medians'][a]-r['group_readout_all']['group_medians'][b] for r in runs]
                other[key]=dict(before=sum(x < -1e-7 for x in delta),tie=sum(abs(x)<=1e-7 for x in delta),after=sum(x>1e-7 for x in delta))
            table['per_start'].append(dict(start_index=i,start_id=start['id'],anchor_order=start['anchor_order'],
                prior_candidate_id=start['prior_candidate_id'],declared_prior_B_signature=target,
                prior_signature_recovered=counts[target],ties_in_either_B_contrast=sum(v for k,v in counts.items() if 0 in k),
                final_joint_B_counts=[dict(signs=k,count=v) for k,v in sorted(counts.items())],
                other_B_direction_counts=other,
                early_counts=dict(before=early[-1],tie=early[0],after=early[1]),
                excess_objective_percent=match['scopes']['all_starts']['summary']['per_start_quality'][i],
                interpretation='Recovery of the declared B-family target; e1 is a mirrored early-prefix variant, not its own verified full-data baseline.'))
        output['settings'].append(table)
    for lam in [.01,.05]:
        ss={s['sampler']:s for s in frequency['settings'] if s['cohort']=='paragraph207' and s['penalty']==lam}
        files={name:HERE/s['frequency_file'] for name,s in ss.items()}
        data={name:json.loads(path.read_text()) for name,path in files.items()}
        output['source_sha256'].update({path.name:sha(path) for path in files.values()})
        assert [m['id'] for m in data['point']['meta']]==[m['id'] for m in data['sheet']['meta']]
        rows=[]
        for i,meta in enumerate(data['point']['meta']):
            a=data['point']['scopes']['all_starts']['rows'][i];b=data['sheet']['scopes']['all_starts']['rows'][i]
            wa=a['quantiles'][-1]-a['quantiles'][0];wb=b['quantiles'][-1]-b['quantiles'][0]
            rows.append(dict(id=meta['id'],group=meta['group'],bifolium=meta['bifolium'],
                native_position_width_p95_p05=wa,sheet_position_width_p95_p05=wb,difference=wb-wa))
        groups={}
        for group in sorted({r['group'] for r in rows}):
            selected=[r for r in rows if r['group']==group]
            groups[group]=dict(rows=len(selected),
                median_native_position_width=float(np.median([r['native_position_width_p95_p05'] for r in selected])),
                median_sheet_position_width=float(np.median([r['sheet_position_width_p95_p05'] for r in selected])),
                median_paired_width_difference=float(np.median([r['difference'] for r in selected])))
        output['page_vs_sheet'].append(dict(penalty=lam,group_widths=groups,
            largest_width_changes=sorted(rows,key=lambda r:abs(r['difference']),reverse=True)[:20],
            interpretation='Descriptive retained-position dispersion difference, not a confidence interval or a formal significance test; broad position variation need not imply an order reversal.'))
    output['interpretation']={
        'arrivals':'800 candidates per setting deliberately balance eight enriched starts; 100 independent masks are shared across starts.',
        'winners':'100 per-mask best-recovered candidates, comparing only identical masks and lambda; finite search, not global optimality or history probability.',
        'gaps':'Positive100*(E_after-E_before)/min(E_after,E_before) favors first-listed group; missing directions are search-limited.',
        'sampling':'Native point sampling and sheet adaptation differ; all masks use Bernoulli(.9), not exact90% sample sizes.',
        'lambda':'Enriched starting permutations are selected per cohort and lambda; contrasts between lambdas include differing start sets.'}
    output['interpretation']['3c_coverage']='Starts enrich3a/3b by4a/4c directions;3c ordering is only the coverage actually present/recovered in those eight starts, not all120 B permutations.'
    payload=json.dumps(output,separators=(',',':'),allow_nan=False);assert len(payload.encode())<2_000_000
    (HERE/'ensemble_findings.json').write_text(payload,encoding='utf-8')
    for s in output['settings']:
        print(json.dumps({k:s[k] for k in ['cohort','sampler','penalty','retained_group_arrival_counts','winners','directional_gaps']}))


if __name__=='__main__':main()
