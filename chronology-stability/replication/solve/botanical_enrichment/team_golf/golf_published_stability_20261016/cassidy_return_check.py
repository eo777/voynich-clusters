"""Check native-H objective preference after the 36 declared test excursions."""
from pathlib import Path
import hashlib
import json
import numpy as np
from threadpoolctl import threadpool_limits
import cassidy_h_continuation as study

HERE=Path(__file__).resolve().parent


def main():
    source=HERE/'cassidy_h_results.json'
    result=json.loads(source.read_text())
    data=study.grid.COMMON.load_data();x,_=study.grid.COMMON.embed(data['pool_ds']['combined'])
    basins=json.loads((HERE.parent/'golf_stability_20261016/basins.json').read_text())
    setting=next(s for s in basins['settings'] if s['cohort']=='pool50' and s['lambda']==.01)
    enriched=json.loads((study.OLD/'expanded_influence_grid.json').read_text())
    extra=next(s for s in enriched['seed_sets'] if s['penalty']==.01)
    seeds={r['start']:np.asarray(r['nodes']) for r in setting['baseline_runs']}
    seeds.update({r['name']:np.asarray(r['nodes']) for r in extra['added_full_models']})
    hi=result['h_index'];meta=result['meta']
    baseline=[];returns=[]
    with threadpool_limits(limits=1):
        for name,nodes in seeds.items():
            objective,signs=study.summary(x,nodes,meta)
            baseline.append(dict(seed=name,objective=objective,signs=signs))
        for case in result['path_refinement_checks']:
            row=next(r for r in result['settings'] if r['plane']==case['plane'])
            index=case['cell_index'];plane=case['plane'];start=x[hi,plane]
            target=np.array([row['xs'][index%25],row['ys'][index//25]])
            original=seeds[case['seed']]
            outbound,a=study.path(x,original,hi,plane,[start,target])
            returned,b=study.path(x,outbound,hi,plane,[target,start])
            objective,signs=study.summary(x,returned,meta)
            baseline_objective=next(r['objective'] for r in baseline if r['seed']==case['seed'])
            returns.append(dict(seed=case['seed'],plane=plane,excursion_cell=index,
                excursion_coordinates=target.tolist(),native_objective=objective,signs=signs,
                baseline_objective=baseline_objective,objective_change=objective-baseline_objective,
                returned_node_rms=float(np.sqrt(np.mean((returned-original)**2))),
                return_node_sha256=hashlib.sha256(np.asarray(returned,dtype='<f8').tobytes()).hexdigest(),
                outbound_events=a['events'],return_events=b['events']))
    original_e=[r['objective'] for r in baseline];original_s=[r['signs'] for r in baseline]
    combined_e=original_e+[r['native_objective'] for r in returns]
    combined_s=original_s+[r['signs'] for r in returns]
    output=dict(status='complete',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='The same36 selected direct-path cases used for half-nudge checking, returned to actualH; two of ten seeds, nine grid points per plane. Not exhaustive excursion search.',
        baseline=baseline,returns=returns,
        native_H_herbal_original=study.preference(original_e,original_s),
        native_H_herbal_with_selected_returns=study.preference(combined_e,combined_s),
        native_H_terminal_original=study.preference(original_e,original_s,2),
        native_H_terminal_with_selected_returns=study.preference(combined_e,combined_s,2),
        improved_return_cases=sum(r['objective_change'] < -1e-12 for r in returns),
        worse_return_cases=sum(r['objective_change'] > 1e-12 for r in returns),
        maximum_objective_improvement=max(-r['objective_change'] for r in returns),
        counts=study.COUNTERS)
    study.save('cassidy_return_check.json',output)
    print(json.dumps({k:v for k,v in output.items() if k not in ['baseline','returns']},indent=2))


if __name__=='__main__':
    main()
