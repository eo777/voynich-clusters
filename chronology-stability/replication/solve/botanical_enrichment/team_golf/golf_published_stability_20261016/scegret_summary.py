"""Summarize original-R scEGRET readouts without choosing a favorable binning."""
from pathlib import Path
from itertools import combinations, product
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
B=['3a','3b','3c','4a','4c']
TOL=1e-12


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def as_list(value):return value if isinstance(value,list) else [value]


def stats(values):
    return dict(mean=float(np.mean(values)),maximum=float(np.max(values)),
        median=float(np.median(values)),sum=float(np.sum(values)),count=len(values)) if values else None


def original_assignment(positions,steps):
    """Auxiliary labels matching R cut; used only when R omitted sparse output."""
    p=np.asarray(positions,dtype=float);lo,hi=float(p.min()),float(p.max())
    width=hi-lo
    if width==0:
        change=abs(lo) if lo else 1.
        edges=np.linspace(lo-change/1000,hi+change/1000,steps+1)
    else:
        edges=np.linspace(lo,hi,steps+1);edges[0]-=width/1000;edges[-1]+=width/1000
    return np.searchsorted(edges,p,side='left').tolist()


def describe(run,meta,candidate):
    occupancy=run['occupancy'];retained=run.get('retained_bins',[])
    if run['status']!='complete':retained=[i+1 for i,n in enumerate(occupancy) if n>3]
    assignment=run.get('assignment',original_assignment(candidate['positions'],run['steps']))
    assert len(assignment)==len(meta)
    groups=sorted(set(m['group'] for m in meta))
    coverage={}
    for group in groups:
        indices=[i for i,m in enumerate(meta) if m['group']==group]
        included=[i for i in indices if assignment[i] in retained]
        excluded=[i for i in indices if assignment[i] not in retained]
        coverage[group]=dict(total=len(indices),retained=len(included),excluded=len(excluded),
            retained_bins=sorted(set(assignment[i] for i in included)),
            excluded_bins=sorted(set(assignment[i] for i in excluded)))
        if group in ['3a','3b','4a','4c']:
            coverage[group]['excluded_ids']=[meta[i]['id'] for i in excluded]
    h_rows=[dict(id=m['id'],group=m['group'],bin=assignment[i],retained=assignment[i] in retained)
            for i,m in enumerate(meta) if m['bifolium']=='H:58-65']
    focal_bins=sorted(set(r['bin'] for r in h_rows if r['id'] in ['f58r','f58v','H:58-65']))
    pair_resolution={}
    for a,b in [('3a','3b'),('4a','4c')]:
        ai=[i for i,m in enumerate(meta) if m['group']==a and assignment[i] in retained]
        bi=[i for i,m in enumerate(meta) if m['group']==b and assignment[i] in retained]
        within=sum(assignment[i]==assignment[j] for i in ai for j in bi)
        before=sum(assignment[i]<assignment[j] for i in ai for j in bi)
        after=sum(assignment[i]>assignment[j] for i in ai for j in bi)
        denominator=len(ai)*len(bi);between=before+after
        original_pairs=sum(m['group']==a for m in meta)*sum(m['group']==b for m in meta)
        pair_resolution[a+'|'+b]=dict(retained_cross_group_pairs=len(ai)*len(bi),
            original_cross_group_pairs=original_pairs,
            within_same_bin=within,first_bin_before=before,first_bin_after=after,
            both_groups_retained=bool(ai and bi),between_bin_pairs=between,
            fraction_retained_pairs_between_bins=between/denominator if denominator else None,
            fraction_original_pairs_between_bins=between/original_pairs if original_pairs else None,
            bin_separates_any_retained_pairs=between>0,
            resolution_status='no_both_group_coverage' if not denominator else
                'all_retained_pairs_same_bin' if not between else
                'partial_between_bin_resolution' if within else 'all_retained_pairs_between_bins',
            mixed_retained_bins=sorted(set(assignment[i] for i in ai)&set(assignment[j] for j in bi)))
    tied_pairs=[[a,b] for a,b in combinations(B,2)
                if abs(candidate['group_medians'][a]-candidate['group_medians'][b])<=1e-8]
    result=dict(candidate=run['candidate'],matrix=run['matrix'],steps=run['steps'],status=run['status'],
        occupancy=occupancy,retained_bins=retained,
        excluded_bins=[i+1 for i,n in enumerate(occupancy) if n<=3],
        retained_observations=sum(occupancy[i-1] for i in retained),
        total_observations=sum(occupancy),
        golf_group_coverage=coverage,H_observations=h_rows,H_focal_bins=focal_bins,
        focal_pair_bin_resolution=pair_resolution,
        B_group_median_order=sorted(B,key=lambda g:candidate['group_medians'][g]),
        B_tied_group_medians=tied_pairs,B_order_is_strict=not tied_pairs,
        B_group_medians={g:candidate['group_medians'][g] for g in B})
    if run['status']!='complete':return result
    composition=[]
    for bin_id in range(1,run['steps']+1):
        counts={g:sum(a==bin_id and m['group']==g for a,m in zip(assignment,meta))
                for g in sorted(set(m['group'] for m in meta))}
        assert sum(counts.values())==occupancy[bin_id-1]
        composition.append(dict(bin=bin_id,count=occupancy[bin_id-1],retained=bin_id in retained,
            golf_counts=counts,observations=[m['id'] for a,m in zip(assignment,meta) if a==bin_id]))
    rows=as_list(run['adjacent']);pairs=run['adjacent_bin_pairs']
    assert len(rows)==len(pairs)
    adjacent=[]
    for p,r,features in zip(pairs,rows,run['top_adjacent_features']):
        adjacent.append(dict(bins=p,value=float(r['Relative_Step_Change']),spans_excluded_bins=p[1]-p[0]>1,
            top_features=list(zip(as_list(features['features']),as_list(features['contribution'])))))
    all_values=[r['value'] for r in adjacent]
    strict_values=[r['value'] for r in adjacent if not r['spans_excluded_bins']]
    largest=max(adjacent,key=lambda r:r['value'])
    all_pair_values={tuple(int(v) for v in p):float(v) for p,v in zip(run['full_pairs'],as_list(run['full_enthalpy']))}
    bypass=[]
    for middle in focal_bins:
        left,right=middle-1,middle+1
        item=dict(H_bin=middle,bins=[left,middle,right],
            memberships=[composition[i-1] for i in [left,middle,right] if 1<=i<=run['steps']])
        if all(i in retained for i in [left,middle,right]):
            incoming=all_pair_values[(left,middle)];outgoing=all_pair_values[(middle,right)]
            direct=all_pair_values[(left,right)]
            item.update(status='available',incoming_change=incoming,outgoing_change=outgoing,bypass_change=direct,
                both_H_neighbor_changes_exceed_bypass=incoming>direct+TOL and outgoing>direct+TOL,
                top_incoming_features=next(r['top_features'] for r in adjacent if r['bins']==[left,middle]),
                top_outgoing_features=next(r['top_features'] for r in adjacent if r['bins']==[middle,right]))
        else:item.update(status='unavailable_nonretained_or_edge_bin')
        bypass.append(item)
    result.update(feature_count=run['feature_count'],adjacent=adjacent,
        source_retained_adjacent_summary=stats(all_values),strict_original_neighbor_summary=stats(strict_values),
        skipped_bin_spans=sum(r['spans_excluded_bins'] for r in adjacent),
        strict_neighbor_coverage=len(strict_values)/(run['steps']-1),
        maximum_step=largest,bin_composition=composition,
        H_bypass_comparisons=bypass,
        full_pairs=[[int(a),int(b)] for a,b in run['full_pairs']],full_pair_changes=as_list(run['full_enthalpy']))
    return result


def comparison(cases,metric):
    keys=[f'h{h}_t{t}' for h,t in product([-1,1],repeat=2)]
    wins={k:0 for k in keys};ties=0;complete_cases=[];all_same_bins=0
    for case in cases:
        values={k:(r.get(metric) or {}).get('mean') for k,r in case['directions'].items()}
        if any(v is None for v in values.values()):continue
        best=min(values.values());winning=[k for k,v in values.items() if abs(v-best)<=TOL]
        for k in winning:wins[k]+=1
        ties+=len(winning)>1
        same=len({tuple(r['retained_bins']) for r in case['directions'].values()})==1
        all_same_bins+=same
        complete_cases.append(dict(matrix=case['matrix'],steps=case['steps'],means=values,
            lowest=winning,all_four_same_retained_bin_indices=same,
            focal_retained_counts={k:{g:r['golf_group_coverage'][g]['retained'] for g in ['3a','3b','4a','4c']}
                for k,r in case['directions'].items()},
            focal_excluded_ids={k:{g:r['golf_group_coverage'][g].get('excluded_ids',[]) for g in ['3a','3b','4a','4c']}
                for k,r in case['directions'].items()},
            H_retained_ids={k:[v['id'] for v in r['H_observations'] if v['retained']] for k,r in case['directions'].items()},
            focal_pair_bin_resolution={k:r['focal_pair_bin_resolution'] for k,r in case['directions'].items()},
            strict_neighbor_coverage={k:r.get('strict_neighbor_coverage') for k,r in case['directions'].items()}))
    pairwise=[]
    for a,b in combinations(keys,2):
        differences=[c['means'][a]-c['means'][b] for c in complete_cases]
        matched=[c['means'][a]-c['means'][b] for c in complete_cases
                 if c['focal_retained_counts'][a]==c['focal_retained_counts'][b]
                 and c['H_retained_ids'][a]==c['H_retained_ids'][b]]
        matched_ids=[c['means'][a]-c['means'][b] for c in complete_cases
                     if c['focal_excluded_ids'][a]==c['focal_excluded_ids'][b]
                     and c['H_retained_ids'][a]==c['H_retained_ids'][b]]
        resolution_checks={}
        for focal in ['3a|3b','4a|4c']:
            eligible=[c for c in complete_cases if all(
                c['focal_pair_bin_resolution'][k][focal]['both_groups_retained']
                and c['focal_pair_bin_resolution'][k][focal]['bin_separates_any_retained_pairs'] for k in [a,b])]
            values=[c['means'][a]-c['means'][b] for c in eligible]
            resolution_checks[focal]=dict(cases=len(values),first_lower=sum(d < -TOL for d in values),
                tied=sum(abs(d)<=TOL for d in values),second_lower=sum(d>TOL for d in values),
                interpretation='Both candidate readouts retain both groups and separate at least one pair into different bins; this is only minimum bin resolution, not isolated evidence for that order.')
        pairwise.append(dict(first=a,second=b,first_lower=sum(d < -TOL for d in differences),
            tied=sum(abs(d)<=TOL for d in differences),second_lower=sum(d>TOL for d in differences),
            matched_focal_group_counts_and_H_ids=dict(cases=len(matched),first_lower=sum(d < -TOL for d in matched),
                tied=sum(abs(d)<=TOL for d in matched),second_lower=sum(d>TOL for d in matched)),
            matched_focal_observation_ids_and_H_ids=dict(cases=len(matched_ids),first_lower=sum(d < -TOL for d in matched_ids),
                tied=sum(abs(d)<=TOL for d in matched_ids),second_lower=sum(d>TOL for d in matched_ids)),
            both_routes_some_bin_resolution=resolution_checks,
            differences=differences,all_six_views_agree_first_lower=len(differences)==6 and all(d < -TOL for d in differences),
            all_six_views_agree_second_lower=len(differences)==6 and all(d>TOL for d in differences)))
    return dict(metric=metric,complete_cases=len(complete_cases),lowest_counts=wins,tied_lowest_cases=ties,
        all_four_same_retained_bin_indices_cases=all_same_bins,cases=complete_cases,pairwise=pairwise)


def main():
    result_path=HERE/'scegret_results.json';inputs_path=HERE/'scegret_inputs.json'
    data=json.loads(result_path.read_text());inputs=json.loads(inputs_path.read_text())
    assert data['status']=='complete' and data['run_count']==648
    out=dict(status='complete',source_sha256={p.name:sha(p) for p in [Path(__file__),result_path,inputs_path]},
        method='Original scEGRET R outputs, with declared descriptive aggregation only.',
        interpretation=['Larger relative step change means greater transition-feature contrast between bins, not physical energy.',
            'Mean and maximum step changes are descriptive summaries added here; not a published chronology-selection statistic.',
            'Rank frequencies count six dependent diagnostic views (word/line times three resolutions), not probabilities or p-values.',
            'Original scEGRET uses unweighted arithmetic bin means. Page-unit and pooled-sheet means therefore have different representation weights.',
            'Bins with at most3observations are excluded by the original R strict >3 rule. IDs and gap spans remain explicit.',
            'Matching retained-bin indices does not imply matching observation membership across different trajectories.',
            'Candidate scores must be read with group coverage: excluding a sparse Stars4c bin can reduce the reported discontinuity while discarding that group.',
            'Equal coverage is not sufficient: groups can both be absent or all their cross-group pairs can lie within one bin. Per-pair resolution flags and between-bin fractions are separate.',
            'A mixed bin cannot resolve the order of its members. H bypass contrasts compare actual source all-pair scores, without removing or refitting observations.',
            'Raw event counts are a length-sensitivity control; per-million counts are primary. Scales are not averaged together.'],
        settings=[],B_comparisons=[],original_golf=[],H_bypass_robustness=[],coverage_warnings=[])
    for cohort in data['cohorts']:
        original=next(c for c in inputs['cohorts'] if c['cohort']==cohort['cohort'])
        candidates={c['id']:c for c in original['candidates']}
        detailed=[describe(r,cohort['meta'],candidates[r['candidate']]) for r in cohort['runs']]
        index={(r['candidate'],r['matrix'],r['steps']):r for r in detailed}
        # Retain detailed bins/pairs only for the original-Golf plots; compact other runs.
        for r in detailed:
            out['settings'].append(dict(cohort=cohort['cohort'],**{k:v for k,v in r.items()
                if k not in ['bin_composition','full_pairs','full_pair_changes','adjacent']}))
            if r['candidate']=='original_golf_spline':out['original_golf'].append(dict(cohort=cohort['cohort'],**r))
        for penalty in [.01,.05]:
            for scale in ['per_million','raw']:
                cases=[]
                for channel in ['word','line']:
                    matrix=channel+'_'+scale
                    for steps in inputs['design']['steps'][cohort['cohort']]:
                        directions={f'h{h}_t{t}':index[(f'B120_l{penalty:.2f}_h{h}_t{t}',matrix,steps)]
                                    for h,t in product([-1,1],repeat=2)}
                        cases.append(dict(matrix=matrix,steps=steps,directions=directions))
                out['B_comparisons'].append(dict(cohort=cohort['cohort'],penalty=penalty,scale=scale,
                    source_adjacent=comparison(cases,'source_retained_adjacent_summary'),
                    strict_neighbor=comparison(cases,'strict_original_neighbor_summary'),
                    orders={k:dict(groups=v['B_group_median_order'],tied_pairs=v['B_tied_group_medians'],
                        strict=v['B_order_is_strict']) for k,v in cases[0]['directions'].items()}))
        if any(r.get('skipped_bin_spans',0)>0 for r in detailed):
            out['coverage_warnings'].append(dict(cohort=cohort['cohort'],
                runs_spanning_excluded_bins=sum(r.get('skipped_bin_spans',0)>0 for r in detailed),
                runs_with_excluded_bins=sum(bool(r['excluded_bins']) for r in detailed),
                instruction='Source adjacent means use neighboring retained bins, which may span missing original bins; also inspect strict-neighbor summaries and coverage.'))
        for scale in ['per_million','raw']:
            cases=[]
            for r in detailed:
                if r['candidate']!='original_golf_spline' or not r['matrix'].endswith('_'+scale):continue
                for h in r.get('H_bypass_comparisons',[]):
                    cases.append(dict(matrix=r['matrix'],steps=r['steps'],**h))
            available=[r for r in cases if r['status']=='available']
            out['H_bypass_robustness'].append(dict(cohort=cohort['cohort'],scale=scale,
                evaluated_cases=len(cases),available_three_consecutive_retained_bins=len(available),
                both_sides_exceed_bypass_cases=sum(r['both_H_neighbor_changes_exceed_bypass'] for r in available),
                cases=cases,interpretation='These compare bin-average profiles, including all other observations sharing H bins; not an H-only effect.'))
    (HERE/'scegret_summary.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n')
    print(json.dumps(dict(status=out['status'],settings=len(out['settings']),B_comparisons=len(out['B_comparisons']),
        primary=[dict(cohort=c['cohort'],penalty=c['penalty'],lowest_counts=c['source_adjacent']['lowest_counts'],
            same_bin_indices=c['source_adjacent']['all_four_same_retained_bin_indices_cases']) for c in out['B_comparisons'] if c['scale']=='per_million'],
        coverage_warnings=out['coverage_warnings']),indent=2))


if __name__=='__main__':main()
