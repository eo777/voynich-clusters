"""Audit the calibration accounting and render compact numerical evidence."""
from pathlib import Path
from collections import Counter
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import run_calibration as cal

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    nulls = [read(p) for p in sorted(HERE.glob('null_*.json'))]
    powers = [read(p) for p in sorted(HERE.glob('power_*.json'))]
    matched = [read(p) for p in sorted(HERE.glob('matched_*.json'))]
    bypass = [read(p) for p in sorted(HERE.glob('bypass_*.json'))]
    assert len(nulls) == 24 and len(powers) == 8 and len(matched) == 8 and len(bypass) == 16
    assert len({(r['cohort'],r['channel'],r['model'],r['bandwidth']) for r in nulls}) == 24
    assert all(r['code_sha256'] == cal.sha(HERE/'run_calibration.py') and
               r['design_sha256'] == cal.sha(HERE/'DESIGN.md') for r in nulls+powers)
    checks = 0
    for r in nulls:
        assert len(r['calibration_max_z']) == 999 and len(r['heldout_max_z']) == 500
        maximum = np.array(r['calibration_max_z'])
        threshold = np.quantile(maximum, .95, method='higher')
        assert threshold == r['threshold95_max_standardized_score']
        assert np.count_nonzero(np.array(r['heldout_max_z']) > threshold) == r['heldout_false_alarm']['events']
        for b in r['boundaries']:
            z = (b['observed']-b['null_mean'])/b['null_sd']
            assert abs(z-b['z']) < 1e-12
            p = (1+sum(float(v) >= z for v in maximum))/1000
            assert p == b['scan_tail'] and min(1,4*p) == b['whole_study_tail']
            checks += 1
    for r in powers:
        assert len(r['conditions']) == 9
        for row in r['conditions']:
            assert row['calibration_evaluable'] == 399
            assert row['gap_evaluable'] == row['random_evaluable'] == 300
    for r in matched:
        assert len(r['conditions']) == 9
        assert r['code_sha256'] == cal.sha(HERE/'matched_power.py')
        assert r['design_sha256'] == cal.sha(HERE/'MATCHED_DESIGN.md')
        ref = HERE/f"null_{r['cohort']}_{r['channel']}_{r['model']}_h0.12.json"
        assert r['source_reference_sha256'] == cal.sha(ref)
        assert all(v['detection']['trials'] == 300 for v in r['conditions'])
    for r in bypass:
        assert r['code_sha256'] == cal.sha(HERE/'bypass_calibration.py')
        assert r['generator_sha256'] == cal.sha(HERE/'run_calibration.py')
        assert r['design_sha256'] == cal.sha(HERE/'BYPASS_DESIGN.md')
        assert r['verification_max_error'] < 1e-12
        assert len(r['calibration_max_z']) == 999 and len(r['heldout_max_z']) == 500
        for b in r['boundaries']:
            p = (1+sum(v >= b['z'] for v in r['calibration_max_z']))/1000
            assert b['two_family_adjusted_tail'] == min(1,8*p)
    cases = cal.load()
    compositions = []
    for c in cases:
        if c['mode'] != 'word':
            continue
        for steps in c['steps']:
            assignment, _, occupancy = cal.audit.bins(c['t'], steps)
            compositions.append(dict(cohort=c['cohort'], steps=steps, bins=[dict(bin=b,
                n=int(occupancy[b-1]), groups=dict(Counter(c['meta'][i]['group'] or 'unassigned'
                    for i in range(len(assignment)) if assignment[i] == b)),
                members=[c['ids'][i] for i in range(len(assignment)) if assignment[i] == b])
                for b in range(1,steps+1)]))
    primary = [r for r in nulls if r['model'] == 'sheet_residual' and r['bandwidth'] == .12]
    primary_h = [dict(cohort=r['cohort'], channel=r['channel'], **b) for r in primary
        for b in r['boundaries'] if b['touches_H'] and b['scale']=='per_million'
        and b['steps'] == (12 if r['cohort']=='paragraph207' else 8)]
    flagged = [dict(cohort=r['cohort'],channel=r['channel'],bandwidth=r['bandwidth'],**b)
        for r in nulls if r['model']=='sheet_residual'
        for b in r['boundaries'] if b['whole_study_tail'] <= .05]
    all_h = [dict(cohort=r['cohort'],channel=r['channel'],model=r['model'],bandwidth=r['bandwidth'],**b)
        for r in nulls for b in r['boundaries'] if b['touches_H']]
    count_null = 24*(400+999+500)
    count_delete = 8*9*(399+300)
    count_matched = 8*9*300
    count_bypass = 16*(400+999+500)
    out = dict(status='complete', simulation_counts=dict(continuous=count_null,
        deletion_controls=count_delete, matched_gap_controls=count_matched,
        supplementary_bypass_controls=count_bypass,
        total=count_null+count_delete+count_matched+count_bypass), primary_H=primary_h,
        primary_model_flagged_boundaries=flagged, all_H_boundaries=all_h,
        bin_compositions=compositions, diagnostics=[{k:r[k] for k in
            ['cohort','channel','model','bandwidth','diagnostics','heldout_false_alarm','overall_whole_study_tail']}
            for r in nulls], literal_deletion=powers, matched_gap=matched,
        bypass_H=[dict(cohort=r['cohort'],channel=r['channel'],model=r['model'],bandwidth=r['bandwidth'],**b)
            for r in bypass for b in r['boundaries'] if b['touches_H']])
    cal.write('summary.json', out)
    audit = dict(status='passed', continuous_reference_models=24, deletion_conditions=72,
        matched_conditions=72, supplementary_bypass_models=16, observed_boundary_checks=checks,
        numerical_scorer_verification=read(HERE/'scorer_verification.json'),
        simulation_counts=out['simulation_counts'],
        result_hashes={p.name:cal.sha(p) for pattern in ['null_*.json','power_*.json','matched_*.json','bypass_*.json'] for p in sorted(HERE.glob(pattern))})
    cal.write('audit.json', audit)
    plot(nulls, matched)
    tables(out, nulls)
    feature_attribution(cases)
    print(json.dumps(dict(status='complete', simulations=out['simulation_counts'], flagged_boundaries=len(flagged))), flush=True)


def plot(nulls, matched):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.labelcolor':'#243c50','text.color':'#243c50'})
    fig, axes = plt.subplots(1,3,figsize=(13.4,4.8),gridspec_kw={'width_ratios':[1,1,1.3]})
    fig.patch.set_facecolor('white')
    rows = [('word',.08),('word',.12),('word',.20),('line',.08),('line',.12),('line',.20)]
    colors = {'word':'#137c91','line':'#ad5e19'}
    for ax,pair,title in zip(axes[:2],[(7,8),(8,9)],['Entering the f58 bin','Leaving the f58 bin']):
        for j,(channel,h) in enumerate(rows):
            r=next(r for r in nulls if r['cohort']=='paragraph207' and r['channel']==channel
                and r['bandwidth']==h and r['model']=='sheet_residual')
            b=next(b for b in r['boundaries'] if b['steps']==12 and b['scale']=='per_million'
                and (b['left'],b['right'])==pair)
            lo,hi=b['reference95']
            ax.plot([lo,hi],[j,j],color=colors[channel],lw=5,alpha=.32,solid_capstyle='round')
            ax.scatter([b['null_mean']],[j],color=colors[channel],s=26,zorder=3)
            ax.scatter([b['observed']],[j],color='#262b31',marker='D',s=40,zorder=4)
        ax.set_title(title,loc='left',fontweight='bold',pad=14)
        ax.set_yticks(range(6),[f"{ch.title()} · h={h:.2f}" for ch,h in rows])
        ax.invert_yaxis();ax.set_ylim(5.65,-.65)
        ax.set_xlim(.10,.38);ax.set_xlabel('scEGRET feature-change score')
        ax.grid(axis='x',color='#e5eaee');ax.set_axisbelow(True)
    axes[1].set_yticklabels([])
    ax=axes[2]
    for cohort,marker in [('paragraph207','o'),('pool50','s')]:
        for channel in ['word','line']:
            r=next(r for r in matched if r['cohort']==cohort and r['channel']==channel and r['model']=='sheet_residual')
            points=[v for v in r['conditions'] if v['center']==.65]
            offset=(-.0075 if cohort=='paragraph207' else .0075)+(-.002 if channel=='word' else .002)
            xs=np.array([v['width'] for v in points])+offset
            yy=np.array([v['detection']['rate']*100 for v in points])
            lower=yy-np.array([v['detection']['interval95'][0]*100 for v in points])
            upper=np.array([v['detection']['interval95'][1]*100 for v in points])-yy
            label=('Pages' if cohort=='paragraph207' else 'Bifolia')+' · '+channel
            ax.errorbar(xs,yy,yerr=[lower,upper],color=colors[channel],marker=marker,
                linestyle='-' if cohort=='paragraph207' else '--',ms=5,capsize=2,label=label)
    ax.set_title('Known gaps near the H region',loc='left',fontweight='bold',pad=14)
    ax.set_xticks([.1,.2,.3],['10%','20%','30%']);ax.set_xlim(.075,.325)
    ax.set_xlabel('Omitted span of latent feature path')
    ax.set_ylabel('Controls detected (%)');ax.set_ylim(-.5,20)
    ax.grid(axis='y',color='#e5eaee');ax.set_axisbelow(True)
    ax.legend(loc='upper left',fontsize=8,frameon=False)
    fig.suptitle('scEGRET calibration · f58 and known gaps',
        x=.055,y=.99,ha='left',fontsize=15,fontweight='bold')
    fig.text(.055,.075,'Left: continuous-control 95% range and mean; black diamond = observed score.  h = smoothing bandwidth.',fontsize=9)
    fig.text(.055,.028,'Right: fixed-count gap controls, center 0.65; exact scan-adjusted rule; 300 draws per point, 95% Wilson intervals.  Sheet-resampling model.',fontsize=9)
    fig.subplots_adjust(left=.12,right=.985,top=.86,bottom=.25,wspace=.38)
    for ext in ['png','svg']:
        fig.savefig(HERE/f'calibration.{ext}',dpi=200,facecolor='white')
    plt.close(fig)


def tables(summary, nulls):
    lines=['# Calibration numerical tables','',
        'All tails below are conditional simulation references. “Adjusted” scans boundaries, three bin counts, two scales, and applies a factor of four across channels/cohorts. They are not historical probabilities.','',
        '## H boundaries at primary bandwidth 0.12, sheet-resampling model','',
        '| Cohort | Channel | Bins | Observed | Continuous mean | Continuous 95% range | Local tail | Adjusted tail |',
        '|---|---|---|---:|---:|---|---:|---:|']
    for r in summary['primary_H']:
        lo,hi=r['reference95']
        lines.append(f"| {r['cohort']} | {r['channel']} | {r['left']} → {r['right']} | {r['observed']:.3f} | {r['null_mean']:.3f} | {lo:.3f}–{hi:.3f} | {r['local_tail']:.3f} | {r['whole_study_tail']:.3f} |")
    lines+=['','## Supplementary two-sided excursion diagnostic','',
        'B = min(E entering, E leaving) − E bypass. This is our diagnostic, not a published scEGRET test. Adjusted tails also include a factor of two for the additional statistic family.','',
        '| Cohort | Channel | Model | Bandwidth | Bins | Scale | Middle bin | Observed B | Local tail | Adjusted tail |',
        '|---|---|---|---:|---:|---|---:|---:|---:|---:|']
    for r in summary['bypass_H']:
        if r['bandwidth'] != .12 or r['scale'] != 'per_million':continue
        lines.append(f"| {r['cohort']} | {r['channel']} | {r['model']} | {r['bandwidth']:.2f} | {r['steps']} | {r['scale']} | {r['middle']} | {r['observed']:.3f} | {r['local_tail']:.3f} | {r['two_family_adjusted_tail']:.3f} |")
    lines+=['','## Other boundaries flagged under sheet resampling','',
        'All qualifying results are retained, including scale and bandwidth sensitivity. These are deviations from a specified smooth reference, not proof of missing text.','',
        '| Cohort | Channel | Bandwidth | Bins | Scale | Original bins | Score | Adjusted tail |',
        '|---|---|---:|---:|---|---|---:|---:|']
    for r in summary['primary_model_flagged_boundaries']:
        lines.append(f"| {r['cohort']} | {r['channel']} | {r['bandwidth']:.2f} | {r['steps']} | {r['scale']} | {r['left']} → {r['right']} | {r['observed']:.3f} | {r['whole_study_tail']:.3f} |")
    lines+=['','## Model-sensitive H-neighborhood flags','',
        'These are all H-touching boundaries that pass the adjusted screen under any sensitivity model. They must be read alongside results under empirical sheet resampling and per-million scaling.','',
        '| Cohort | Channel | Model | Bandwidth | Bins | Scale | Boundary | Adjusted tail |',
        '|---|---|---|---:|---:|---|---|---:|']
    for r in summary['all_H_boundaries']:
        if r['whole_study_tail']>.05:continue
        lines.append(f"| {r['cohort']} | {r['channel']} | {r['model']} | {r['bandwidth']:.2f} | {r['steps']} | {r['scale']} | {r['left']} → {r['right']} | {r['whole_study_tail']:.3f} |")
    lines+=['','## Exact decision-rule power: gap join at 0.65','',
        '| Cohort | Channel | Noise model | Gap width | Detected / 300 | Detection rate (95% interval) | Continuous false alarms / 500 |',
        '|---|---|---|---:|---:|---|---:|']
    for r in summary['matched_gap']:
        for row in r['conditions']:
            if row['center'] != .65:
                continue
            v=row['detection'];lo,hi=v['interval95']
            lines.append(f"| {r['cohort']} | {r['channel']} | {r['model']} | {row['width']:.0%} | {v['events']} | {v['rate']:.1%} ({lo:.1%}–{hi:.1%}) | {r['continuous_heldout_false_alarm']['events']} |")
    lines+=['','## Literal removal versus same-size random thinning','',
        'This uses the paper-style largest adjacent CPM score at the middle bin count, with a separate same-size random-deletion reference. It does not use the standardized multi-view decision rule above.','',
        '| Cohort | Channel | Model | Center | Span | Removed | Gap detection | Random deletion false alarms |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for r in summary['literal_deletion']:
        for row in r['conditions']:
            lines.append(f"| {r['cohort']} | {r['channel']} | {r['model']} | {row['center']:.2f} | {row['width']:.0%} | {row['removed']} | {row['gap_detection']['rate']:.1%} | {row['random_deletion_false_alarm']['rate']:.1%} |")
    lines+=['','## Null-model diagnostics','',
        'Zero fractions refer to the channel’s features observed at least once anywhere; all-zero features have no effect on EGRET. Diagnostics are essential: a conditional calibration is only as credible as its generator.','',
        '| Cohort | Channel | Model | Bandwidth | Detected features: data / median simulated | Zero fraction: data / median simulated | Residual squared error: data / median simulated | Held-out 5% scan false alarms |',
        '|---|---|---|---:|---|---|---|---:|']
    for r in summary['diagnostics']:
        d=r['diagnostics'];f=d['detected_features'];z=d['zero_fraction'];e=d['mean_squared_residual']
        lines.append(f"| {r['cohort']} | {r['channel']} | {r['model']} | {r['bandwidth']:.2f} | {f['observed']} / {f['simulated_median']:.0f} | {z['observed']:.3f} / {z['simulated_median']:.3f} | {e['observed']:.5f} / {e['simulated_median']:.5f} | {r['heldout_false_alarm']['rate']:.1%} |")
    (HERE/'TABLES.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def feature_attribution(cases):
    results=[]
    for c in cases:
        if c['cohort']!='paragraph207':continue
        x=c['x'];selected=np.count_nonzero(x,axis=0)>len(x)*.1
        values=(x/x.sum(axis=1)[:,None]*1e6)[:,selected]
        denom=np.log1p(values.mean(axis=0));names=np.array(c['features'])[selected]
        a,_,_=cal.audit.bins(c['t'],12)
        sets={'bin7':a==7,'bin8':a==8,'bin9':a==9,
            'f58rv':np.array([v in ['f58r','f58v'] for v in c['ids']]),
            'other_bin8':(a==8)&np.array([v not in ['f58r','f58v'] for v in c['ids']])}
        means={k:np.log1p(values[v].mean(axis=0)) for k,v in sets.items()}
        pairs=[('bin7','bin8'),('bin8','bin9'),('bin7','bin9'),('f58rv','bin7'),
            ('f58rv','bin9'),('other_bin8','bin7'),('other_bin8','bin9'),('f58rv','other_bin8')]
        contrasts=[]
        for left,right in pairs:
            contributions=np.abs(means[left]-means[right])/denom/len(denom)
            top=np.argsort(contributions)[-15:][::-1]
            contrasts.append(dict(left=left,right=right,score=float(contributions.sum()),
                top_features=[dict(feature=names[j],contribution=float(contributions[j])) for j in top]))
        results.append(dict(channel=c['mode'],groups={k:np.array(c['ids'])[v].tolist() for k,v in sets.items()},
            contrasts=contrasts))
    cal.write('feature_attribution.json',dict(status='complete',runs=results,
        warning='Custom subgroup contrasts retain the full-cohort feature filter and denominator. Two- and three-profile subgroups do not meet original EGRET bin eligibility; their scores are descriptive attribution only, without a separate significance test.'))


if __name__=='__main__':
    main()
