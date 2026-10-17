"""Inspect retained observations and joint orders; no new fits or geometry."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
RF = HERE.parent
ZL = RF.parent/'golf_published_stability_20261016'
TOL = 1e-7


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def aggregate(runs, design):
    n = len(design['meta'])
    masks = np.asarray(design['masks'], dtype=bool)
    before = np.zeros((n, n), dtype=np.int32)
    ties = np.zeros_like(before)
    together = np.zeros_like(before)
    hist = np.zeros((n, 100), dtype=np.int32)
    samples = [[] for _ in range(n)]
    for run in runs:
        keep = masks[run['replicate']]
        pos = np.asarray(run['positions'])
        joint = keep[:, None] & keep[None, :]
        delta = pos[:, None]-pos[None, :]
        before += joint & (delta < -TOL)
        ties += joint & (np.abs(delta) <= TOL)
        together += joint
        for i in np.flatnonzero(keep):
            samples[i].append(pos[i])
            hist[i, min(99, int(pos[i]*100))] += 1
    assert np.array_equal(before+before.T+ties, together)
    assert np.array_equal(hist.sum(1), together.diagonal())
    upper = np.triu_indices(n, 1)
    return dict(fits=len(runs), masks=len({r['replicate'] for r in runs}),
        hist=hist.tolist(), retained=hist.sum(1).tolist(),
        quantiles=[np.quantile(v, [.05, .5, .95]).tolist() if v else None for v in samples],
        before=before[upper].tolist(), ties=ties[upper].tolist(), together=together[upper].tolist())


def pair(data, i, j, n):
    if i > j:
        b, t, total = pair(data, j, i, n)
        return total-b-t, t, total
    k = n*i-i*(i+1)//2+j-i-1
    return data['before'][k], data['ties'][k], data['together'][k]


def main():
    sources = {str(p.relative_to(RF)): sha(p) for p in [Path(__file__), HERE/'template.html']}
    settings = []
    for rep, root in [('RF', RF), ('ZL', ZL)]:
        index = json.loads((root/'ensemble_index.json').read_text())
        competition = json.loads((root/'ensemble_competition.json').read_text())
        for name in index['setting_files']:
            path = root/name
            assert sha(path) == index['result_sha256'][name]
            source = json.loads(path.read_text())
            d = source['design']
            rows = source['runs']
            eligible = [r for r in rows if r['group_readout_all']['contrasts']['1a|1b']['sign'] == -1]
            best = [min([r for r in eligible if r['replicate'] == mask], key=lambda r:(r['objective'], r['start_index']))
                    for mask in sorted({r['replicate'] for r in eligible})]
            scopes = {'all': aggregate(rows, d), 'best': aggregate(best, d)}
            freq_path = root/name.replace('ensemble_', 'ensemble_frequency_', 1)
            expected = json.loads(freq_path.read_text())['scopes']['all_starts']['rows']
            assert scopes['all']['hist'] == [r['histogram'] for r in expected]
            check = next(s for s in competition['settings'] if all(s[k] == d[k] for k in ['cohort', 'sampler', 'penalty']))
            eligible_scope = next(v for k,v in check['scopes'].items() if '1a' in k)
            expected_winners = [(r['replicate'],r['winner']['start_index']) for r in eligible_scope['mask_records'] if r['winner']]
            assert [(r['replicate'],r['start_index']) for r in best] == expected_winners
            settings.append(dict(rep=rep, cohort=d['cohort'], sampler=d['sampler'], penalty=d['penalty'],
                                 meta=d['meta'], scopes=scopes))
            sources[rep+'/'+name] = sha(path)
            sources[rep+'/'+freq_path.name] = sha(freq_path)
            sources[rep+'/ensemble_competition.json'] = sha(root/'ensemble_competition.json')
            print(rep, name, 'verified', flush=True)
    palette = json.loads((RF.parent/'nonlinear_sequence_20261010/plot_data.json').read_text())['palette']
    findings = []
    for rep in ['RF', 'ZL']:
        for cohort, focal, group, relation in [('paragraph207','f86v4','3b','before'),('pool50','I:67-68','2a','after')]:
            ss = [s for s in settings if s['rep'] == rep and s['cohort'] == cohort and s['sampler'] == 'point']
            meta = ss[0]['meta']; n = len(meta); i = next(i for i,m in enumerate(meta) if m['id'] == focal)
            for scope in ['all', 'best']:
                comparisons=[]
                for j,m in enumerate(meta):
                    if j == i or m['group'] != group: continue
                    counts = np.sum([pair(s['scopes'][scope],i,j,n) for s in ss],axis=0)
                    b,t,total = map(int,counts)
                    comparisons.append(dict(id=m['id'],before=b,tied=t,after=total-b-t,joint_fits=total,
                        direction_fraction=(b if relation=='before' else total-b-t)/total))
                findings.append(dict(rep=rep,cohort=cohort,focal=focal,scope=scope,relation=relation,
                    comparisons=comparisons,min_fraction=min(r['direction_fraction'] for r in comparisons),
                    max_fraction=max(r['direction_fraction'] for r in comparisons)))
    data = dict(settings=settings,palette=palette,tolerance=TOL,
        method='Jointly retained comparisons; ties preserved; equal weighting of selected representation/penalty settings. Best scope uses existing within-mask objective winners with 1a-before-1b.',
        source_sha256=sources,findings=findings)
    payload=json.dumps(data,separators=(',',':'),allow_nan=False)
    (HERE/'data.json').write_text(payload,encoding='utf-8')
    template=(HERE/'template.html').read_text(encoding='utf-8')
    assert template.count('__DATA__')==1
    (HERE/'index.html').write_text(template.replace('__DATA__',payload.replace('</','<\\/')),encoding='utf-8')
    (HERE/'verification.json').write_text(json.dumps(dict(status='passed',new_fits=0,settings=len(settings),
        original_histograms_exact=True,original_winners_exact=True,pair_partition_exact=True,
        findings=findings,outputs={n:sha(HERE/n) for n in ['data.json','index.html']}),indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','settings':len(settings),'html_MiB':(HERE/'index.html').stat().st_size/2**20}),flush=True)


if __name__ == '__main__': main()
