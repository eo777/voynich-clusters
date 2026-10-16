"""Prepare observed transition features and saved fitted orders for original scEGRET R."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json

import numpy as np

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    common = module('egret_input_common', TEAM/'whole_manuscript_topology_20261015/common.py')
    driver = module('egret_input_counts', TEAM/'b_horseshoe_20261015/drivers.py')
    data = common.load_data()
    alphabet = [r['unit'] for r in common.rows(common.DATA/'alphabet.csv')]
    basins = json.loads((TEAM/'golf_stability_20261016/basins.json').read_text())
    searches = json.loads((TEAM/'golf_influence_20261016/b_permutations.json').read_text())
    original = json.loads((TEAM/'nonlinear_sequence_20261010/plot_data.json').read_text())
    out = dict(design=dict(published_source='scegret_source/R/egret_fly.R',
        inputs='Separate word and line transition features; TPM-like counts per million per profile as primary, raw event counts as a length sensitivity control.',
        pseudotime='Saved segment-projected normalized arclength; original Golf spline also included.',
        steps={'paragraph207':[10,12,14],'pool50':[6,8,10]},
        minCellCount=3,min_detect=.1,strict_source_comparisons=True,
        source_feature_formula='mean_g(abs(log1p(mean_bin1(C_g))-log1p(mean_bin2(C_g)))/log1p(mean_all(C_g)))',
        interpretation='Source-code scEGRET continuity screening applied to text features; no biological cutoff transferred.'),cohorts=[])
    for cohort,meta in [('paragraph207',data['meta']),('pool50',data['pool_meta'])]:
        obsids=[m['id'] for m in meta]
        matrices={}
        for mode in ['word','line']:
            counts=driver.read_counts(mode,[m['profile_id'] for m in meta],alphabet).reshape(len(meta),-1)
            assert np.all(counts.sum(1)>0)
            for scale in ['per_million','raw']:
                vals=counts/counts.sum(1)[:,None]*1e6 if scale=='per_million' else counts
                filename=f'scegret_{cohort}_{mode}_{scale}.csv'
                with (HERE/filename).open('w',encoding='utf-8',newline='') as f:
                    w=csv.writer(f);w.writerow(['feature']+obsids)
                    for j,(a,b) in enumerate((a,b) for a in alphabet for b in alphabet):
                        w.writerow([a+'__to__'+b]+[format(float(v),'.12g') for v in vals[:,j]])
                matrices[mode+'_'+scale]=filename
        candidates=[]
        def add(ident,positions,description,penalty=None):
            pp=np.asarray(positions,dtype=float)
            assert pp.shape==(len(meta),) and np.all(np.isfinite(pp))
            medians={g:float(np.median(pp[[i for i,m in enumerate(meta) if m['group']==g]]))
                     for g in sorted(set(m['group'] for m in meta))}
            candidates.append(dict(id=ident,description=description,penalty=penalty,
                positions=pp.tolist(),group_medians=medians))
        for setting in basins['settings']:
            if setting['cohort']!=cohort:continue
            for r in setting['baseline_runs']:
                add(f"native_l{setting['lambda']:.2f}_{r['start']}",r['positions'],
                    'Original nine-start elastic principal-curve fit',setting['lambda'])
        for setting in searches['settings']:
            if setting['cohort']!=cohort:continue
            for h in [-1,1]:
                for t in [-1,1]:
                    eligible=[r for r in setting['runs'] if r['early_relation']==-1
                        and r['b_relations']['3a|3b']==h and r['b_relations']['4a|4c']==t]
                    best=min(eligible,key=lambda r:r['objective'])
                    add(f"B120_l{setting['penalty']:.2f}_h{h}_t{t}",best['positions'],
                        'Best recovered B120 curve for specified joint Herbal/terminal direction',setting['penalty'])
        if cohort=='paragraph207':
            add('original_golf_spline',data['global_t'],'Original Golf spline positions in shared geometry')
        else:
            lookup={r['bifolium']:r['path_position'] for r in original['sheets']}
            add('original_golf_spline',[lookup[m['bifolium']] for m in meta],
                'Original Golf spline pooled positions in shared geometry')
        out['cohorts'].append(dict(cohort=cohort,meta=meta,matrices=matrices,candidates=candidates))
    source_paths=[Path(__file__),TEAM/'golf_stability_20261016/basins.json',
                  TEAM/'golf_influence_20261016/b_permutations.json',
                  TEAM/'nonlinear_sequence_20261010/plot_data.json',
                  common.DATA/'transitions-word-counts.csv',common.DATA/'transitions-line-counts.csv']
    out['source_sha256']={str(p.relative_to(TEAM.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    (HERE/'scegret_inputs.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({c['cohort']:len(c['candidates']) for c in out['cohorts']}),flush=True)


if __name__=='__main__':main()
