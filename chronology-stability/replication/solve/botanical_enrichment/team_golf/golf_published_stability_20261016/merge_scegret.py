"""Merge disjoint original-R batches; never substitute or average duplicate cases."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
PARTS=['first_page','page_word_raw','page_line_per_million','page_line_raw','pool_words','pool_lines']


def main():
    inputs=json.loads((HERE/'scegret_inputs.json').read_text())
    merged={}; provenance={}; runtimes=[]
    for part in PARTS:
        path=HERE/f'scegret_job_{part}.json'
        obj=json.loads(path.read_text())
        assert obj['status']=='complete', part
        provenance[path.name]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                   retained_cases=sum(len(c['runs']) for c in obj['cohorts']))
        runtimes.append(obj['runtime'])
        for c in obj['cohorts']:
            for r in c['runs']:
                key=(c['cohort'],r['matrix'],r['candidate'],r['steps'])
                assert key not in merged, ('overlapping partition',key)
                assert r['status'] in ['complete','insufficient_occupied_bins'], key
                merged[key]=r
    fullruntime=[r for r in runtimes if 'packages' in r]
    assert fullruntime, 'At least one completed process must record package versions.'
    assert all(r['R']==fullruntime[0]['R'] and r['webr']==fullruntime[0]['webr'] for r in runtimes)
    packages=lambda r: sorted((p['Package'],p['Version']) for p in r['packages'])
    assert all(packages(r)==packages(fullruntime[0]) for r in fullruntime)
    out=dict(status='complete',run_count=len(merged),
        method='Unmodified egret.flap and egret.fly from pinned scEGRET source, executed in R/WebR',
        runtime=fullruntime[0],design=inputs['design'],cohorts=[],partition_provenance=provenance,
        execution_note='Disjoint feature batches ran in parallel. Completed batches from two serial workers were preserved before stopping redundant continuation; partially evaluated overlap is not retained. One separate instrumented original-function case diagnosed runtime. Main scientific results contain648 unique evaluations, not duplicate runs.')
    for c in inputs['cohorts']:
        runs=[]
        for matrix in c['matrices']:
            for candidate in c['candidates']:
                for steps in inputs['design']['steps'][c['cohort']]:
                    runs.append(merged.pop((c['cohort'],matrix,candidate['id'],steps)))
        out['cohorts'].append(dict(cohort=c['cohort'],meta=c['meta'],runs=runs))
    assert not merged and out['run_count']==648
    (HERE/'scegret_results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='complete',runs=648,partitions=len(PARTS))))


if __name__=='__main__':main()
