"""Retain completed original-R partitions before their authorized worker stop."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
JOBS=[dict(pid=68692,segment='pages_remaining',source='scegret_job_pages_remaining.json',
           destination='scegret_job_page_word_raw.json',cohort='paragraph207',matrices=['word_raw']),
      dict(pid=50548,segment='pools',source='scegret_job_pools.json',
           destination='scegret_job_pool_words.json',cohort='pool50',matrices=['word_per_million','word_raw'])]


def main():
    inputs=json.loads((HERE/'scegret_inputs.json').read_text())
    reports=[]
    for job in JOBS:
        raw=(HERE/job['source']).read_bytes()
        try:
            data=json.loads(raw)
        except json.JSONDecodeError:
            reports.append(dict(pid=job['pid'],status='checkpoint_write_in_progress'))
            continue
        source_cohort=next(c for c in data['cohorts'] if c['cohort']==job['cohort'])
        runs=[r for r in source_cohort['runs'] if r['matrix'] in job['matrices']]
        input_cohort=next(c for c in inputs['cohorts'] if c['cohort']==job['cohort'])
        wanted={(c['id'],m,s) for c in input_cohort['candidates'] for m in job['matrices']
                for s in inputs['design']['steps'][job['cohort']]}
        actual={(r['candidate'],r['matrix'],r['steps']) for r in runs}
        assert len(actual)==len(runs) and actual<=wanted
        if actual!=wanted:
            reports.append(dict(pid=job['pid'],status='running',completed=len(runs),expected=len(wanted)))
            continue
        assert all(r['status'] in ['complete','insufficient_occupied_bins'] for r in runs)
        destination=HERE/job['destination']
        if destination.exists():
            retained=json.loads(destination.read_text())
            assert retained['cohorts'][0]['runs']==runs
        else:
            retained=dict(data)
            retained.update(status='complete',run_count=len(runs),
                assignment_encoding='original_cut_bin_labels',
                partition=dict(cohort=job['cohort'],matrices=job['matrices'],
                    source_checkpoint=job['source'],source_checkpoint_sha256=hashlib.sha256(raw).hexdigest(),
                    source_checkpoint_saved_cases=sum(len(c['runs']) for c in data['cohorts']),
                    other_saved_cases_not_retained=sum(len(c['runs']) for c in data['cohorts'])-len(runs),
                    original_worker_pid=job['pid'],original_worker_segment=job['segment'],
                    retention_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    scope='Completed original-R partition only; other partitions are required for all648 cases.'))
            retained['cohorts']=[dict(source_cohort,runs=runs)]
            temporary=destination.with_suffix('.json.tmp')
            temporary.write_text(json.dumps(retained,indent=2)+'\n')
            temporary.replace(destination)
        reports.append(dict(pid=job['pid'],segment=job['segment'],status='partition_retained_worker_may_stop',
            cases=len(runs),destination=job['destination'],sha256=hashlib.sha256(destination.read_bytes()).hexdigest()))
    print(json.dumps(reports))


if __name__=='__main__':
    main()
