"""Verify production records against the two pre-production RF check fits."""
from pathlib import Path
import hashlib
import json
import numpy as np
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    before=json.loads((HERE/'execution_verification.json').read_text())
    index=json.loads((HERE/'ensemble_index.json').read_text())
    assert index['status']=='complete' and index['completed_fits']==index['converged_fits']==4800
    checks=[]
    for pilot in before['rf_checks']:
        name='ensemble_'+pilot['setting']+'.json'
        obj=json.loads((HERE/name).read_text())
        r=next(r for r in obj['runs'] if r['replicate']==0 and r['start_index']==0)
        assert r['full_node_sha256']==pilot['full_node_sha256']
        assert r['objective']==pilot['objective']
        np.testing.assert_array_equal(r['positions'],pilot['positions'])
        checks.append(dict(setting=pilot['setting'],full_dimensional_node_hash_exact=True,
            positions_exact=True,objective_exact=True,result_sha256=sha(HERE/name)))
    out=dict(status='passed',production_records_checked=len(checks),checks=checks,new_growth_calls=0,
        source_sha256={n:sha(HERE/n) for n in ['verify_completed.py','ensemble_index.json','execution_verification.json']})
    (HERE/'production_replay_verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',production_records_checked=len(checks),new_growth_calls=0)))

if __name__=='__main__':main()
