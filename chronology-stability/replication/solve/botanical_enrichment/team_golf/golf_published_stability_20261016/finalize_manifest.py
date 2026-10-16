"""Write a compact release manifest after numerical and presentation verification."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
TEMP_DIRS={'_deps','_r_runtime','__pycache__'}


def main():
    checks={}
    for name in ['ensemble_audit.json','ensemble_independent_audit.json','cassidy_independent_audit.json',
                 'scegret_independent_audit.json','presentation_verification.json',
                 'calibration/audit.json','calibration/presentation_verification.json',
                 'calibration/replay_verification.json','calibration/omission_test/audit.json',
                 'calibration/omission_test/presentation_verification.json']:
        obj=json.loads((HERE/name).read_text(encoding='utf-8'))
        assert obj.get('status') in ['passed','complete'], (name,obj.get('status'))
        checks[name]=dict(status=obj['status'],sha256=hashlib.sha256((HERE/name).read_bytes()).hexdigest())
    assert json.loads((HERE/'ensemble_index.json').read_text())['status']=='complete'
    assert json.loads((HERE/'scegret_results.json').read_text())['run_count']==648
    files=[]
    for path in sorted(HERE.rglob('*')):
        rel=path.relative_to(HERE)
        if not path.is_file() or any(p in TEMP_DIRS for p in rel.parts) or path.name=='artifact_manifest.json':
            continue
        if path.suffix in ['.log','.pyc']:
            continue
        files.append(dict(path=rel.as_posix(),bytes=path.stat().st_size,
                          sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    out=dict(status='complete',study='Published stability methods for Golf chronology',
             verification=checks,files=files,retained_bytes=sum(f['bytes'] for f in files),
             privacy='No private or unpublished study is an input.',
             temporary_runtime_note='Downloaded runtime directories are temporary; their dependency versions and research source hashes are retained separately.')
    (HERE/'artifact_manifest.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='complete',files=len(files),retained_MiB=round(out['retained_bytes']/2**20,2))))


if __name__=='__main__':main()
