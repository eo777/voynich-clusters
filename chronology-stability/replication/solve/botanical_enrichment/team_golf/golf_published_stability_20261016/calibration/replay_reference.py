"""Deterministic replay of two complete reference models without overwriting them."""
import json
from pathlib import Path
import run_calibration as cal

HERE=Path(__file__).resolve().parent


def main():
    cases=cal.load();checks=[]
    with cal.threadpool_limits(limits=1):
        for case_i in [0,3]:
            c=cases[case_i]
            path=HERE/f"null_{c['cohort']}_{c['mode']}_sheet_residual_h0.12.json"
            saved=json.loads(path.read_text())
            replay=cal.calibrate(c,.12,'sheet_residual',[case_i,0,1])
            assert all(replay[k]==saved[k] for k in replay),path.name
            checks.append(dict(file=path.name,sha256=cal.sha(path),result='exact equality of every reference, boundary, diagnostic and held-out field'))
    cal.write('replay_verification.json',dict(status='passed',checks=checks,
        additional_validation_draws=2*(400+999+500),not_in_production_count=True,
        code_sha256=cal.sha(Path(__file__))))
    print(json.dumps(dict(status='passed',replayed_models=len(checks))),flush=True)


if __name__=='__main__':main()
