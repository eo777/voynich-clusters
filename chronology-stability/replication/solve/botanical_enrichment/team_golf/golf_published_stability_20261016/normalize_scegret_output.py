"""Label-only correction for the initial R wrapper; scientific values are unchanged.

The executing wrapper exported factor integer codes. Its corrected source now
exports numeric level labels. This script preserves the old codes and supplies
the original bin labels, including occupied bins excluded for small size.
"""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    path = HERE / 'scegret_results.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    assert data['status'] == 'complete', 'Wait until the R process finishes.'
    if data.get('assignment_encoding') == 'original_cut_bin_labels':
        print('Already normalized; no change.')
        return
    before = digest(path)
    corrected = 0
    for cohort in data['cohorts']:
        for run in cohort['runs']:
            if run['status'] != 'complete':
                continue
            if 'factor_codes' in run:
                continue  # New worker wrappers already export original labels.
            occupied = [i+1 for i, n in enumerate(run['occupancy']) if n > 0]
            codes = run['assignment']
            assert all(1 <= k <= len(occupied) for k in codes)
            labels = [occupied[k-1] for k in codes]
            corrected += labels != codes
            run['factor_codes'] = codes
            run['assignment'] = labels
    data['assignment_encoding'] = 'original_cut_bin_labels'
    data['output_label_normalization'] = dict(
        before_sha256=before, changed_runs=corrected,
        reason='The executing wrapper used as.integer(TimeBin), which exports compressed factor codes after droplevels. Map through every occupied original bin; preserve codes separately. No score or membership changed.',
        executing_export='assignment=as.integer(obj@Metadata$TimeBin)',
        corrected_export='assignment=as.integer(as.character(obj@Metadata$TimeBin))',
        script_sha256=digest(Path(__file__)))
    manifest = json.loads((HERE/'scegret_source/source_manifest.json').read_text())
    data['published_source_manifest'] = manifest
    data['runtime']['verified_webr_package_version'] = json.loads(
        (HERE/'_r_runtime/node_modules/webr/package.json').read_text())['version']
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status='complete', changed_run_labels=corrected,
                         output_sha256=digest(path))))


if __name__ == '__main__':
    main()
