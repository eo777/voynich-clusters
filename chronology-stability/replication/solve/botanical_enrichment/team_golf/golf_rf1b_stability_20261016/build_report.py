"""Write the comparison report, verify rendered artifacts and freeze the release."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import hashlib
import json
import xml.etree.ElementTree as ET
from PIL import Image
from markdown_it import MarkdownIt

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'golf_published_stability_20261016'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pct(n,d):return f'{100*n/d:.0f}%' if d else 'unavailable'
def key(s):return (s['cohort'],s['sampler'],s['penalty'])

def main():
    c=json.loads((HERE/'comparison.json').read_text());index=json.loads((HERE/'ensemble_index.json').read_text())
    assert c['status']==index['status']=='complete'
    assert json.loads((HERE/'ensemble_audit.json').read_text())['status']=='passed'
    assert json.loads((HERE/'production_replay_verification.json').read_text())['status']=='passed'
    inputs=json.loads((HERE/'inputs.json').read_text())
    rows=[];reads=[]
    for s in c['settings']:
        name='Pages, individual sampling' if s['cohort']=='paragraph207' and s['sampler']=='point' else 'Pages, whole-sheet sampling' if s['cohort']=='paragraph207' else 'Pooled bifolia'
        z,r=[s['winners'][t]['final_1a_first'] for t in ['zl','rf']]
        def count(d,p):return f"{d['directions'][p]['before']}/{d['masks']}"
        rows.append(f"| {name} | {s['penalty']:.2f} | {count(z,'3a|3b')} | {count(r,'3a|3b')} | {count(z,'4a|4c')} | {count(r,'4a|4c')} |")
        for pair,desc in [('3a|3b','3a before 3b'),('4a|4c','Bio before Stars 4c')]:
            q=s['group_contrasts'][pair]
            reads.append(f"| {name} | {s['penalty']:.2f} | {desc} | {pct(q['zl_before'],800)} | {pct(q['rf_before'],800)} | {q['strict_paired_reversals']}/800 |")
    changes=[]
    for s in c['settings']:
        if s['cohort']=='paragraph207' and s['sampler']=='point':
            for r in s['rows']:
                if r['id'] in ['f58r','f58v','f86v4']:
                    z,q,j=r['zl_quantiles'],r['rf_quantiles'],r['joint_quantiles']
                    changes.append(f"| {r['id']} | {s['penalty']:.2f} | {z[0]:.3f}–{z[4]:.3f} | {q[0]:.3f}–{q[4]:.3f} | {j[0]:.3f}–{j[4]:.3f} |")
    main_pairs=[]
    for s in c['settings']:
        if s['sampler']=='sheet':continue
        label='pages' if s['cohort']=='paragraph207' else 'bifolia'
        z,r=[s['winners'][t]['final_1a_first'] for t in ['zl','rf']]
        main_pairs.append(f"- **{label.capitalize()}, λ={s['penalty']:.2f}:** 3a-before-3b wins {z['directions']['3a|3b']['before']}/{z['masks']} ZL samples and {r['directions']['3a|3b']['before']}/{r['masks']} RF samples. Bio-before-Stars 4c wins {z['directions']['4a|4c']['before']}/{z['masks']} and {r['directions']['4a|4c']['before']}/{r['masks']}, respectively.")
    elapsed=index['elapsed_seconds']/60
    md=f'''# RF1b-er / basic EVA: ordering stability

Completed 16 October 2026. **4,800 new ElPiGraph fits**, compared with the frozen 4,800-fit ZL ensemble. All 207 paragraph pages and 50 pooled bifolia match. [Design](DESIGN.md) · [Reproduction](REPRODUCE.md) · [Numerical comparison](comparison.json) · [Fit audit](ensemble_audit.json) · [Original ZL report](../golf_published_stability_20261016/report.html).

This is the requested RF1b-er transcription with **individual EVA characters**: for example, `ch` becomes separate `c` and `h` transitions, and `ee` contains two `e` characters. The existing ZL analysis uses **Zattera glyph groups**. Both transcription and unitization change, so differences cannot be attributed to transcription alone. Golf labels and sheet membership remain fixed. The same paragraph cohort has {sum(inputs['rf_census']['readable_words'].values()):,} readable RF words versus {sum(inputs['zl_census']['readable_words'].values()):,} ZL words; source readings and word boundaries are retained rather than forced into agreement.

## RF1b-er separately

![RF1b-er basic EVA ordering sensitivity](ordering_rf1b.png)

[PNG](ordering_rf1b.png) · [SVG](ordering_rf1b.svg) · [Chart data](ordering_rf1b_data.json)

## Both representations together

![Equal-weight combined ordering sensitivity](ordering_combined.png)

[PNG](ordering_combined.png) · [SVG](ordering_combined.svg) · [Chart data](ordering_combined_data.json)

The combined chart averages the two retained-fit frequency histograms **50:50 for each observation**. It preserves separate modes, rather than replacing alternatives with a bell curve. Equal weight is an explicit sensitivity-design choice, not a probability that either representation is correct. The shared masks mean these are paired experiments, not twice as many independent samples.

## Where the representations differ

![Side-by-side representation comparison](ordering_comparison.png)

[PNG](ordering_comparison.png) · [SVG](ordering_comparison.svg) · [Chart data](ordering_comparison_data.json)

Each horizontal strip is one page or pooled bifolium. Its darker portions are positions reached more often among fits retaining that observation. Multiple separated dark portions indicate alternative placements. The same Golf colors, row order, position bins and opacity scale are used in every panel. Rows are grouped by their fixed Golf labels, not sorted into a new chronology. Positions are normalized arc length in each fitted curve: 0 and 1 indicate the oriented ends, not dates. Coordinates from different representations need not have the same intermediate metric meaning; the combined view deliberately includes that coordinate sensitivity.

## B-section ordering

The following compares the best recovered candidate **within each sample and penalty**, restricted to fitted outcomes whose group medians put 1a before 1b. This is an outcome filter, not a hard ordering constraint. The denominator is the number of samples with an eligible candidate. It is not a posterior probability. “Before” compares group medians; individual pages may interleave. These ordering tables use projections of all reference profiles to define group medians, as in the original competition analysis. The frequency strips themselves use retained observations only.

{chr(10).join(main_pairs)}

| Sampling unit | λ | ZL: 3a before 3b | RF: 3a before 3b | ZL: Bio before Stars 4c | RF: Bio before Stars 4c |
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

The frequency charts retain all starting solutions, including less competitive local optima. The table above supplies their score context. The complete [within-sample competition results](ensemble_competition.json) also report when only one of the two directions was recovered; absence from this finite search is not evidence that an alternative is impossible.

The next table uses **all eight starts**, matching the frequency charts, and separately counts paired runs whose group-median direction reverses between representations. Matching the mask and starting permutation does not guarantee that optimization reaches corresponding local minima. Ties remain in the denominators and are recorded separately in the JSON; the complement of a “before” percentage can include ties as well as the opposite direction.

| Sampling unit | λ | Relation | ZL frequency | RF frequency | Paired strict reversals |
|---|---:|---|---:|---:|---:|
{chr(10).join(reads)}

## f58 and f86v4

These are retained-fit 5th–95th percentile position spans, using individual-page sampling and all starts. They describe sensitivity, not statistical confidence intervals. A span can bridge low-frequency gaps between modes; use the strips above to see those modes.

| Page | λ | ZL position span | RF position span | Equal-weight combined span |
|---|---:|---:|---:|---:|
{chr(10).join(changes)}

## Whole-bifolium sampling control

Pages from a physical sheet are removed together in this control. The pool fits are unchanged because pools already are sheets. Native page sampling permits remaining pages on the same sheet and is therefore a different uncertainty design. The controls use the original shared sheet masks.

![RF whole-sheet control](ordering_rf1b_sheet_control.png)

[PNG](ordering_rf1b_sheet_control.png) · [SVG](ordering_rf1b_sheet_control.svg)

![Combined whole-sheet control](ordering_combined_sheet_control.png)

[PNG](ordering_combined_sheet_control.png) · [SVG](ordering_combined_sheet_control.svg) · [Side-by-side control](ordering_comparison_sheet_control.png)

## Method and checks

The ensemble follows Albergante et al. (2020), [section 2.10 and Figure 1C](https://doi.org/10.3390/e22030296), using ElPiGraph 0.3.4. For each representation, each page/pool setting has 100 independent Bernoulli 90% inclusion masks, eight starting group orders and λ=0.01/0.05; μ=0.1 and the unbranched curve grows to 28 nodes before the original fixed-topology refinement. Pages receive equal total weight per physical sheet; pools have equal weights. The original package's weighted-wrapper issue is handled using exactly the same verified explicit-mask procedure. The full design includes page-point, page-sheet and pooled-sheet ensembles.

RF word and joined-line conditional transition matrices use individual EVA letters. Counts are pooled before normalization. The two base-2 JS matrices are averaged equally, with the same at-least-three-outgoing-observations eligibility rule and unrenormalized row weights as ZL. Classical MDS keeps all positive eigenvalues above 1e-10, giving {inputs['geometry']['page']['positive_dimensions']} RF page dimensions and {inputs['geometry']['pool']['positive_dimensions']} pool dimensions. Each representation's geometry and group-centroid anchors remain fixed during resampling. This measures fit sensitivity conditional on those inputs; it does not resample transcription uncertainty or refit MDS in each draw.

The eight group permutations are the **same ZL-selected starts** as the original ensemble, with centroids recomputed in RF geometry. No fitted ZL curve is passed into RF. This holds search coverage constant; it is not a new exhaustive search for all RF-specific routes. Golf labels are retained for comparison rather than inferred again from RF.

The same counting implementation reproduces every frozen ZL transition count and both distance matrices exactly. An independent JS implementation agrees with RF within 1.2×10⁻¹⁶. The native/explicit-mask RF pilot agrees to numerical precision, and a saved ZL ensemble fit replays exactly. Two RF production records also match their pre-production checks exactly, including their full-dimensional node hashes. All {index['completed_fits']:,} production RF fits converged, and masks, weights, labels, starting permutations and run identities match across representations. The [input checks](inputs.json), [execution verification](execution_verification.json), [production replay checks](production_replay_verification.json), [ensemble audit](ensemble_audit.json) and [paired comparison](comparison.json) retain details. Six single-BLAS-thread workers completed the production fitting in {elapsed:.1f} minutes. No private or unpublished source is used.

The [RF1b-er source](https://www.voynich.nu/data/RF1b-er.txt), source checksums, compact result records and replication scripts are retained. This is an ElPiGraph stability comparison; it does not rerun the scEGRET continuity calibration or H-coordinate influence map.
'''
    (HERE/'RESULTS.md').write_text(md,encoding='utf-8')
    body=MarkdownIt('commonmark',{'html':False}).enable('table').render(md)
    style='body{margin:0;color:#193244;background:#f4f7f9;font:16px/1.55 system-ui,sans-serif}main{max-width:1400px;margin:auto;padding:30px 4vw;background:white}h1{font-size:30px}h2{margin-top:2.4em}img{max-width:100%;height:auto}a{color:#086b9d}table{border-collapse:collapse;font-size:14px;width:100%;display:block;overflow-x:auto}th,td{border-bottom:1px solid #dce4e8;padding:8px 12px;text-align:left}th{background:#eef3f6}code{background:#eef3f6;padding:2px 4px}'
    html='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>RF1b-er ordering stability</title><style>'+style+'</style><main>'+body+'</main></html>'
    (HERE/'report.html').write_text(html,encoding='utf-8')
    (HERE/'index.html').write_text(html,encoding='utf-8')
    class Links(HTMLParser):
        def __init__(self):super().__init__();self.refs=[]
        def handle_starttag(self,tag,attrs):
            for k,v in attrs:
                if k in ['href','src']:self.refs.append(v)
    links=Links();links.feed(html);checked=0
    for ref in links.refs:
        p=urlsplit(ref)
        if p.scheme or not p.path:continue
        assert (HERE/unquote(p.path)).is_file(),ref;checked+=1
    images=[]
    for png in sorted(HERE.glob('ordering_*.png')):
        with Image.open(png) as im:im.verify()
        ET.parse(png.with_suffix('.svg'));images.append(png.name)
    verify=dict(status='passed',local_links=checked,png_svg_pairs=len(images),images=images,
        source_sha256={n:sha(HERE/n) for n in ['comparison.json','ensemble_audit.json','inputs.json','execution_verification.json','build_report.py','plot_comparison.py']})
    (HERE/'presentation_verification.json').write_text(json.dumps(verify,indent=2)+'\n',encoding='utf-8')
    files=[]
    for p in sorted(HERE.rglob('*')):
        rel=p.relative_to(HERE)
        if not p.is_file() or any(part in ['_deps','__pycache__'] for part in rel.parts) or p.name=='artifact_manifest.json' or p.suffix in ['.log','.pyc','.tmp']:continue
        files.append(dict(path=rel.as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    manifest=dict(status='complete',files=files,retained_bytes=sum(f['bytes'] for f in files),
        scientific_checks_passed=True,visual_inspection_recorded_separately=True,excluded_temporary_runtime=True)
    (HERE/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',links=checked,figure_pairs=len(images),retained_MiB=manifest['retained_bytes']/2**20)),flush=True)

if __name__=='__main__':main()
