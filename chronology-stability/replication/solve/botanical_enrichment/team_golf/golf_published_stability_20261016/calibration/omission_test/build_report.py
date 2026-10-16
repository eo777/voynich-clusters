"""Summarize the requested targeted removals without changing earlier results."""
from pathlib import Path
import sys
import json
import re
from markdown_it import MarkdownIt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
import build_gallery as gallery


def main():
    results={}
    for mode in ['word','line']:
        for model in ['sheet_residual','dirichlet_multinomial']:
            d=json.loads((HERE/f'results_{mode}_{model}.json').read_text())
            assert d['status']=='complete' and d['full_data_reference_reproduced_exactly']
            results[mode,model]=d
    audit=json.loads((HERE/'audit.json').read_text());assert audit['status']=='passed'
    def rows(mode,model,ident):
        case=next(r for r in results[mode,model]['cases'] if r['id']==ident)
        return [b for b in case['boundaries'] if b['focal_join'] and b['steps']==12 and b['scale']=='per_million']
    primary=results['word','sheet_residual']['cases']
    assert all(b['targeted_adjusted_tail']>.05 for mode in ['word','line']
        for c in results[mode,'sheet_residual']['cases'] for b in c['boundaries'] if b['focal_join'] and b['scale']=='per_million')
    table=['| Removal | Profiles removed | Word score at H join | Line score at H join |',
           '|---|---:|---:|---:|']
    tails=['| Case | Word: sheet model | Line: sheet model | Word: independent model | Line: independent model |',
           '|---|---:|---:|---:|---:|']
    rawtails=['| Case | Word: targeted tail | Word: whole-path tail | Line: targeted tail | Line: whole-path tail |',
              '|---|---:|---:|---:|---:|']
    for case in primary:
        ident=case['id']
        table.append(f"| {case['title']} | {len(case['removed'])} | "+
            ' / '.join(f"{r['observed']:.3f}" for r in rows('word','sheet_residual',ident))+' | '+
            ' / '.join(f"{r['observed']:.3f}" for r in rows('line','sheet_residual',ident))+' |')
        tails.append('| '+case['title']+' | '+' | '.join(f"{min(b['targeted_adjusted_tail'] for b in rows(mode,model,ident)):.3f}"
            for mode,model in [('word','sheet_residual'),('line','sheet_residual'),('word','dirichlet_multinomial'),('line','dirichlet_multinomial')])+' |')
        rv=[]
        for mode in ['word','line']:
            cc=next(v for v in results[mode,'sheet_residual']['cases'] if v['id']==ident)
            bounds=[b for b in cc['boundaries'] if b['focal_join'] and b['scale']=='raw']
            rv.extend([min(b['targeted_adjusted_tail'] for b in bounds),min(b['whole_path_adjusted_tail'] for b in bounds)])
        rawtails.append('| '+case['title']+' | '+' | '.join(f'{v:.3f}' for v in rv)+' |')
    md=f'''# Removing f58 and the astro/stars material

**Completed 16 October 2026.** Seven cases including the full-data baseline; word and line features; 10/12/14 bins; normalized and raw counts; two previously calibrated noise models. No new curve fits. [Design](DESIGN.md) · [Code](run_omissions.py) · [Audit](audit.json).

## Answer

**Some versions of the test do flag a break after removing f58, but the finding depends on count scaling and the noise model.** With per-million frequencies and empirical sheet resampling, no tested removal produces a flag. Raw counts do produce H-region flags under sheet resampling, and the independent-profile Dirichlet–multinomial model also flags some normalized comparisons after the narrower removals. This is evidence of a sensitive comparison, not a robust demonstration that f58 is the surviving bridge preventing a discontinuity.

The **per-million feature contrast** actually decreases after removing f58. At 12 bins, the new comparison is approximately **0.168 word / 0.177 line**, versus **0.271 word / 0.279–0.286 line** entering and leaving the original focal bin. These are unitless EGRET feature-change scores, not gap sizes or probabilities. “Raw counts” below means a different input scaling, not simply an unadjusted test statistic.

## Observed comparisons

These scores use per-million transition frequencies and the 12-bin view.

{chr(10).join(table)}

The five-profile focal bin contains **f58r, f58v, f67v2, f69r and f70r2**. Golf 5a comprises f67r1, f67v2, f68v2, f69r, f70r2 and fRos. The broader Golf removal additionally removes Stars groups 3c and 4c. The separate A/C/S case uses the original source-section labels, so it is an alternative definition, not a nested extension of the Golf-group removal. [Exact removed IDs and resulting bin membership](results_word_sheet_residual.json) are retained.

The broad exclusions remove the last part of the fitted sequence, shrinking the surviving normalized position range from 0–1 to about 0–0.902. Native scEGRET consequently moves its bin boundaries even though every surviving profile keeps its original pseudotime. Their 12-bin scores compare different mixtures from the narrow removals and should not be interpreted as the isolated contribution of f58.

## The minimum-bin-size rule matters

Removing f58 leaves three profiles in the focal bin. scEGRET requires **more than three**, so it excludes that bin and directly compares the original bins 7 and 9. Removing all five focal profiles produces essentially the same comparison. Conversely, removing the three other profiles while **retaining f58** leaves only two profiles, also excluded, and produces essentially the same comparison again.

Thus this test cannot isolate whether f58 itself is an indispensable bridge. The source method operates on sufficiently populated bins; it does not evaluate a two-page bridge directly. An excluded bin is not equivalent to physical absence of every text profile it contained.

## What the calibrated comparisons say

The following are targeted simulation-tail comparisons, adjusted over the H-region boundaries in all removal cases, resolutions and scales, and across the two feature channels. Smaller values mean the surviving join is unusually contrasting under that particular reference model. Values at or below 0.05 are flags, not historical probabilities. For the baseline's two boundaries the smaller tail is shown.

{chr(10).join(tails)}

For the per-million inputs shown above, none of these joins passes the adjusted screen under the empirical sheet model. Under the independent-profile model, some narrow-removal joins do. Importantly, this also happens in the control **retaining f58 but removing its three bin companions**. Because the source's occupancy rule then excludes f58 itself, the test cannot identify f58 specifically as the material preventing a break.

A lower observed EGRET score can nevertheless become more unusual under a model: the new comparison averages larger neighboring bins, so a model that treats profiles as independent expects much less variation between those means. The empirical residual model retains broader feature covariance and sheet variation, giving a less decisive result. The earlier calibration's mismatch in feature detection under the independent model remains relevant; that model cannot simply be selected because it gives the stronger finding.

### Raw counts do produce flags under sheet resampling

The following table reports the smallest adjusted tail across the three bin resolutions for raw counts. The adjustments already include scanning those resolutions. “Whole-path” additionally scans **every manuscript boundary** in every case, rather than only the H-region joins.

{chr(10).join(rawtails)}

For example, removing f58 gives a raw word-feature flag at 14 bins (targeted tail 0.004; whole-path tail 0.042). Removing f58 and all Golf 5a also passes the whole-path raw word screen (0.024). These are actual positive screening results and should not be omitted. They do not persist under per-million normalization. Raw input reflects absolute transition totals and gives longer profiles greater influence on bin means; it tests a different mixture of text amount and composition. That makes the scale dependence a material limitation, not grounds for silently choosing either result. The source method supports both inputs; the earlier study designated per-million frequencies as primary for the text comparison.

All views and their actual bin members are retained in the JSON results, including findings outside H. The unchanged 0.12 reference bandwidth is the primary setting from the prior calibration; other smoothing bandwidths were not rerun for this quick omission test.

## Does this support partially missing material around f58?

**Partial loss around a surviving f58 is compatible with these results, but it is not demonstrated by them.** The calculation does not estimate the quantity of missing text. It also does not prove that the gradient lacks a break: the previous controls showed weak detection under empirical sheet noise.

The practical finding is that the currently observed regions on either side can be compared without a larger raw feature jump when f58 is removed. Whether they belong to one historical progression remains a modeling question. This screen neither establishes an actual historical connection nor shows that it is impossible.

## Method and reproduction

We reused the frozen bandwidth-0.12 generators and the exact existing deterministic continuous-control streams. Because synthetic matrices had not been retained, **7,596 control datasets were regenerated**, each scored under all seven cases. These are paired replays, not seven independent sets of controls. Each channel/model uses 400 standardization draws, 999 reference draws and 500 held-out draws. The full-data means and standard deviations reproduce the earlier calibration exactly.

The same omissions are applied to observed and simulated data. Feature detection, normalization, native binning and occupancy exclusions are recomputed after each deletion. This asks whether the remaining join is unusually contrasting **relative to a smooth reference with the same omissions**. It conditions on the missing-data pattern rather than estimating whether or how much material was historically lost. The reference feature curve and noise estimates stay fixed; the original Golf ordering is not refitted.

Every observed removal calculation was checked against the independently R-validated EGRET implementation, with maximum discrepancy below 10⁻¹². This is the same published feature-change statistic as the main study, with the manuscript-specific calibration documented in the [parent report](../report.html). No pooled-bifolium deletion is substituted for removing f58 pages: the H pool also contains f65v.

From this directory:

```
python run_omissions.py
python build_report.py
```

Numerical results: [word / sheet](results_word_sheet_residual.json), [line / sheet](results_line_sheet_residual.json), [word / independent](results_word_dirichlet_multinomial.json), [line / independent](results_line_dirichlet_multinomial.json).

[Return to the calibration report](../report.html) · [Main study](../../report.html)
'''
    (HERE/'RESULTS.md').write_text(md,encoding='utf-8')
    html=gallery.doc('f58 and astro/stars removal tests',MarkdownIt('commonmark',{'html':False}).enable('table').render(md),True)
    (HERE/'report.html').write_text(html,encoding='utf-8')
    targets=[]
    for link in re.findall(r'(?:href|src)="([^"]+)"',html):
        if re.match(r'https?://|#',link):continue
        assert (HERE/link.split('#')[0]).is_file(),link
        targets.append(link)
    (HERE/'presentation_verification.json').write_text(json.dumps(dict(status='passed',local_links=len(targets),targets=targets),indent=2)+'\n')
    print(json.dumps(dict(status='complete',links=len(targets))),flush=True)


if __name__=='__main__':main()
