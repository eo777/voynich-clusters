# Removing f58 and the astro/stars material

**Completed 16 October 2026.** Seven cases including the full-data baseline; word and line features; 10/12/14 bins; normalized and raw counts; two previously calibrated noise models. No new curve fits. [Design](DESIGN.md) · [Code](run_omissions.py) · [Audit](audit.json).

## Answer

**Some versions of the test do flag a break after removing f58, but the finding depends on count scaling and the noise model.** With per-million frequencies and empirical sheet resampling, no tested removal produces a flag. Raw counts do produce H-region flags under sheet resampling, and the independent-profile Dirichlet–multinomial model also flags some normalized comparisons after the narrower removals. This is evidence of a sensitive comparison, not a robust demonstration that f58 is the surviving bridge preventing a discontinuity.

The **per-million feature contrast** actually decreases after removing f58. At 12 bins, the new comparison is approximately **0.168 word / 0.177 line**, versus **0.271 word / 0.279–0.286 line** entering and leaving the original focal bin. These are unitless EGRET feature-change scores, not gap sizes or probabilities. “Raw counts” below means a different input scaling, not simply an unadjusted test statistic.

## Observed comparisons

These scores use per-million transition frequencies and the 12-bin view.

| Removal | Profiles removed | Word score at H join | Line score at H join |
|---|---:|---:|---:|
| Full data | 0 | 0.271 / 0.270 | 0.279 / 0.286 |
| Remove f58r/v | 2 | 0.168 | 0.177 |
| Remove the five-profile focal bin | 5 | 0.168 | 0.177 |
| Remove f58 and Golf 5a | 8 | 0.182 | 0.186 |
| Remove f58 and Golf 5a, 3c, 4c | 34 | 0.257 | 0.287 |
| Remove source sections A, C, S | 38 | 0.154 | 0.159 |
| Retain f58; remove its three bin companions | 3 | 0.170 | 0.177 |

The five-profile focal bin contains **f58r, f58v, f67v2, f69r and f70r2**. Golf 5a comprises f67r1, f67v2, f68v2, f69r, f70r2 and fRos. The broader Golf removal additionally removes Stars groups 3c and 4c. The separate A/C/S case uses the original source-section labels, so it is an alternative definition, not a nested extension of the Golf-group removal. [Exact removed IDs and resulting bin membership](results_word_sheet_residual.json) are retained.

The broad exclusions remove the last part of the fitted sequence, shrinking the surviving normalized position range from 0–1 to about 0–0.902. Native scEGRET consequently moves its bin boundaries even though every surviving profile keeps its original pseudotime. Their 12-bin scores compare different mixtures from the narrow removals and should not be interpreted as the isolated contribution of f58.

## The minimum-bin-size rule matters

Removing f58 leaves three profiles in the focal bin. scEGRET requires **more than three**, so it excludes that bin and directly compares the original bins 7 and 9. Removing all five focal profiles produces essentially the same comparison. Conversely, removing the three other profiles while **retaining f58** leaves only two profiles, also excluded, and produces essentially the same comparison again.

Thus this test cannot isolate whether f58 itself is an indispensable bridge. The source method operates on sufficiently populated bins; it does not evaluate a two-page bridge directly. An excluded bin is not equivalent to physical absence of every text profile it contained.

## What the calibrated comparisons say

The following are targeted simulation-tail comparisons, adjusted over the H-region boundaries in all removal cases, resolutions and scales, and across the two feature channels. Smaller values mean the surviving join is unusually contrasting under that particular reference model. Values at or below 0.05 are flags, not historical probabilities. For the baseline's two boundaries the smaller tail is shown.

| Case | Word: sheet model | Line: sheet model | Word: independent model | Line: independent model |
|---|---:|---:|---:|---:|
| Full data | 1.000 | 0.486 | 0.936 | 0.224 |
| Remove f58r/v | 1.000 | 0.640 | 0.010 | 0.012 |
| Remove the five-profile focal bin | 1.000 | 0.602 | 0.010 | 0.008 |
| Remove f58 and Golf 5a | 1.000 | 0.504 | 0.004 | 0.004 |
| Remove f58 and Golf 5a, 3c, 4c | 1.000 | 1.000 | 1.000 | 0.244 |
| Remove source sections A, C, S | 1.000 | 0.882 | 0.156 | 0.178 |
| Retain f58; remove its three bin companions | 1.000 | 0.660 | 0.010 | 0.012 |

For the per-million inputs shown above, none of these joins passes the adjusted screen under the empirical sheet model. Under the independent-profile model, some narrow-removal joins do. Importantly, this also happens in the control **retaining f58 but removing its three bin companions**. Because the source's occupancy rule then excludes f58 itself, the test cannot identify f58 specifically as the material preventing a break.

A lower observed EGRET score can nevertheless become more unusual under a model: the new comparison averages larger neighboring bins, so a model that treats profiles as independent expects much less variation between those means. The empirical residual model retains broader feature covariance and sheet variation, giving a less decisive result. The earlier calibration's mismatch in feature detection under the independent model remains relevant; that model cannot simply be selected because it gives the stronger finding.

### Raw counts do produce flags under sheet resampling

The following table reports the smallest adjusted tail across the three bin resolutions for raw counts. The adjustments already include scanning those resolutions. “Whole-path” additionally scans **every manuscript boundary** in every case, rather than only the H-region joins.

| Case | Word: targeted tail | Word: whole-path tail | Line: targeted tail | Line: whole-path tail |
|---|---:|---:|---:|---:|
| Full data | 1.000 | 1.000 | 0.660 | 1.000 |
| Remove f58r/v | 0.004 | 0.042 | 0.016 | 0.182 |
| Remove the five-profile focal bin | 0.004 | 0.038 | 0.016 | 0.168 |
| Remove f58 and Golf 5a | 0.002 | 0.024 | 0.014 | 0.146 |
| Remove f58 and Golf 5a, 3c, 4c | 1.000 | 1.000 | 0.870 | 1.000 |
| Remove source sections A, C, S | 0.004 | 0.044 | 0.008 | 0.076 |
| Retain f58; remove its three bin companions | 0.002 | 0.022 | 0.010 | 0.098 |

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
