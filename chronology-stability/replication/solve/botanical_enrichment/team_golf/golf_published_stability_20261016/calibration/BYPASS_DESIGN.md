# Supplementary calibration of the two-sided excursion diagnostic

Specified after inspecting some adjacent-score controls, before running any bypass simulations. The original report's “pocket” language depended on both adjacent EGRET scores exceeding the score directly between the two outer bins. Testing one adjacent score at a time does not fully test that assertion.

For each retained consecutive triplet (a,b,c), define `B = min(E(a,b), E(b,c)) - E(a,c)`. Positive B is the previously described two-sided contrast; its magnitude measures the weaker side's excess over the direct outer-bin contrast. This is our diagnostic built from published EGRET scores, **not a statistic or discontinuity criterion introduced by Zheng et al.** A smoothly bending path or a noisy middle bin can also produce positive B.

Repeat continuous-control calibration for all cohort/channel combinations at bandwidths 0.08, 0.12 and 0.20 under the sheet-residual model, and bandwidth 0.12 under the Dirichlet–multinomial sensitivity model. Use 400/999/500 independent simulations for standardization/calibration/held-out false alarms. Scan all triplets, all three bin counts, both scales. The target H diagnostic requires the f58/H bin to be the middle bin. Apply a factor of eight to scan tails: four cohort/channels times two statistic families (adjacent-change and supplementary bypass). The adjacent score remains the primary test. Do not select whichever diagnostic is more favorable and call it a preplanned published test.

The generator, source score, positions and event totals are unchanged. Seed streams have a separate namespace. No rerouting or historical branch inference is performed. The direction of time still cannot be determined from these absolute contrasts.

Additional resource estimate: 30,384 streamed synthetic datasets, under 4 MB retained results; no dependency downloads.
