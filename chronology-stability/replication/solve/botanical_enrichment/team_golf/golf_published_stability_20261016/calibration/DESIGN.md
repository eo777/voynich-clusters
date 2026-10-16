# scEGRET manuscript calibration — design fixed before simulation

16 October 2026. This follow-up supplies the controls missing from the original scEGRET screening. The original 648 R results remain unchanged.

## Question and relationship to the paper

Zheng et al. (2023), https://doi.org/10.1093/bib/bbad356, validate the EGRET feature-change score using simulated continuous trajectories, removal of intermediate states, and random downsampling controls (Methods; Figures 1–2). We repeat that validation logic with the observed manuscript dimensions, transition totals, pseudotimes and binning. This is a text-specific calibration experiment, not a method or statistical significance test supplied by the authors. The published 0.5 reference flag is not used as a universal threshold.

The original Golf order is held fixed. This tests feature continuity conditional on that order; it does not calibrate uncertainty from learning the order, identify a historical missing bifolium, or test branch topology. No finite set of observations can rule out every continuous curve: an arbitrarily fast smooth excursion can interpolate any finite data. Consequently the smoothness and noise assumptions below are part of every inference.

## Generators

Use the actual 207 paragraph profiles and 50 independently analysed pooled profiles; word and line transition channels are separate. Preserve all observed transition totals and positions. Normalize position to [0,1]. Fit a Gaussian-kernel mean of per-profile transition proportions, with bandwidths 0.08, 0.12 and 0.20 of the full path, fixed in advance. All kernels are smooth; no intermediate state is removed in baseline controls. Fit includes H, rather than estimating a null after selectively discarding the suspect region. Bandwidth 0.12 is the primary reference, with both other choices reported, not selected by significance.

Two noise models expose the consequences of assumptions:

1. **Sheet residual bootstrap (primary).** Compute square-root proportion residuals against a leave-one-physical-sheet-out smooth prediction. Center the residual library. Resample whole donor sheets; map within-sheet ordered profiles by relative rank when sheet sizes differ. Add residuals to the square root of the smooth target mean, truncate negative roots at zero, square and renormalize. Match observed event totals by deterministic largest-remainder rounding. This retains empirical feature covariance and sheet dependence approximately. Adjust residual amplitude for donor/target count totals using the overdispersion moment estimate below. No second multinomial draw is added: residuals already contain sampling noise. Truncation, donor mapping, estimated residuals and global noise exchangeability are limitations, checked through diagnostics.
2. **Dirichlet–multinomial (sensitivity).** Estimate a single excess-variance parameter per cohort/channel/bandwidth from leave-sheet-out squared proportion residuals, including lack of smooth fit as noise. Draw each profile's probabilities from a Dirichlet distribution around its smooth mean, then multinomial counts at its observed total. This reproduces count noise and overdispersion but assumes independent profiles and a restricted covariance structure. Its results alone cannot establish a robust finding.

No generator is a known ground truth for the manuscript. Compare simulated and observed feature detection, zero fraction and squared residual magnitude to expose obvious mismatch. Smooth generators may retain gradual excursions; continuity does not mean a straight line or a constant rate of change.

## Exact statistic and multiplicity

Use scEGRET's published source calculation unchanged: retain features detected in >10% of profiles, divide each bin-mean log1p feature contrast by its whole-dataset log1p mean, average over retained features; retain bins with >3 observations. Bins use R cut semantics. Primary input is counts per million; raw counts are a scale sensitivity control. Page bin counts are 10,12,14; pooled counts are 6,8,10. Retain original bin labels and sparse-bin exclusions.

The production implementation is the independent NumPy translation already validated against all 648 original-R runs to floating-point precision. Before calibration, recheck all original runs. Also compare a fast vectorized simulation scorer against that translation, including sparse-bin cases. This accelerates simulation without replacing the EGRET statistic. The supplementary two-sided bypass pattern is reported descriptively, not relabelled as the authors' test.

For each generator/bandwidth/cohort/channel use disjoint random streams: 400 continuous simulations estimate each boundary's mean and standard deviation; 999 further simulations calibrate the largest standardized adjacent EGRET score across every retained boundary, all three bin resolutions and both count scales. A third set of 500 continuous simulations estimates the out-of-sample false-alarm rate at the resulting 95% threshold. This accounts for scanning boundaries and bin choices within the channel. Apply a conservative factor of four across the two channels and two cohorts for a whole-study tail fraction; no independence assumption is required for that factor. Generator/bandwidth choices are sensitivity models, not six opportunities to select a favourable result. All observed boundary scores receive both local and scan-adjusted simulation tails. H was selected after inspecting the manuscript, so the scan-adjusted result governs any claim of detection.

Finite simulation tails are (1 + exceedances)/(999 + 1). Report held-out false-alarm rates with binomial Wilson intervals. These are conditional simulation references, not calibrated historical probabilities or exact p-values after estimating a null and trajectory from the same data.

## Missing-intermediate positive controls

At primary bandwidth 0.12, for both noise models, cohorts and channels, remove contiguous intervals centered at normalized path positions 0.35, 0.65 and 0.85, with widths 0.10, 0.20 and 0.30 (intervals clipped to [0,1]; endpoints preserved). Use the resulting number of removals to generate randomly thinned continuous calibration controls, preserving endpoints. For each condition, 399 random-thinning simulations set the 95% reference threshold for the maximum adjacent CPM EGRET score at the middle bin count (12 pages, 8 pools). In 300 additional paired simulations compare deliberate-gap and same-size random deletion. Recompute feature detection, binning and sparse-bin exclusions after every deletion. Report detection rate, random-removal false-alarm rate, Wilson intervals, retained sample size and bin counts. These powers are for the specified synthetic gaps, not a probability of detecting every historical loss. Fixed original pseudotime is used; trajectory reconstruction is not repeated. This is the screening stage of the paper's validation, not a repeat of all five trajectory-inference pipelines.

## Attribution checks

For the observed H region, recompute the original score after excluding f58r/v, after excluding the other three profiles in its 12-bin neighborhood, and after excluding the complete H sheet. Preserve original positions; document loss of bin eligibility. No apparent disappearance caused by the >3 rule counts as evidence that a particular profile caused the signal. List influential transition features and compare the H profiles with their neighbors directly as a descriptive check.

## Reproducibility and resources

Fixed seed 20261016, deterministic NumPy SeedSequence streams. Save configuration, input/source hashes, compact per-boundary results and simulation reference distributions, controls and diagnostics. Do not retain synthetic feature matrices. Estimated peak memory below 350 MB and additional disk below 30 MB; no model arrays or new dependencies required. Run the workspace storage check after completion. Private/unpublished material is not an input.
