# Golf ordering uncertainty and H-position sensitivity

**Completed 16 October 2026.** The study runs the original ElPiGraph fitting implementation, Cassidy's derivative-based continuation calculation, and the original scEGRET assessment functions, with the manuscript-specific adaptations documented below. All numerical audits passed.

The main findings are that **3a/3b remains competitively unstable**, observed H lies very near an equal-score frontier, and pooling can hide feature discontinuities rather than resolve them. Bio/Stars alternatives persist, but the best recovered elastic scores usually favor Stars before Bio. scEGRET adds useful continuity information while exposing a serious resolution limitation: pooled Bio and Stars often fall in the same bin, so its aggregate score cannot settle their internal order.

[Open the figures and downloads](index.html). The two requested replacements are [overall ordering sensitivity](ordering_ensemble.png) and the [H-position map at lambda 0.01](cassidy_H_regions.png). The additional [scEGRET figure](scegret_original_golf.png) shows the original Golf route's feature continuity. The [whole-sheet sampling control](ordering_ensemble_sheet_control.png) is provided separately.

## Purpose

The two requested outputs are (1) an overall ordering-sensitivity chart covering individual paragraph profiles and pooled bifolia, including B-region alternatives, and (2) an MDS1/2 and MDS1/3 map of how the position of observed H:58–65 affects preference between Herbal 3a and 3b. Here “Herbal A/B” refers to Golf 3a/3b; it does not relabel Currier A/B.

The original Golf sequence used an eleven-knot cubic spline. This study uses the published elastic-principal-curve model to examine alternative orderings. Its numerical variability belongs to this model. Independent historical evidence is not an input to the scores. No private or unpublished research is used.

| Question | Published calculation used | What is varied | Output |
|---|---|---|---|
| Which fitted positions and orders recur? | ElPiGraph principal-graph ensembles | Random observation inclusion, initial route and penalty | Position frequencies and joint order outcomes |
| How does moving H change competing fitted solutions? | Cassidy implicit-function continuation | Two coordinates of the observed H pool | Score-preference frontier and route dependence |
| Where does a proposed route join sharply different text profiles? | scEGRET continuity assessment | Proposed order, bin width, word/line features and count scaling | Adjacent-bin and all-bin feature changes |

## Published techniques actually implemented

### Random principal-graph ensembles: ElPiGraph

Albergante et al. (2020), section 2.10 and Figure 1C, describe repeated fits to random 90% subsets, with 100 realizations, to inspect which regions of a principal graph vary. The same section describes varying initializations and parameters as additional ensemble constructions. We reproduce that random-resampling procedure using ElPiGraph 0.3.4 for fresh graph growth, rather than reusing the earlier leave-one-bifolium-out runs.

The Python package implements the 90% fraction as independent Bernoulli inclusion, so individual subsample sizes vary. A separate physical-sheet control includes or excludes all paragraph profiles from a bifolium together. The block control changes the sampling unit, because profiles on a sheet are related; it is explicitly an adaptation rather than an additional published example.

The native wrapper subsets observations without subsetting supplied point weights. We therefore generate exactly the same Bernoulli masks outside the wrapper, subset the weights consistently, and call the unmodified graph-growth algorithm. An unweighted native-wrapper comparison reproduces the explicitly masked fit to maximum node difference 9.71×10⁻¹⁷. The wrapper's extra averaged graph is not counted as another random realization. `AvoidResampling` controls initialization-point reuse, not data subsampling; the earlier fits had no random data subsampling because they used `ProbPoint=1` and `nReps=1`.

There are 100 independent masks per sampling scheme, eight starts per mask and two penalties. Four starts cover the joint 3a/3b × Bio/Stars alternatives recovered by the earlier B-permutation search; their mirrored early-A prefixes supply the other four. The representatives are selected separately by cohort and penalty. This is a balanced exploration of sampling, starts and regularization, not a pure paired experiment isolating lambda. Eight starts share each mask and must not be treated as eight independent data samples. Main point ensembles comprise 3,200 fits; the paragraph whole-sheet control adds 1,600.

Every retained observation has a position histogram on the normalized curve. The figure combines starts and penalties with equal design weight and preserves multiple modes. These position-frequency strips are our readout of the published ensemble procedure, not a confidence-interval formula established by that paper. The physical-sheet control is shown separately. Group-median and individual-pair order frequencies distinguish actual reversals from coherent movement along the curve. All-reference projections are retained as diagnostics; omitted observations do not enter the main retained-only histograms.

Sources: [paper, section 2.10](https://pmc.ncbi.nlm.nih.gov/articles/PMC7516753/), [published Python implementation](https://github.com/j-bac/elpigraph-python).

### Implicit-function predictor and corrector: Cassidy

Cassidy (2023), section 2.2, equations 7–8, obtains the sensitivity of an optimum to changed data from the inverse parameter Hessian and the mixed data/parameter derivative. Remark 1 permits a general twice continuously differentiable objective; the method is not restricted to a biological likelihood.

For fitted vertex coordinates Y and data X, the elastic objective within fixed nearest-node assignments is

`G(Y,X) = sum_i w_i ||x_i - y_a(i)||² + lambda ||B Y||² + mu ||C Y||²`.

Here B takes adjacent-vertex differences, C takes each interior vertex minus the mean of its two neighbors, and a(i) is the assigned nearest vertex. Define `A = diag(assigned masses) + lambda B.T B + mu C.T C`. The Hessian for each coordinate is `2A`. Moving H in that coordinate has mixed derivative `-2 w_H e_a(H)`. Cassidy's equation 8 gives

`dY / dx_H = solve(A, w_H e_a(H))`.

This predictor is exactly affine as long as assignments do not change. The implementation detects nearest-node equalities for every observation and competing vertex along a proposed move. At an assignment boundary it corrects the fitted solution and recomputes the derivative. The smooth theorem is not asserted across those events. Endpoint checks compare straight, x-then-y, y-then-x and return paths; numerical derivative and event-step checks are retained.

Source: [Cassidy (2023), sections 2.1.1–2.3](https://link.springer.com/article/10.1007/s11538-023-01200-0). Exact implementation and diagnostics are described in [CASSIDY_DESIGN.md](CASSIDY_DESIGN.md).

### Continuity screening: scEGRET

scEGRET belongs in this analysis because it evaluates continuity along an already proposed trajectory. Being a screening method rather than a trajectory-fitting method is not a reason to exclude it. The original `egret.flap` and `egret.fly` R functions are executed unchanged from commit `93517e5ce3f09ccdd3c0d604d0a1b53e70c5cc70`. We use a local R/WebR interpreter; the package's plotting wrappers and unrelated Seurat/Monocle convenience dependencies are not needed to execute these computations. Runtime versions, source files and checksums are retained.

The text adaptation treats each transition type as a feature and each paragraph or pool as an observation. Word and line transitions are evaluated separately. Counts per million transitions within each profile are the primary input, analogous to the TPM input the package accepts; raw event counts provide a text-length sensitivity control. These are declared text inputs, not a claim that expression thresholds transfer to manuscript analysis.

The actual source calculation for retained bins k and l is

`E(k,l) = mean_g[ abs(log1p(mean_k C_g) - log1p(mean_l C_g)) / log1p(mean_all C_g) ]`.

Features must occur in strictly more than 10% of observations. Equal-width pseudotime bins must contain strictly more than three observations. The latter differs from the paper's brief prose; the whole-profile mean normalization is also more specific than its displayed formula. We follow the executable source and disclose both details. Bin averages give observations equal weight; they do not use the curve's equal-physical-sheet weighting. Cohorts have different retained feature sets, so absolute scores should not be compared as if based on identical assays.

We examine 27 saved orders per cohort: eighteen original elastic fits, eight B-permutation representatives covering both Herbal and terminal directions at both penalties, and the original Golf spline. Paragraph resolutions are 10/12/14 bins; pooled resolutions are 6/8/10. Each is evaluated for word/line × normalized/raw inputs, giving 648 unique retained evaluations. The middle resolutions follow the approximate Rice-rule guidance in the source. Sparse bins remain visible with original bin numbers; connecting surviving bins across a discarded bin is explicitly marked. Adjacent scores and all-bin-pair scores reproduce the published continuity-screening approach. Mean and maximum adjacent scores are our descriptive summaries, not calibrated chronology likelihoods. We do not transfer the biological example's 0.5 reference threshold.

scEGRET examines the same observed features used to build the geometry. It is a complementary criterion, not independent historical evidence. A jump may reflect missing text, section/topic differences, scribal differences or short samples. A smooth route can still be historically wrong: mixing heterogeneous observations into similar bin averages can also lower between-bin change. We therefore retain bin composition and geometric-fit context rather than minimizing scEGRET as a stand-alone ordering objective. Whole-route reversal leaves the absolute-change calculation unchanged, so the method cannot establish historical direction by itself. Moving H in MDS space does not supply a hypothetical text profile; we do not invent such profiles for scEGRET. Although the original function also performs Wilcoxon comparisons, those p-values are not used as manuscript significance tests here.

Sources: [Zheng et al. (2023), methods and Figure 2](https://academic.oup.com/bib/article/24/6/bbad356/7325910), [original scEGRET source](https://github.com/qingshanni/scEGRET).

## Data and readout

The frozen data are 207 paragraph-text profiles on 50 physical bifolia and their 50 pooled profiles. The 223-unit whole-manuscript display also includes other text-unit kinds; the frozen chronology-test cohort selects `kind=para`. Distances combine the existing word- and line-transition divergences. Classical MDS treats those divergences as squared-distance input and retains all positive dimensions: 196 for paragraph profiles, 49 for pools. These two native geometries are separate. Both model families use 28 vertices and harmonicity penalty mu=0.1; the ordering ensemble includes lambda=0.01 and 0.05, while the main H map uses lambda=0.01 only.

Fits minimize distance to the nearest vertex plus the elastic penalties. Displayed ordering is calculated afterward by projection to the closest curve segment and normalized arc length `t=s/L`. Thus vertex-assignment stability and ordering stability are distinct. Positions are relative to each curve, not dates or a common elapsed-time scale. Final whole-curve orientation is normalized; initialization does not impose the final internal order.

Group-order statements compare median projected positions. “3a before 3b” does not mean that every 3a page precedes every 3b page; groups may interleave. Individual-pair counts are retained separately. Pooling changes the observations, distances, geometry and sometimes Golf labels for mixed sheets, so it can change the preferred route rather than simply narrow uncertainty around the page-level route.

## How to read the figures

### Ordering ensemble

The figure shows empirical position frequencies from random-resampled graph ensembles, retaining separate modes rather than fitting bell curves or averaging competing routes. Each horizontal strip is one observation, with position from 0 to 1 on the horizontal axis. Darker color means more runs placed that observation in that interval. Several dark patches mean several recurring positions. Rows are grouped by Golf label for readability; row order is not a recovered chronology. A shared intensity scale allows comparisons between rows, with saturation above 20% per interval. Page and pooled results remain separate. Starts and penalties receive declared design weights. Frequencies describe this chosen ensemble; they are not posterior probabilities of historical positions or confidence intervals. Joint group-order outcomes remain available in the report because overlapping marginal position distributions alone cannot prove a reversal.

### H-position map

Only the two displayed coordinates of H:58–65 change. H pools the observed f58r, f58v and f65v paragraph profiles. Its other coordinates, all other profiles, weights, MDS basis, topology and penalties remain fixed. The plotted frontier compares minimum recovered scores for the two strict Herbal directions among final-1a-first solutions. Colors describe preference within the explored families and continuation paths. They do not locate missing text or prove that a historical connection exists. Path disagreement and assignment events are reported separately from equal-score crossings.

This map perturbs the observed H pool; it does not reconstruct the missing inner-H bifolia. A hypothetical MDS coordinate need not correspond to a realizable transition profile. The earlier missing-material insertion experiments remain scenarios, not estimated missing texts. scEGRET can flag observed discontinuities consistent with absent intermediates, but it cannot place unobserved bifolia without further feature assumptions or evidence.

## Completed H-continuation results

Ten fitted families are continued over each 25×25 grid along three routes: straight, horizontal-then-vertical and vertical-then-horizontal. Straight routes are also reversed. A score frontier is where competing families exchange preference; it is not automatically a loss of local stability or evidence of a historical branch.

| Result | MDS 1/2 | MDS 1/3 |
|---|---:|---:|
| Grid cells preferring 3a before 3b | 301/625 | 175/625 |
| Cells changing preference versus the earlier independent-refit map | 127/625 | 127/625 |
| Cells where preference depends on the continuation route | 35/625 | 32/625 |
| Cells preferring Stars 4c before Bio 4a | 625/625 | 617/625 |

At observed H the recovered Herbal advantage is only **0.01686% in favor of 3a first**. It remains effectively near the equal-score boundary, not a decisive ordering result. The eight Bio-first cells in MDS1/3 also show why the terminal ordering should not be described as universally invariant to this perturbation. Grid-cell proportions depend on the plotted range and are not probabilities of where H belongs.

Each individually tracked seed retains its focal group-order signs through these moves; coordinates and costs change, and the preferred seed can change. Some round trips return to different fine solutions. A separate check of 36 selected excursions/returns leaves the native-H Herbal advantage unchanged but reduces the recovered Stars-before-Bio advantage from 1.464% to 0.971%. This is evidence of local optimization/path sensitivity, not an exhaustive search over all possible histories.

The minimum fixed-assignment Hessian eigenvalue encountered is positive (0.02340). Positive curvature is expected from this connected quadratic model with positive assigned mass; it supports local uniqueness within those assignments, not global uniqueness or a historical connection across missing material. Sixty derivative checks have maximum relative error 1.97×10⁻¹⁰. Thirty-six smaller event-crossing checks recover identical endpoints. An independent audit checks 30 derivatives, 24 event roots, 72 outward paths and 24 returns; objective differences are at most 2.08×10⁻¹⁷. See [continuation summary](cassidy_h_summary.json), [independent audit](cassidy_independent_audit.json) and [return check](cassidy_return_check.json).

## Ensemble results

All 4,800 fresh fits completed and their final refinements converged. The complete accounting audit verifies the masks, retained weights, source hashes and objective components. Convergence means convergence to a recovered local solution, not proof of the global optimum.

### Herbal reversal persists under random resampling

For each random mask and penalty, we compare the best recovered candidates on that same retained data. The following counts restrict candidates to those whose final median order is 1a before 1b. This respects the requested early-order preference as a readout, without pretending initialization imposed a hard constraint. Each row has 100 independent masks; the eight starts within a mask are competing solutions, not independent samples.

| Data and sampling unit | Lambda | 3a before 3b | Bio 4a before Stars 4c | Stars 3c before Bio 4a | Stars 3c before Stars 4c |
|---|---:|---:|---:|---:|---:|
| Pages, individual profiles | 0.01 | 44/100 | 2/100 | 98/100 | 98/100 |
| Pages, whole physical sheets | 0.01 | 38/100 | 2/100 | 98/100 | 98/100 |
| Pages, individual profiles | 0.05 | 24/100 | 13/100 | 87/100 | 87/100 |
| Pages, whole physical sheets | 0.05 | 25/100 | 24/100 | 76/100 | 76/100 |
| Pooled bifolia | 0.01 | 58/100 | 1/100 | 99/100 | 99/100 |
| Pooled bifolia | 0.05 | 23/100 | 0/100 | 100/100 | 100/100 |

There are no median-order ties among these winners. These are finite-search frequencies under the declared resampling scheme, not chances that a historical order is true. The 3a/3b direction is competitively unsettled: both directions win substantial numbers of samples, and pooled results change their majority between penalties. Pooling does not simply tighten a fixed order.

For scale, define the directional score gap as `100 × (best reversed score − best first-listed score) / min(the two scores)`. Positive gaps favor the first-listed group. Comparisons use only an identical mask and lambda. The 5th/50th/95th percentiles for **3a before 3b**, among final-1a-first candidates, are:

| Data and sampling unit | Lambda | 5th percentile | Median | 95th percentile |
|---|---:|---:|---:|---:|
| Pages, individual profiles | 0.01 | −0.297% | −0.047% | +0.618% |
| Pages, whole physical sheets | 0.01 | −0.504% | −0.136% | +0.489% |
| Pages, individual profiles | 0.05 | −0.752% | −0.182% | +0.386% |
| Pages, whole physical sheets | 0.05 | −0.762% | −0.199% | +0.355% |
| Pooled bifolia | 0.01 | −0.990% | +0.100% | +1.836% |
| Pooled bifolia | 0.05 | −0.721% | −0.253% | +0.856% |

These distributions cross zero in every setting. They quantify competition between recovered routes; the percentages are elastic-objective differences, not probabilities or chronological effect sizes. Different routes can change several B relationships together. They do not isolate the causal effect of a single swap. Starting representatives also differ by lambda, so the cross-lambda contrast includes both regularization and the enriched search design.

### Bio/Stars alternatives persist, with a different score pattern

Stars 4c before Bio 4a has the lower recovered score in most samples, even though Bio-first alternatives frequently persist as local solutions from other starts. At lambda 0.05, including or omitting physical sheets together raises the Bio-first winner count from 13/100 to 24/100. This descriptive difference shows that the sampling unit matters; it is not a formal significance test. At lambda 0.01 the corresponding count is 2/100 under both page schemes.

The wider sheet-control positions are especially visible for Stars 4c at lambda 0.05: the median across pages of the empirical 5th-to-95th position span changes from 0.094 to 0.295 on the normalized curve. For f58r/v at lambda 0.01, those spans change from approximately 0.09–0.10 under page sampling to approximately 0.36 under sheet sampling. These spans describe retained-only positions across all starts. A wide position span can reflect a whole block moving together and is not itself proof of a relative-order reversal.

The eight starts deliberately cover the joint 3a/3b and Bio/Stars alternatives; they do not exhaust all 120 possible five-group B orders. Stars 3c relationships are reported only for the coverage actually explored. Some 3c and Bio/Stars comparisons select the same competing candidates, explaining repeated winner counts; they are not independent corroborating tests.

### Why the combined chart includes more than the winning route

The position strips give equal design weight to the eight starts and two penalties. They show where solutions recur and preserve alternative modes, including locally persistent solutions with somewhat worse scores. Consequently, a mode's height partly reflects our deliberate balance of starting routes. It is neither a random-initialization basin volume nor an objective-weighted posterior. The per-mask winner table above supplies the missing score context; neither display replaces the other. Averaging the modes into one bell curve would put apparent mass between distinct solutions and hide the alternatives the study was designed to expose.

The earlier figures and scientific results remain preserved for comparison. Detailed counts, pairwise orders, omitted-direction availability and per-start objective excesses are retained in [ensemble findings](ensemble_findings.json), [ensemble summary](ensemble_summary.json), [within-mask competition](ensemble_competition.json) and the [accounting audit](ensemble_audit.json).

## scEGRET results

### H is a distinct pocket in some page-level views

Along the original Golf route, the 12-bin page calculation puts f58r and f58v with f67v2, f69r and f70r2 in bin 8. Both transitions into and out of that bin are larger than the direct feature contrast between its two neighbors:

| Primary transition channel | Bin 7 → 8 | Bin 8 → 9 | Bin 7 → 9, bypassing bin 8 |
|---|---:|---:|---:|
| Word | 0.27096 | 0.27023 | 0.16953 |
| Line | 0.27901 | 0.28606 | 0.17686 |

The bypass is an existing all-bin-pair scEGRET score, not a refit after deleting H. This two-sided contrast occurs in four of the six primary word/line × resolution views, and four of six raw-count controls. It holds at 10 and 12 bins; changing the neighboring membership at 14 bins removes this particular pattern. Thus the feature evidence supports a locally distinct pocket in those views, rather than an unambiguous smooth intermediate. It does not establish whether the cause is missing text, topic, scribal differences or short samples, and it concerns the entire mixed bin rather than f58 alone.

Pooling changes the test. At eight bins, H:58–65 shares a bin with four 3a pools, two 3b pools and one 5a pool. Neither channel shows the same two-sided contrast. Across the six pooled primary views, none of the four eligible comparisons satisfies it; the 10-bin comparisons are unavailable because the H-containing bin has only two observations and is excluded by scEGRET. Raw controls have the same pattern. Averaging H into a mixed bin can therefore obscure the page-level contrast. That is a loss or change of resolution, not evidence that pooling has resolved the historical uncertainty.

### Continuity scores do not select a robust B ordering

Among the four B-family representatives at each penalty, no single candidate has the smallest mean adjacent-bin score in all six primary channel/resolution views, for either cohort. The counts below summarize these six dependent diagnostic views; they are not probabilities. A family label reports two median-order signs and can conceal other changes in the five-group B route.

| Cohort and lambda | 3a first, Bio first | 3a first, Stars first | 3b first, Bio first | 3b first, Stars first |
|---|---:|---:|---:|---:|
| Pages, 0.01 | 2 | 0 | 4 | 0 |
| Pages, 0.05 | 2 | 2 | 2 | 0 |
| Pools, 0.01 | 0 | 2 | 4 | 0 |
| Pools, 0.05 | 0 | 2 | 2 | 2 |

Raw-count controls often change these rankings. More decisively, no pairwise comparison between the pooled B-family candidates has a view in which **both candidates retain any between-bin Bio 4a–Stars 4c resolution**. In the original Golf eight-bin pooled view, all five Bio pools and both Stars 4c pools occupy bin 8: all ten cross-group pairs are within one bin. Both groups are present, but their internal order is invisible to the score. A pure two-pool Stars bin would instead fail the minimum occupancy rule. Pooled aggregate score preferences therefore cannot adjudicate the Bio/Stars reversal here.

The page result also depends on the broader B route, not just a final swap. At lambda 0.01, the descriptive primary winners at 10 and 12 bins follow median order Bio 4a → Stars 4c → Stars 3c → Herbal 3b → Herbal 3a; at 14 bins the last two switch. Raw-count winners change that broader preference. These comparisons are useful diagnostics of feature continuity and its dependence on binning and length scaling, not grounds for replacing the Golf chronology with a uniquely preferred sequence.

### Reading the scEGRET figure

The triangle compares every pair of retained position bins; stronger color means greater feature contrast. The bottom adjacent-bin tracks show where the proposed route encounters those changes, separately for word and line transitions. Occupancy numbers and Golf-colored composition bars explain what was averaged. An asterisk identifies the f58-containing page bin or H:58–65 pool bin. Hatched cells are unavailable because a bin failed the minimum occupancy rule. Dashed connections span an excluded original bin and must not be read as an observed continuous intermediate.

All 648 original-R outputs pass an independent reconstruction of filtering, bin labels, membership, adjacent scores, all-pair scores and top feature contributions. Maximum absolute score discrepancies are 1.33×10⁻¹⁵ for adjacent scores and 4.44×10⁻¹⁶ for all-pair scores. Seven legacy run exports were corrected from compressed factor codes to original bin numbers; memberships and numerical scores did not change. See [scEGRET summary](scegret_summary.json), [independent audit](scegret_independent_audit.json), and the retained [source manifest](scegret_source/source_manifest.json).

## Verification and reproducibility

The principal-graph calculation comprises 4,800 production fits, two native-wrapper equivalence fits and two independent weighted replay fits. Every production fit is accounted for; none failed or was discarded. Both independent replays reproduce the saved full-node hashes exactly, with independent projection differences at most 3.33×10⁻¹⁶. The [replay audit](ensemble_independent_audit.json) is separate from the complete [mask, weight and objective audit](ensemble_audit.json).

The scEGRET release contains 648 unique evaluations from the authors' unchanged R functions, plus runtime versions and source hashes. A probe and redundant intermediate work are not counted as independent evaluations. All 648 summary means and 1,296 pair-bin resolution records also passed independent review. The Cassidy derivative, assignment-event, outward-path and return checks are reported above.

All four PNG figures were visually inspected; their SVG versions parse and match their recorded provenance hashes. Source code, compact result data, feature inputs, captions and reproduction instructions are retained. Temporary research runtimes and duplicate partition checkpoints were removed after verification, with their hashes recorded. The final workspace storage check passed at 18.29 GiB with no prohibited completed-study arrays. See [REPRODUCE.md](REPRODUCE.md) for the execution sequence and [DESIGN.md](DESIGN.md) for the declared choices.

## What these procedures do not establish

Random-resampling ensembles can reveal instability without assigning historical probabilities. Derivative-based continuation can track sensitivity and local solution persistence without making a score frontier a loss-of-stability boundary. Neither procedure proves global optimization, independently dates the text, or discovers branches when constrained to a chain. The previous learned-tree and connectivity studies address topology separately.

The frozen distance representation and MDS basis are not re-estimated inside the ensembles. This study measures fitting and resampling sensitivity conditional on that representation; it does not encompass transcription errors, alternative distance definitions, all missing-text models or every possible fitted route. This study does not claim a jackknife confidence interval, a Bayesian chronology posterior or a temporal Lyapunov exponent.
