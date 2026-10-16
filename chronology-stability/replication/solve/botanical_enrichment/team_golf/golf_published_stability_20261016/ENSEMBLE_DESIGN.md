# Random-resampled principal-curve ensemble

Declared before fitting, 16 October 2026. This implements the random-resampled principal-graph ensemble in Albergante et al. (2020), section 2.10 and Figure 1C: 100 realizations with a 90% observation-inclusion probability. That section also allows multiple initializations and parameter choices. The figure's subsequent optional consensus-graph construction is not used here.

Primary sources: https://www.mdpi.com/1099-4300/22/3/296 and https://github.com/j-bac/elpigraph-python/blob/master/elpigraph/_BaseElPiWrapper.py ; API https://elpigraph-python.readthedocs.io/en/stable/elpigraph.computeElasticPrincipalCurve.html .

## Balanced design

- Native point subsampling: paragraph207 and pool50, 100 Bernoulli(0.9) masks, eight starts, lambda 0.01 and 0.05: 3,200 fresh fits.
- Physical-sheet adaptation: paragraph207 only, 100 Bernoulli(0.9) sheet masks, eight starts, both lambdas: 1,600 fresh fits. Pool observations already are physical sheets, so no duplicate pooled control is needed.
- Total: 4,800 fresh ensemble fits, plus two pilot growth calls that compare one native unweighted package fit against an explicit-mask equivalent. The pilot is excluded from the ensemble.
- Masks are shared across starts and penalties within each cohort/sampler. Sheet masks use the same sheet-level draws as the pool point masks, mapped by bifolium identity. The paragraph point and physical-sheet masks use different fixed random seeds. No redraw to force an exact sample size or desired outcome.
- Each of four cohort/penalty settings selects the lowest-objective prior B120 candidate in each strict joint final direction of 3a versus 3b and 4a versus 4c, among final-1a-first candidates. Its original nine-centroid initialization is used, not its fitted nodes. Each of those four B permutations is paired with initial prefixes [1a,1b,1.5a,2a] and [1b,1a,1.5a,2a]. These are eight declared enriched starts, not an exhaustive search. Final 1a-first is an outcome/readout filter, not a fitting constraint.

## Fixed choices and package implementation

Each cohort keeps its full positive-eigenvalue MDS coordinates fixed. Initial group centroids also remain fixed at their full-data, equal-physical-sheet-weighted values. This is resampling of the point cloud conditional on that geometry and those anchors; it neither repeats distance estimation/MDS nor claims held-out prediction validation.

All fits grow an actual ElPiGraph unbranched curve to 28 nodes with mu=0.1, followed by the previously package-verified fixed-topology refinement to convergence. The fitting objective is nearest-node weighted squared error plus lambda times the sum of squared edge lengths plus mu times squared harmonicity deviations. Readout positions are closest-segment, normalized arc length, which is a separate projection and is not historical time. Chain direction is aligned using the same full-reference early and late group means on every run.

Installed ElPiGraph 0.3.4 samples X[SelPoints] with np.random.random(N)<=ProbPoint but forwards the original PointWeights without subsetting in its wrapper. Therefore the exact native Bernoulli mask law is applied explicitly outside the wrapper; selected original weights are subset and renormalized before unmodified package growth with nReps=1 and ProbPoint=1. The package itself is not patched. Native equal-weight versus explicit-mask equivalence is checked in the two-fit pilot. EmbPointProb remains 1. AvoidResampling controls initial-condition-point reuse and does not disable data subsampling. The native nReps>1 wrapper additionally constructs an average curve; no such extra curve is counted as an ensemble realization here.

Paragraph original weights assign equal total mass to each physical sheet and equal weights within a sheet. Selected weights retain their original relative values and are normalized to unit mass. In point subsampling, remaining sheet masses can differ; in sheet-block subsampling, remaining complete sheets retain equal mass. Pool original weights are equal. This distinction is explicit and is why the block control is included.

## Readout, uncertainty and storage

Each fitted curve projects every original observation, including omitted rows. Raw projected positions and inclusion masks are saved, allowing either all-reference readouts or included-only readouts. Primary position frequencies use only realizations that retained that row; pair frequencies use only realizations retaining both rows, with ties reported separately at 1e-7 normalized arc length. All-reference group medians and retained-only group medians are both saved and labeled. A missing retained group has a null median. The equal weight over starts and lambdas is a sensitivity-design choice, not a posterior or probability of historical order. Random sampling variability, initialization variation, and parameter variation remain separately inspectable.

No ensemble of full-dimensional models is retained. Each realization keeps positions, first-three-coordinate nodes, energies, convergence, and group summaries; masks are stored once per setting. Six 800-fit compact JSON files plus summaries are expected to occupy approximately 40-60 MB, below the study's 100 MB permanent-data allowance. Working dependencies are shared with the H-continuation study and removed only after all verification. Two worker processes use one BLAS thread each. Expected wall time is 20-40 minutes, to be updated using the pilot. Source and input hashes and exact mask seeds are saved. No confidential source is used.
