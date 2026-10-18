# Sequence / path collection v2

**Default: `v2_default` — overall pairwise consensus.** Golf **2.1** identifies
the classification scheme; path collection **v2** identifies the ordering data.
The classification memberships have not changed.

## Named paths

| ID | Meaning |
| --- | --- |
| `v2_default` | Overall consensus; no ending-order filter |
| `v2_3a_3b__4a_4c` | 3a before 3b; Bio 4a before Stars-S 4c |
| `v2_3a_3b__4c_4a` | 3a before 3b; Stars-S 4c before Bio 4a |
| `v2_3b_3a__4a_4c` | 3b before 3a; Bio 4a before Stars-S 4c |
| `legacy_proposed` | Original proposed curve, available in the selector and archive only |

The default's group-median scores put **3b before 3a and 4c before 4a** in both
cohorts. The fourth conditional combination is omitted because it has these
same directions. The default still combines every eligible ending direction;
it is **not** filtered to that fourth combination. These pair constraints do
not specify the order of every B group or prevent interleaving of their pages.

## Calculation

The inputs are the frozen ElPiGraph ensembles from the published stability
study: 100 90%-inclusion samples, eight initializations, lambda 0.01 and 0.05,
and mu 0.1. Pages use individual-page sampling; bifolia use pooled-bifolium
sampling. All sections, including unassigned observations, remain in the peer
population. Both ZL/Zattera and RF1b-er/individual basic EVA are included.
This changes both transcription and glyph grouping, not transcription alone.

Within each representation, penalty and sample, first retain converged fits
whose full-reference group medians place 1a before 1b. For a conditional
variant, also require its two strict ending directions. Ties fail a strict
condition. Then select the minimum-objective eligible fit, resolving exact
objective ties by start index. The objective is compared **only within one
representation/penalty/sample**. An unavailable sample is omitted and counted;
another route is never substituted to fill it.

For each pair of observations, count the selected fits in which both were
retained. Compute the fraction in which the focal observation is later, plus
half the tied fraction (tolerance 1e-7). Normalize these fractions separately
for each representation/penalty setting, then average the four settings with
weight 1/4 each. The observation's score is the mean of these values over all
other members of its cohort. Lower scores mean earlier relative positions.
The score is not stretched to force its endpoints to 0 and 1.

This exactly reproduces the individual inspector's presets: both equal
weight, both penalties equal weight, best per sample, pairwise consensus,
all sections. The conditional variants filter **before** selecting the best
fit, so they describe the best available fit within the named scenario, not
only the subset of unconstrained winners with that scenario.

The datasets contain 207 page scores and 50 pooled-bifolium scores per path.
They are separately calculated peer-relative scores, not positions along one
shared arc-length axis. Their difference is descriptive; it is not a geometric
page-to-curve residual. There is no v2 score for a page without a paragraph
observation. Such fields are blank, not zero, and are not filled from the old
circular/radial projections. In the assignment CSV, a page's paragraph score
is repeated across that page's text-kind rows; its pool score is repeated
across the constituent pages. These do not constitute label-specific estimates.

## Files and support

- [pages.csv](pages.csv) and [bifolia.csv](bifolia.csv): IDs, Golf 2.1 labels,
  scores, ranks and retained-fit counts for the four v2 paths.
- [collection.json](collection.json): full definitions, source hashes,
  per-setting sample support, group-median scores and pairwise summaries.
- [selected-fits.csv](selected-fits.csv): exact sample/start identifiers and
  objectives used in each conditional or default ensemble.
- [reproduce.py](reproduce.py): regenerate these files from the published
  frozen ensembles with NumPy; no ElPiGraph fitting or downloads are required.

Run from the site checkout:

```text
python path-v2/reproduce.py --site-root . --out reproduced-v2
```

Each setting has 100 eligible samples for the default. Conditional scenarios
have unequal support; consult `eligible_samples`, `omitted_samples`, and
`minimum_joint_samples`. Even the sparsest published scenario has a jointly
retained comparison for every pair in every setting. Equal representation and
penalty weights are maintained despite unequal conditional sample counts.
Ranks resolve exactly equal scores by the original observation order; they
should not be read as evidence distinguishing tied observations.

## Curves and chart interpretation

Sequence and Curves keep their existing MDS coordinates and chart layouts.
The selector changes the scores, section chart, bifolium table, offsets and
within-bifolium score-variation statistic together. A consensus is not itself
one fitted curve. V2 draws smooth **display guides** through the existing
eight-dimensional ZL coordinates as functions of the fixed consensus scores.
The guides use Gaussian local linear regression, bandwidth 0.06, evaluated at
160 positions. Separate guides use pools, pages, or both. They do not refit
ElPiGraph, project observations, or alter the consensus ordering. RF/ZL
coordinates are never averaged or overlaid as though they shared a basis.

The held-out statistic and score confidence interval from the original fit
do not apply to v2. They are shown as unavailable in v2, and their original
values return when the legacy path is selected. V2 within-bifolium variation
is recomputed from the v2 page scores. Historical fit overlays remain explicitly
labeled as legacy comparisons. The old data are preserved under
`archive/path-v1/`; the current assignment CSV contains only v2 path columns.

These are sensitivity-derived consensus scores, not chronology probabilities,
dates, calibrated confidence intervals, or proof of direct connections.
Pairwise preferences can contain cycles. A consensus ranking need not equal
any single fitted route, and filtering group medians does not force every
individual page to obey a block order.

The ensemble fitting is motivated by [Albergante et al., ElPiGraph (2020)](https://doi.org/10.3390/e22030296).
The conditional selection and pairwise consensus score above are explicitly
documented downstream summaries, not a claimed ElPiGraph-native probability.
See the [stability report](../chronology-stability/index.html) and
[individual inspector](../chronology-stability/individuals.html) for the source
methods and distributions.
