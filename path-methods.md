# The curved manuscript sequence

The Curves and Sequence tabs explore the manuscript's transition-profile gradient
and a candidate ordering along it. They use a shared MDS containing 207 paragraph
page profiles, 16 non-paragraph profiles and 50 pooled bifolia. Hiding a point does
not remove it from the embedding or refit the curve.

## Representations and fits

Word transitions preserve word boundaries. Joined-line transitions remove spaces
within each source line and use line boundaries. Combined divergence gives the two
Jensen–Shannon divergences equal weight. Classical MDS is applied to the square
root of divergence. Axis percentages refer to positive eigenvalue mass.

The Curves tab compares three representations on the same 273 profiles. Each
curve is fitted to a different subset of that shared geometry. The four subsets
are 50 pools, 207 paragraph pages, their union, and all 273 profiles. These are
full-space fits with the same 11-knot spline specification, not curves fitted to
the two displayed axes. The supplementary figures also compare embeddings made
separately from pages alone and bifolia alone.

The Sequence tab uses the proposed Early A return in the combined representation.
Its initialization allows this return before freely fitting in all 254 positive
MDS dimensions, using 11 knots and smoothing penalty 0.05. The default route fits
almost as well. The proposed turn is also available in the combined Curves view.
It was not fitted in the word-only or joined-line comparison displayed there.

Positions are full-space projections onto the fitted path, normalized to the
range of the 50 pool projections. Zero begins at the group 1 end and one reaches
the Bio side. A page can project beyond either end of that pool range. These are
relative positions, not dates, elapsed time or an outlier score. Distance from
the curve is retained separately in the selected-point details and downloads.

## Bifolium agreement

For each coordinate, subtract the mean of a bifolium's paragraph pages from each
of its pages and sum the squared deviations. Divide that within-bifolium sum by
the total squared deviation of all 207 paragraph pages from their overall mean.
Every page has equal weight. The resulting percentage measures how much page
variation remains within bifolia, with the rest lying between them.

MDS1 has 2.287% within-bifolium variation. The proposed path has 0.758%, a 66.9%
reduction. MDS2 has 15.07% and MDS3 has 22.66%. The chart includes all eight saved
axes on a common linear scale. The CSV retains their unrounded values.

The held-out calculation omits each tested bifolium from curve fitting. The
shared MDS stays fixed, including the tested profiles. Each fold's squared
within-bifolium deviations are normalized by the variance of all page projections
onto that fold's curve, then summed over folds and divided by 207. Its result is
1.549%. This tests curve-fitting dependence, not a fully independent refit of the
embedding. It is not an ANOVA of stitched positions from different folds.

The per-bifolium table uses a different descriptive measure, the SD of its page
positions divided by overall page SD. The offset chart shows each page minus its
pooled profile, in overall page SD units. The main percentage instead uses the
mean of a bifolium's pages, so it does not reward closeness to a pool containing
the same counts.

## Sections and classifications

Manuscript section labels come from the cluster map's manuscript metadata.
Herbal A and Herbal B additionally use the Currier assignment. They are not
inferred from the curve. The section chart contains every paragraph page and
shows a mixed pool in each section containing one of its pages.

Golf colors retain the app's palette. A pool takes its dominant constituent
membership by paragraph word count. Its individual assignments remain visible
when selected. Diagonal fills indicate tentative or unverified support on page
profiles. These classification labels supply colors, not inputs to MDS. The
proposed path's initialization did use the Early A grouping to specify the turn.

## Interpretation and downloads

Endpoint identity is sensitive to the fitted route. D:25-32 remains first across
twenty smoothing settings with the current initialization, but falls to eighth
when its own pool is omitted from curve fitting. An alternative initialization
captures almost the same pooled variation and starts at C:19-22 instead. The
Bio-side terminal bifolium also changes with representation and smoothing.
The broad ends are more stable than the identity of the first or last bifolium.
These checks hold the upstream MDS fixed and do not provide date uncertainty.

Tight bifolium agreement is evidence of structure linked to the physical
manuscript. It does not by itself distinguish chronology from other sources of
production-linked variation. The sequence is a candidate for evaluating a
proposed chronology, not an independently established historical order.

- [Bifolia in sequence](sequence-bifolia.csv)
- [Every page and profile](sequence-profiles.csv), including full-space residuals
- [Held-out projections](sequence-heldout.csv)
- [Variation across all eight MDS axes and the path](sequence-variation.csv)
- [Representation comparison](path-representations.csv)
- [Complete chart data and source hashes](path-data.json)

The source studies were saved on 10 October 2026. The exported views reuse their
coordinates and full-space projections. They do not estimate new dates or change
the Golf assignments.
