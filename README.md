# Voynich Clusters

Expanding Currier's voynich languages beyond A / B

Voynich Clusters explores similarities between manuscript pages using
Jensen–Shannon divergence between within-word transition matrices.

The Golf regimes come from npcompl33t's
[analysis of languages beyond Currier A/B](https://www.voynich.ninja/thread-6084.html).
Golf is the project's internal name for those classifications. The viewer also
includes [René Zandbergen's language classifications](https://www.voynich.nu/extra/rz_lang.html).

## Explore the manuscript

- Compare pages in the matrix, select individual pairs, and zoom into patterns.
- Regroup the overview by section, Currier language, Golf regime, RZ language,
  or LFD hand while following the bifolium and quire tracks.
- Inspect outlier scores and external distances, and download charts or data.

## Page classifications and scores

[Download the page classifications CSV](golf-classifications.csv) for all 225
pages and panels in the viewer, in manuscript order. Each row gives the Golf
classification and bifolium, followed by separate `paragraph_` and `all_` scores.
Unassigned pages are retained. These are the viewer's existing classifications.

Nearest and second-nearest external JS values include small samples and exclude
the page's own bifolium. The two neighbors may share a bifolium with each other.
JS values are in bits. The outlier score uses the nearest external page that
meets the sample rules, divided by internal split-half JS. Its reference page
and distance are included because it can differ from the unrestricted nearest
neighbor. Blank outlier scores mean the page does not meet the sample rules or
has zero internal JS.

The isolation percentile range compares equally short samples from Early A,
using random words and contiguous windows. Higher values mean greater isolation
for the word count. Comparisons with fewer than five source bifolia are marked
as limited. Blank percentiles mean the page was not studied or had no qualifying
sources. Entirely blank paragraph fields mean the page has no paragraph profile.
For f57v, `all_no_singletons` identifies the 68-word isolation comparison that
excludes 124 single-unit words and filters comparison pages too. Its other
all-text metrics still use all 192 words.

## Text and method

Paragraph text is the default and covers 207 pages and panels. All text adds
circular and radial text, covering 225. Both use the ZL transcription and exclude
labels. Text categories follow the transcription's locus codes.

Golf and RZ classifications remain fixed. Ward clustering orders pages within
each group. The diagonal compares alternating line halves within a page.
Outlier scores compare the nearest page on another bifolium against this
within-page variation. Ranked pages require at least 60 words and 20 per half.

For studied pages, a separate ranking compares isolation with equally short
samples from longer pages. Choose a reference regime to inspect the percentile
range from random words and contiguous windows. This measure requires neighbors
on two distinct external bifolia. Comparisons with fewer than five source
bifolia are displayed separately without a rank. Unstudied pages have no score.
For f57v, this ranking uses 68 words after excluding single-unit words from the
page and its comparison material. This exception is labeled in the viewer;
the matrix and scatterplot retain the full selected text.

Sparse text and the distance geometry limit interpretation. Clustering and
outlier scores are exploratory rather than tests establishing distinct languages.
The viewer's expandable method notes describe the calculations in more detail.
