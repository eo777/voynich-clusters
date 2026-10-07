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

Sparse text and the distance geometry limit interpretation. Clustering and
outlier scores are exploratory rather than tests establishing distinct languages.
The viewer's expandable method notes describe the calculations in more detail.
