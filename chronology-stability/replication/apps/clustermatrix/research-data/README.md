# Transition matrices and MDS data

These CSVs contain the inputs and coordinates behind Voynich Clusters. They
preserve the published analysis as of 12 October 2026. No classifications or
curve fits were changed for this export.

## Start here

- [Profiles](profiles.csv) links every matrix identifier to its page, text kind,
  bifolium and classification.
- [Word transition counts](transitions-word-counts.csv) and
  [conditional probabilities](transitions-word-probabilities.csv) include every
  cell of every profile's 29-by-29 matrix, including zero counts.
- [Joined-line transition counts](transitions-line-counts.csv) and
  [conditional probabilities](transitions-line-probabilities.csv) give the
  second representation used by Golf 2.1.
- [Combined pairwise divergence](pairwise/shared-combined-js.csv) contains the
  273-by-273 matrix used for the shared page-and-bifolium geometry.
- [Combined shared MDS coordinates](mds/shared-combined.csv) includes all 254
  positive dimensions, with the first eight oriented as in Curves and Sequence.
- [Reproduction script](reproduce.py) verifies the matrices and MDS using only
  these CSVs and NumPy. Run `python reproduce.py` from this directory.

## Profiles and text selection

The 592 rows of `profiles.csv` comprise 302 individual page/text-kind units,
225 derived whole-page profiles excluding labels, 15 combined circular/radial
profiles used in the MDS, and 50 paragraph-pooled bifolia. These are overlapping
representations, not 592 independent manuscript pages.

The 302 individual units are 207 paragraph, 27 circular, 11 radial and 57 label
profiles. Their `page|text_kind` identifiers join directly to `unit_id` in the
[Golf 2.1 assignments](../golf-v21-assignments.csv). Derived profiles list their
constituents explicitly. A whole-page profile does not receive a new joint
classification. Pool classifications are descriptive labels from the saved
viewer and their support is recorded as constituent assignments.

`shared_mds_index` identifies the 273 profiles jointly embedded for Curves and
Sequence. These are 207 paragraph profiles, 15 circular/radial profiles, the
Rosette label profile and 50 pools. Supplementary text types are jointly fitted
in this shared geometry. The proposed curve itself was fitted to the 50 pools.
The separate comparison view described below projects non-prose profiles into
a paragraph-only fit instead.

`matrix_paragraph_index` and `matrix_all_index` identify the 207-page and
225-page cluster-map cohorts. Indices start at zero. Empty fields mean absent
or inapplicable, not zero. `shared_mds_id` preserves the viewer's original
identifier, such as `f57v*` or `N:85-86`.

[Pool members](pool-members.csv) lists the paragraph units summed into each
bifolium. Counts are summed before normalization. A pool is not the mean of
its pages' conditional probabilities. Non-paragraph text is excluded from
these pools, even where the same bifolium also has circular or label profiles.

## Transition data

The text source is the ZL transcription, using the existing canonical
`all_spaces` and `preferred_first` parsing policy. Unreadable tokens are
excluded. The matrices use Zattera glyph units, not individual EVA characters.
For example, `ch`, `sh`, `ee` and `ii` are single units. Rare units collapse to
`other`. [Alphabet](alphabet.csv) gives the full ordered inventory.

Each transition file has one row per profile and source unit. The `to_` columns
give all possible destination units. `row_total` is the number of outgoing
transitions from that source. `row_weight` is this total divided by the profile's
total transition count. `row_eligible_ge3` records whether the source has at
least three outgoing observations.

Probability cells are conditional next-unit frequencies, calculated by dividing
the count by `row_total`. An unobserved source has blank probability cells.
Its count cells and row weight are zero. All observed rows sum to one. No
pseudocounts or smoothing are applied.

In the word representation, `^` and `$` mark the beginning and end of each word.
In the line representation, spaces within each transcription line are removed
and the concatenated string is unitized again. These markers then represent
true line boundaries. Unreadable words break a line into separate readable
segments. There are no transitions across those gaps, no artificial boundary
markers at a gap, and no joining across transcription lines. Paragraph-initial
lines are included. Single-unit words are included. The f57v singleton-removal
exception in the separate isolation ranking does not alter these matrices.

## Pairwise divergence

All files in `pairwise/` are square matrices. The first column gives the row
identifier and the remaining headers give column identifiers in the same order.

| File | Comparison |
| --- | --- |
| [shared-word-js.csv](pairwise/shared-word-js.csv) | Word representation for all 273 shared profiles |
| [shared-line-js.csv](pairwise/shared-line-js.csv) | Joined-line representation for the same profiles |
| [shared-combined-js.csv](pairwise/shared-combined-js.csv) | Arithmetic mean of those two divergence matrices |
| [matrix-paragraph-word-js.csv](pairwise/matrix-paragraph-word-js.csv) | The cluster map's 207 paragraph pages |
| [matrix-all-nonlabel-word-js.csv](pairwise/matrix-all-nonlabel-word-js.csv) | The cluster map's 225 pages with paragraph, circular and radial text |

For each source unit, the two conditional next-unit distributions are compared
with base-2 Jensen–Shannon divergence. Their weight is the mean of the two
profiles' `row_weight` values. A row contributes only when both profiles have
at least three outgoing observations. The weighted contributions are summed.
Excluded weights are **not renormalized**. The accompanying
[word eligible weights](pairwise/shared-word-eligible-weight.csv) and
[line eligible weights](pairwise/shared-line-eligible-weight.csv) record the
total contributing weight for each pair. Low weight means little comparable
transition evidence.

The combined score is `0.5 * word_JS + 0.5 * line_JS`. The transition matrices
themselves remain separate. Scores are in bits. These files contain divergence,
not its square root. The diagonal is zero. The cluster-map display's split-half
internal divergence is a separate statistic and does not replace this diagonal.
The page matrices reproduce the viewer to its saved rounding precision.

## MDS coordinates

Classical MDS double-centers the divergence matrix, equivalently using the
square root of divergence as its input distance. It eigendecomposes
`B = -0.5 * J * D * J`, where `J` centers the observations. A positive coordinate
axis is its eigenvector multiplied by the square root of its eigenvalue.
All eigenvalues are retained in the spectrum export. Coordinate files keep
eigenvalues greater than `1e-12`.

| File | Jointly fitted profiles | Positive dimensions |
| --- | --- | --- |
| [shared-word.csv](mds/shared-word.csv) | 273 pages/text profiles and pools | 249 |
| [shared-line.csv](mds/shared-line.csv) | The same 273 profiles | 251 |
| [shared-combined.csv](mds/shared-combined.csv) | The same 273 profiles | 254 |
| [viewer-pages-combined.csv](mds/viewer-pages-combined.csv) | 223 page/text profiles, without pools | 210 |
| [viewer-bifolia-combined.csv](mds/viewer-bifolia-combined.csv) | 50 paragraph pools, without page points | 49 |

[Eigenvalues](mds/eigenvalues.csv) includes positive, near-zero and negative
eigenvalues for all five fits. [Models](mds/models.csv) records their dimensions
and negative eigenvalue mass. The shared combined fit has about 1.15% negative
eigenvalue mass. Retaining only its positive coordinates therefore approximates
the supplied dissimilarities even when all 254 axes are used.

Axis signs match the saved figures for the first eight shared axes and first
four separate-view axes. Remaining signs follow the eigensolver and have no
intrinsic meaning. Axes from separately fitted models are not interchangeable.
Pools and their constituent pages overlap in evidence when fitted together.

[Viewer branches](mds/viewer-branches.csv) exports the four saved coordinates
for every branch, for page, bifolium and joint views. These are independent
branch-specific fits. A branch pool contains only that branch's paragraph
members. `omitted_members` identifies other pages on the physical bifolium.
The `golf21` rows are the reviewed Golf 2.1 branch study. The `golf20` rows retain
the earlier comparison study, including its experimental three-way Early A
labels. Those experimental labels are not official Golf 2.0 assignments.
[Branch statistics](mds/viewer-branch-statistics.csv) gives the saved axis shares.

[Comparison coordinates](mds/viewer-comparison-saved.csv) preserves the earlier
whole-manuscript and Early A close-up charts and their experimental two/three
group labels. `fit_role` distinguishes fitted observations from the 16
non-prose projections into the 207-paragraph context fit.
[Original Golf 2.0 coordinates](mds/original-golf20-saved.csv) preserves the
original 163-unit word-transition map at its saved four-decimal precision.
These saved comparison files contain four axes rather than newly fitted axes.

[Path nodes](mds/path-nodes.csv) preserves the eight-dimensional display
projections of every curve currently available in Curves and Sequence.
`node_index` is only the vertex order, not the normalized sequence coordinate.
[Path projections](mds/path-projections.csv) gives the saved full-space path
coordinate of each profile. The `combined` / `proposed` rows reproduce the
published Sequence path. Other rows identify their representation and fit.
Residual distance is provided where the published payload stores it. A blank
does not mean zero. The curves were fitted in the full positive MDS space.
Their eight display coordinates are insufficient to repeat the full-space fit.
Path positions are not dates or measures of elapsed time.

## Provenance and verification

[Manifest](manifest.json) records input hashes, CSV hashes, row counts and
column counts. Source paths there are relative provenance identifiers from
the research workspace. They are not additional files required to read these
exports. The reproduction script uses only the published CSVs and NumPy.

The export checks every original MDS transition profile against the canonical
corpus loader, all 302 readable-word counts against the assignment export,
both cluster-map matrices within `5.1e-9` bits, and the saved MDS coordinates
within `1e-9`. Probability and count tables preserve full numeric precision.
All CSVs use UTF-8, comma delimiters and one header row.
