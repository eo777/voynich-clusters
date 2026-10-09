# Voynich Clusters

Expanding Currier's voynich languages beyond A / B

Voynich Clusters explores similarities between manuscript pages using
Jensen–Shannon divergence between within-word transition matrices.

The Golf regimes come from npcompl33t's
[analysis of languages beyond Currier A/B](https://www.voynich.ninja/thread-6084.html).
Golf is the project's internal name for those classifications. The viewer also
includes [René Zandbergen's language classifications](https://www.voynich.nu/extra/rz_lang.html).

Golf 2.0 adds the revised profile assignments from the 8 October 2026 membership
audit. Both Golf versions appear in the overview and are available as matrix
groupings. Diagonal stripes identify six intermediate candidates. They retain
their exported working membership while their candidate boundaries remain unconfirmed.

The **Golf V2** tab reproduces 11 saved faceted transition figures and eight
internal label comparisons as interactive JavaScript plots. Select a point to
inspect its exact counts, denominators and assignment. Page circles use the same
support fills as the overview. Diamonds identify pooled bifolia or within-page
blocks, which do not have page-level support ratings. Grey crosses are unassigned
or outside the comparison. Coordinates and rates remain unchanged. The accompanying tables
show conditional transition rates, TF-IDF word examples and LSA associations.
Append `#golf-v2` to the site URL to link directly to this tab. Switching tabs
updates the URL, and browser Back and Forward restore the selected tab.
The branch view defaults to 225 pages without labels, using each page's largest
P/C/R text kind, with prose winning ties. This is a display rule, not a joint
whole-page classification. Labels and all text kinds are separate options,
covering up to 17 profiles and 227 pages. Dashed relationships remain provisional,
and the tree's spacing does not represent measured distance.

**No audit recommendation** means the audit withheld a recommendation. It does
not mean no working classification or affinity exists. Page details retain the
previous preliminary class, nearest-profile affinity, alternatives, and the
reason for withholding a recommendation. Grouping uses the explicit membership
export's version 2 working assignments. The export retains a previous class with
uncertainty when the audit abstains; nearest affinity never determines membership.
Full-color fills show strong support. Diagonal stripes in the group color and a lighter shade show tentative support,
and lighter fills with white diagonal stripes show unverified membership.
The exported `assignment_support` field
determines these levels. Strong support includes validated bifolium assistance;
unverified does not mean probably wrong. Key counts follow the selected view.
Diagonal stripes in the two candidate colors show intermediate cases without an
additional support pattern. Their support remains in details, exports and key
counts. Grey cells
have no established working class. Each unit appears only in its working group;
candidate alternatives are shown in its details. The strict and
supported analysis masks remain separate.
Repeatability percentages identify their target affinity and are not assignment
confidence scores.

Page details include Golf's authored case discussions where available, with
source attribution. Units without a separate writeup receive a clearly labeled
explanation of their saved audit measurements. The supporting table separates
page evidence, assistance from its bifolium, and the whole-bifolium comparison.

Current identifiers follow the Golf codebook, including 3c Stars-T, 4c Stars-S,
3d-c Rosette circles, 3d-l Rosette labels and 6 f58. V2 uses distinct hues rather
than lightness variations of the original Golf colors. Herbal B uses pink and
brown; Stars uses violet and gold. Lighter fills indicate qualified support.

## Explore the manuscript

- Compare pages in the matrix, select individual pairs, and zoom into patterns.
- Regroup the overview by section, Currier language, Golf regime, RZ language,
  LFD hand, quire, or bifolium, with breaks between the selected groups.
- Inspect outlier scores and external distances, and download charts or data.

## Page classifications and scores

[Download the page classifications CSV](golf-classifications.csv) for all 225
pages and panels in the viewer, in manuscript order. Each row gives the Golf
classification and bifolium, followed by separate `paragraph_` and `all_` scores.
Unassigned pages are retained. These are the viewer's existing classifications.
The `golf_v2_assignment_code` and `golf_v2_assignment_id` columns give authoritative
working membership. `golf_v2_assignment_source` and `golf_v2_assignment_status`
record its provenance and qualification. These match Golf's `page_membership.csv`
for the controlling text kind on each page.
`golf_v2_assignment_support` gives Strong support, Tentative support, Unverified
membership or Unassigned. Detailed assignment and audit statuses remain separate.
The separate `golf_v2_class` and `golf_v2_code` columns retain the audit ID and code,
text kind, assessment status, intermediate-candidate flag and competing profile
IDs. A blank Golf V2 class means no audit recommendation. Previous working class,
nearest affinity, reference membership, recommendation provenance and analysis
masks have separate columns. `golf_v2_display_class` and `golf_v2_placement_source`
record working membership placement separately from the audit recommendation.
Membership labels, support tiers, evidence basis and warnings carry the
qualifications. Affinity repeatability has its own target-profile column and
must not be read as confidence in the retained assignment.

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

Golf 2.0 classifies text kinds separately. Paragraph mode uses the prose
assignment. All-text mode groups a page by its largest included text kind, with
prose winning ties; selected-page metadata retains every text-kind result.
Label text never controls a page's matrix placement. All-text distances still use the full
selected text, so the grouping is not a new classification of combined text.

Assignments remain fixed within each text kind. Ward clustering orders pages within
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
