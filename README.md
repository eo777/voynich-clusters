# Voynich Clusters

Expanding Currier's voynich languages beyond A / B

A static, self-contained matrix viewer for manuscript pages and panels.
Open `index.html` directly or serve it from any static host. No installation,
backend, CDN, analytics, or uploads are required. Downloads are generated locally.

Golf and RZ keep their fixed page classifications; Ward orders pages within each
group. The overview defaults to manuscript order and can group pages by section,
Currier A/B, Golf, RZ or LFD hand. Its gaps follow the selected groups, and its
bifolium and quire tracks follow the displayed page order. Outliers and anomalies shows a ranked
bar chart of nearest external JS / split-half JS, followed by a raw-distance
scatter for broad row distinctiveness, with no mode switching. Ranked pages and
references meet the 60-word / 20-words-per-half cutoffs. The ratio is a descriptive
G14 measure. Paragraph text is the default and reproduces Golf's saved outlier
results for all 187 supported pages. All text adds circular and radial loci.
The matrix, Ward orders, diagonal, neighbors and charts switch together.
Each chart is clickable and has its own SVG download.
Manuscript and physical bifolium matrix orders also offer section gaps. The default
branch gap is 2 pixels. Branch depth is a display cut,
not an inferred language count. RZ codes represent published marker rules.

## GitHub Pages

Commit these files to your repository. In Settings > Pages, choose Deploy from
a branch, then your branch and /(root). Keep `.nojekyll` beside `index.html`.
The site works at a project subpath such as https://OWNER.github.io/REPOSITORY/.
See https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site.

## Data and interpretation

Paragraph text uses ZL P loci across 207 nonempty pages. All text uses P/C/R loci
across 225 pages. Both use context-weighted Jensen-Shannon divergence in bits.
Labels are excluded in both modes. The displayed diagonal measures alternating
line halves; clustering uses the zero diagonal. Small samples and the
non-Euclidean geometry limit interpretation; see the viewer's method notes.

Source categories are transcription tags, not visual judgments of prose.
fRos includes Pb:37, Ca:3, Cc:338 (378 words; line halves 188/190). Its P-only
halves 20/17 fail the original sample rules; 152 readable L words are excluded.
CSV metadata uses P/C/R source-locus headers and retains full source codes.

RZ method: https://www.voynich.nu/extra/rz_lang.html (page-label snapshot 2026-10-04).
Source-matrix SHA-256: 1a6e1c062b44d98e9ac3ff58239e9b7da898f020a9d2b6cdedb450139148d669

Global Ward view: omitted.
This folder is generated. Edit the source app and rerun its Python publisher.
`site-manifest.json` records the exact generated file hashes. Manually supplied
files such as CNAME and LICENSE are preserved by the exporter.
