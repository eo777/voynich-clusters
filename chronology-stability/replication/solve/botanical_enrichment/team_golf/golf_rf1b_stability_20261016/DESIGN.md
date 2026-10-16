# RF1b-er / basic EVA ElPiGraph stability comparison

Declared 16 October 2026 before fitting. Requested comparison: RF1b-er using individual basic EVA characters, against the completed ZL / Zattera-unit ensemble. This changes both transcription and unitization. It is not an isolated estimate of the transcription effect.

## Matched design

Use exactly the original 207 paragraph pages and 50 physical bifolia, with the same Golf labels, pool membership, equal-sheet weights, 100 Bernoulli(0.9) masks, eight starting group orders per cohort/penalty, lambda 0.01 and 0.05, mu 0.1, 28 nodes, orientation and retained-row readout. The page-point, page-sheet and pooled-sheet settings total 4,800 new fits. Masks and initial group permutations are read from the frozen ZL ensemble. RF centroids are recomputed from RF geometry. No fitted ZL nodes are used as RF initial nodes.

The eight starting permutations were originally selected using the ZL B-permutation search. Holding them fixed makes this a matched experiment conditional on that search coverage; it is not a new exhaustive RF search. Existing ZL fits are reused without modification. No private or unpublished paper is an input.

The public source is https://www.voynich.nu/data/RF1b-er.txt (IVTFF Eva- 2.0, D9). Its already downloaded text is copied into this study and hashed. Use the existing IVTFF parser, all-spaces / preferred-first policy, P loci including opening lines, and readable tokens only. Each plain EVA character is one transition unit, including uncommon plain letters present in the source. RF contains no extended or ligature atoms among readable tokens. Word transitions have true word boundaries. Joined-line transitions remove word spaces, retain true line boundaries and never cross unreadable tokens or separate transcription lines.

As in the original study, build unsmoothed conditional transition probabilities; retain a source row for a pair only when both have at least three outgoing observations; weight its base-2 JS by the mean outgoing-row share, without renormalizing excluded shares. Average word and joined-line JS equally. Pool counts before normalization. Fit separate paragraph and pool classical MDS, treating JS as squared distance and keeping eigenvalues >1e-10. Recompute RF distances and geometry once, then hold them fixed during resampling. Do not recluster the Golf labels.

## Verification

Before fitting, independently reconstruct the frozen ZL word/line transition counts and distances with the same parser/counting implementation in Zattera mode; require equality of counts and distance error <1e-12. Check all 207 RF pages are available and pools contain the same pages. Require finite symmetric distances and no invented/mapped observations. A native-versus-explicit-mask package pilot checks growth equivalence. A matched ZL replay and RF replays verify execution and weighted objective accounting. After fitting, check all masks, weights, energies, convergence, projections and run counts against the original audit, plus cross-representation mask/start equality.

## Charts and interpretation

Render the original retained-fit frequency strips for RF separately. A combined frequency chart gives each representation one-half weight for every observation; preserve multimodality rather than fitting a Gaussian. Also provide a side-by-side view so between-representation differences remain visible. Reproduce the whole-sheet sampling control separately. Starts, penalties and representations have equal design weight, not posterior weights. Shared masks and the two representations are dependent; doubling fits does not double independent samples. Normalized arc positions are relative coordinates in each independently fitted geometry, not common physical distances or dates. Pair/group reversal frequencies and within-mask best-recovered comparisons provide ordering context.

## Storage and execution

Workspace baseline: 18.31 GiB. Plan <250 MiB additional peak disk (temporary pinned runtime, <100 MiB fits/summaries, <20 MiB plots and inputs). Retain compact positions, energies and first-three-coordinate nodes; no ensemble of full-dimensional models. Count/distance NPZ stays <10 MB. Remove the temporary dependency runtime after verification, retaining versions and source hashes. Use six worker processes with one BLAS thread each; the scientific design is unchanged by worker count. Run the workspace storage check after completion.

Method: Albergante et al. (2020), section 2.10 / Figure 1C, https://doi.org/10.3390/e22030296 . ElPiGraph 0.3.4. Original implementation/design: ../golf_published_stability_20261016/ENSEMBLE_DESIGN.md .
