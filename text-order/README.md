# Reproduce scores, or fit your own text with ElPiGraph

The published Path 1 scores are **conditional on a declared design**. They are
not a discovery of historical direction from text alone. `path-v2/reproduce.py`
recalculates the consensus from saved fits. It does not fit curves. There are
now two additional commands below, both using the documented ElPiGraph engine.

| Command | Starting data | What it does |
| --- | --- | --- |
| `python path-v2/reproduce.py --site-root . --out replay` | Published fitted ensembles | Exactly rebuilds the published score tables, without fitting |
| `python text-order/refit_published.py --out fresh-published` | Raw ZL and RF text, explicit published design | Recounts transitions, rebuilds MDS, grows/refines 6,400 fresh ElPiGraph curves, then calculates the four v2 paths |
| `python text-order/fit.py --input my-text.csv --out my-order` | Your text | Fits new ElPiGraph curves without Golf labels, saved curves, or group-order filters; outputs an order and its equally valid reverse |

Run these commands from a clone of the **whole repository**, not a downloaded
script alone. The scripts resolve their dependencies relative to themselves,
so the working directory otherwise does not matter. They never overwrite an
existing output directory. They run locally and do not upload text.

## Install

The published runtime used Python 3.9. Use a virtual environment with Python
3.9 through 3.12 and the pinned dependencies:

```text
python -m venv .venv
```

Activate it (`.venv\Scripts\activate` on Windows;
`source .venv/bin/activate` on macOS/Linux), then:

```text
python -m pip install -r chronology-stability/requirements.txt
python text-order/test_text_order.py
```

NumPy, BLAS, operating system and package versions can change floating-point
rounding and occasionally the recovered local solution. The pinned Python 3.9
runtime is the closest replay environment. Frozen-fit consensus replay is the
way to check the published numbers exactly without those numerical differences.

## Refit the published scores from raw text

```text
python text-order/refit_published.py --out fresh-published
python text-order/verify_refit.py --output fresh-published
```

Inputs are the included `ZL_ivtff_raw.txt`, `RF1b-er.txt`, and
[`published-design.json`](published-design.json). The command reads no saved
transition matrices, MDS coordinates, curve nodes, positions, or fit objectives.
The design contains cohort membership, physical bifolia, Golf labels, the ZL
unit inventory and the original eight anchor permutations for each setting.
It records raw-text hashes and refuses changed transcriptions; use `fit.py`
for a new text. These anchor permutations were selected in the earlier
B-permutation exploration. Here they are explicit **inputs**, not rediscovered
findings. Rebuilding that exploratory choice is outside this command.

Processing is the documented paragraph-text pipeline: canonical IVTFF
`all_spaces` / `preferred_first`, unreadable-word exclusion with line barriers,
ZL/Zattera and RF/basic EVA unitization, word and joined-line transition counts,
row-weighted Jensen-Shannon divergence with at least three outgoing events in
both profiles, equal word/line weighting, and all positive classical-MDS axes
above `1e-10`. Pooled counts are summed before normalization. MDS is computed
once per cohort/representation and held fixed during resampling, as in the
published study.

It uses the original 100 Bernoulli 90%-inclusion masks, equal physical-sheet
weights, eight starts, 28 nodes, lambda 0.01/0.05 and mu 0.1. Growth calls
**`elpigraph.computeElasticPrincipalCurve` from elpigraph-python 0.3.4**.
The existing `core.py` then performs the documented fixed-graph refinement
and closest-segment normalized arc-length projection. Source hashes and
runtime versions are written in `refit-audit.json`.

The 6,400 fits are the eight settings actually consumed by Path v2:
two representations × two penalties × two cohorts × 100 samples × eight
starts. The additional 3,200 whole-sheet page-control fits in the original
stability studies are not needed for these scores.

Assumptions retained for published-score reproduction:

- Starts pass through specified Golf group centers.
- Curve direction is chosen using Early A (1a/1b) and Bio/Stars (4a/4c).
- All v2 paths require the fitted group median of 1a to precede 1b.
- Paths 2–4 additionally require their named 3a/3b and 4a/4c directions.
- Eligible fits are selected by minimum objective **within** each
  representation/penalty/sample. Pairwise fractions are normalized per setting
  before equal-weight averaging. No objectives are compared across penalties.

Output includes fresh ensembles, fit identifiers/objectives, audit information,
and `fresh-published/path-v2/{pages,bifolia}.csv`. Expect roughly one to two hours
on a single CPU process, depending on hardware, and under 100 MB of result data
for the published cohort. No model checkpoints or large arrays are retained.
A short installation check fits one start for one sample in each setting:

```text
python text-order/refit_published.py --replicates 1 --start-limit 1 --out refit-smoke
python text-order/verify_refit.py --output refit-smoke
```

This executes eight real growth/refinement fits. It is explicitly a smoke test
and does not claim to reproduce a 100-sample consensus.
The verifier reads the frozen release only **after fitting**, compares every
requested position and objective component at absolute tolerance `1e-10`, and
checks all four score/rank tables for a full run. It exits unsuccessfully on a
mismatch. See [the recorded verification scope](VERIFICATION.md).

## Run your own text

Supply a UTF-8 CSV with `id,text`. Each row is one original text line and
repeated IDs are pooled into one observation. CSV quoting allows commas in
text; preserve original line boundaries as separate rows. Optional `block` supplies a physical sheet/document for
equal-block weighting; absent blocks default to individual observations.
Optional `group` is used **only after fitting** to summarize the results.
Repeated rows for an ID must have identical metadata. The input row order is
not used as an ordering proposal.

```csv
id,text,block,group
page-a,"the words on the first line",sheet-x,section-one
page-a,"the next line of this page",sheet-x,section-one
page-b,"words on another page",sheet-x,section-two
```

Use at least four observations; this three-row snippet only illustrates the
format. A small runnable artificial fixture is [example.csv](example.csv):

```text
python text-order/fit.py --input text-order/example.csv --samples 3 --starts 2 --nodes 8 --inclusion 1 --out example-order
python text-order/fit.py --input my-text.csv --out my-order
```

The default unit is a Unicode character, preserving case and punctuation.
Whitespace separates words. No language normalization is silently applied.
`?` marks an unreadable token; that whole token is excluded and breaks the
joined-line stream. `^` and `$` are reserved boundary symbols. `--units raw`
uses the EVA atomizer; `--units zattera` uses Zattera's grouped units. Those
two choices are for EVA transcriptions, not arbitrary natural-language text.

This uses the same JS, full positive MDS, actual ElPiGraph growth, fixed-graph
refinement, projection and pairwise-consensus calculation. Its **initialization
design differs** from the published label-based design: eight starts comprise
a weighted PCA line, a PCA-bent line, and six seeded random-direction lines.
These are declared initializations supplied to ElPiGraph, not a claim that
ElPiGraph specifies this particular set. Each starts with three nodes and
grows to 28. `--starts`, `--nodes`, `--samples`, `--seed`, `--penalties` and
`--inclusion` expose the design. The defaults retain the published two penalties,
mu 0.1, 100 samples and 90% inclusion. No labels enter initialization,
orientation, fitting, or selection; no ending order is filtered.

A full-data minimum-objective fit at the **first specified penalty** provides
a reference for aligning arbitrary reversals. Its sign follows PCA1 as a
display convention. Selected resampled fits are reversed when that reduces
squared position disagreement with this reference. Alignment margins and
every start's objective/convergence are retained in `run.json`. In particular,
this procedure cannot determine historical direction or prove that the text
has a unique one-dimensional order. A curve is a model assumption.

`scores.csv` contains score, reverse score, average tied rank and reverse rank.
`pairwise.csv` contains later-plus-half-tie fractions (row relative to column).
`run.json` records masks, all selected positions, candidate objectives,
convergence, raw-input hashes, versions and post-fit group comparisons. It
fails rather than inventing comparisons when a pair is never jointly sampled.
Runs on small datasets may require `--inclusion 1`; those runs have no
subsampling variation. No arbitrary-text score is called the published Path 1.

## Q13 and Q20

Fit the included ZL paragraph text without Golf labels or order filters:

```text
python text-order/fit.py --input chronology-stability/replication/ZL_ivtff_raw.txt --format ivtff --units zattera --metadata text-order/published-design.json --out zl-unlabelled-pages
python text-order/fit.py --input chronology-stability/replication/ZL_ivtff_raw.txt --format ivtff --units zattera --metadata text-order/published-design.json --pool --out zl-unlabelled-bifolia
```

Here the metadata supplies **only page IDs, physical bifolium membership and
quire labels for post-fit reporting**. Golf labels and anchor permutations in
that same file are ignored. Q13 is `M` and Q20 is `T` in the source quire codes;
their pages/bifolia appear in `scores.csv`, and `run.json` contains their median
scores and the `M|T` pairwise comparison. The sign reverses with the equally
valid reversed output. This can test an unoriented relative arrangement; it
cannot establish that Q13 historically preceded Q20. It need not match the
published conditioned route. To analyze only Q13/Q20, export their text lines
to the CSV format and run that subset; changing the comparison population
changes the geometry and peer-relative scores.

## Source method and limits

See [Albergante et al., ElPiGraph (2020), section 2.10](https://doi.org/10.3390/e22030296)
for ensembles of graph fits, and the
[authors' Python implementation](https://github.com/Albluca/ElPiGraph-Python).
The unchanged engine is included at
[`core.py`](../chronology-stability/replication/solve/botanical_enrichment/team_golf/golf_stability_20261016/core.py),
with the original [ensemble design](../chronology-stability/replication/solve/botanical_enrichment/team_golf/golf_published_stability_20261016/ENSEMBLE_DESIGN.md).
Text features, initializations, orientation, selection and consensus scoring
are this project's documented adaptations. Neither ElPiGraph nor a good curve
fit makes these scores dates, chronology probabilities, or evidence of causal
direction. The original [Path v2 definitions](../path-v2/README.md) remain the
contract for the published scores.
