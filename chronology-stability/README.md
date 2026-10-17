# Golf chronology stability: report, data and replication

[Full writeup](index.html) · [Figures and downloads](figures.html) · [Forum methodology](FORUM_METHODS.md)

This directory contains the completed 16 October 2026 analysis, including the four original PNG/SVG figure pairs, the scEGRET calibration figure, the full reports, final metrics, original analysis scripts, numerical audits, scEGRET's pinned R source and license, and the frozen inputs and prior fitted candidates needed by those scripts. No private or unpublished paper is included or used.

The follow-ups are included in full:

- [Individual bifolia and pages](individuals.html): directly linkable distributions, same-run pair comparisons, ties and a descriptive consensus ordering. Defaults are pooled bifolia, both representations equally weighted, all sections, pairwise consensus, and the best recovered candidate per sample with 1a before 1b.
- [RF1b-er / basic EVA comparison](rf1b.html): 4,800 additional matched ElPiGraph fits, separate RF charts, equal-weight combined charts, side-by-side comparisons and whole-sheet controls. Both transcription and glyph grouping differ from ZL/Zattera; these are sensitivity frequencies, not historical probabilities.
- [scEGRET calibration](replication/solve/botanical_enrichment/team_golf/golf_published_stability_20261016/calibration/report.html): 147,888 production simulated datasets, continuous references, missing-state controls, false-alarm and power checks, and calibration of the supplementary two-sided diagnostic.
- [Actual f58 and astro/stars removals](replication/solve/botanical_enrichment/team_golf/golf_published_stability_20261016/calibration/omission_test/report.html): seven cases applied to the observed pages and 7,596 regenerated continuous controls, with exact omitted IDs, bin memberships, normalized and raw-count results, and both noise models. Each control is reused across the seven masks; these are paired comparisons, not seven independent sets of controls.

The `replication/` directory preserves the original relative workspace layout so the scientific source files can remain byte-identical. Its study directory is:

`replication/solve/botanical_enrichment/team_golf/golf_published_stability_20261016/`

The surrounding folders contain only selected dependencies: transition counts, pairwise distances, metadata, original Golf positions, prior starting candidates and shared mathematical routines. The calculation starts from these frozen observed inputs and saved starting candidates; regenerating the entire preceding classification or B-permutation search is outside this package. The full site also provides the separate `research-data/` export.

## Check the package

From this directory, run `python verify.py`. This uses only the standard library, checks all retained file sizes/hashes, and validates local HTML links. No model is fitted and no files are changed. The package manifest records the original source path and checksum of every scientific file.

Create a Python virtual environment and install `requirements.txt` for scientific checks and figures. These pins match the original Python 3.9 environment; platform-dependent floating-point or plotting differences may occur elsewhere. Run `python check_inputs.py` to independently reconstruct the cohort distances from the included raw transition counts and check the exported scEGRET feature matrices. This is read-only.

## Reproduce the calculations

Copy `replication/` to a separate working directory before rerunning: the original scripts write results beside themselves. Do not overwrite the frozen reference release. In the copied study directory:

1. Install ElPiGraph's exact source where the original code records its hashes: `python -m pip install --no-deps --target _deps elpigraph-python==0.3.4`. The other dependencies are provided by `requirements.txt` in your virtual environment.
2. Run `python ensemble.py --workers 2`, then `ensemble_summary.py`, `ensemble_competition.py`, `ensemble_findings.py`, and `ensemble_audit.py`. The full 4,800-fit calculation took approximately 46 minutes on the original machine. `ensemble_independent_audit.py` adds two verification growth fits. It is not another uncertainty sample.
3. Run `python cassidy_h_continuation.py`, then `cassidy_summary.py`, `cassidy_independent_audit.py`, and `cassidy_return_check.py`. This executes the H-coordinate derivative prediction, assignment-event correction and path checks.
4. Run `python prepare_scegret.py` to regenerate the feature inputs and proposed orders. Install Node.js (original version 24.15.0), then run `npm install --prefix _r_runtime webr@0.6.0`. Run `node --import ./webr_windows_register.mjs ./run_scegret.mjs`. The adapter also works on non-Windows paths. This executes all 648 cases in R; WebR downloads its R packages, while the manuscript data are processed locally. Exact original R package versions are recorded in `scegret_results.json`; future package repositories may supply different versions.
5. Run `python scegret_independent_audit.py` and `python scegret_summary.py`. The original run used parallel partitions; the merge/retention scripts are retained for provenance but are unnecessary for a new single-process run. Deleted intermediate partitions are not included.
6. For the calibration, read `calibration/DESIGN.md`, `MATCHED_DESIGN.md` and `BYPASS_DESIGN.md`. Run `python calibration/run_calibration.py`, `python calibration/matched_power.py`, `python calibration/bypass_calibration.py`, `python calibration/summarize.py`, then `python calibration/build_report.py`. The latter two design files identify additions made after preliminary results were inspected. `python calibration/replay_reference.py` separately verifies two complete reference streams. The fast NumPy score is checked against the original R outputs; no R runtime or curve refitting is needed for these follow-ups. Synthetic matrices are streamed rather than retained.
7. For the actual-page omissions, run `python calibration/omission_test/run_omissions.py`, then `python calibration/omission_test/build_report.py`. Its `DESIGN.md` specifies the seven masks, the fixed original Golf ordering, native rebinning, and targeted versus whole-path comparisons. The full-data reference means and standard deviations must reproduce the prior calibration exactly.

Historical JSON provenance keys contain Windows path separators. The portable package verifier normalizes these when checking frozen files. For the original provenance-checking scripts on Unix, first regenerate the relevant result/input file in the copied workspace so its path keys use the host's separators. The scientific source files themselves are unchanged.

## Rebuild the figures

For the RF comparison, use the sibling `replication/solve/botanical_enrichment/team_golf/golf_rf1b_stability_20261016/` directory and follow its [reproduction instructions](replication/solve/botanical_enrichment/team_golf/golf_rf1b_stability_20261016/REPRODUCE.md). Both studies, the raw public transcriptions and the required preprocessing modules are included. Run `plot_comparison.py` and `build_report.py` there to rebuild all six figure pairs and the report without new fits. The storage command in the original workspace instructions is an internal housekeeping step and is not required in this portable export.

From a copied study directory, run:

```text
python plot_ordering_ensemble.py
python plot_ordering_ensemble.py sheet_control
python plot_cassidy_h.py
python plot_scegret.py
python build_gallery.py
python finalize_manifest.py
```

Figures can be rebuilt from the retained results without another fit or an R installation. `python calibration/build_report.py` also rebuilds the calibration figure and report from retained summaries; `python calibration/omission_test/build_report.py` rebuilds the removal report. Use the matching runtime versions for the closest numerical/visual reproduction. The parent website's `#stability` route opens the full writeup; `chronology-stability/index.html` is its standalone entry point.

## Interpretation

Equal weighting of starts and penalties is a declared exploration design, not a posterior distribution. The report separates all-start position frequencies from best-recovered-score comparisons. MDS coordinates stay fixed during resampling; H maps change only the displayed coordinates of the observed pool. scEGRET was adapted from expression features to text-transition features, without transferring its biological threshold. Its bins can mix or exclude the very groups whose relative order is being compared.

The calibrated results do not robustly identify f58 as a discontinuity, and low power under empirical sheet noise prevents a negative screen from establishing continuity. Removing f58 produces positive flags in some raw-count and alternative-noise-model comparisons, but none with per-million frequencies under sheet resampling. Removing either f58 or its three bin companions makes the original focal bin fail scEGRET's strict minimum occupancy rule, so this experiment cannot isolate f58 as an indispensable bridge. The full reports preserve both positive and negative results. The manuscript-specific reference generators, tail calibration and two-sided diagnostic are our adaptations, not additional methods attributed to the scEGRET paper.

See the full report for the source-paper mapping, numerical results and limitations. Source-specific licenses remain with their code; inclusion does not imply a new blanket license for third-party material.
