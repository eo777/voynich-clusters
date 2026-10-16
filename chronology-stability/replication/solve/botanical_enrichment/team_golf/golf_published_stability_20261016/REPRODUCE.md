# Reproducing this study

All paths below are relative to this study directory. Keep the surrounding repository because the input files and earlier candidate curves are referenced by their frozen paths and hashes. A rerun writes result files at these paths; preserve a completed release before rerunning.

## Python calculations

The exact base interpreter/library versions are in `base_runtime_versions.json`; optional research dependencies are in `python_dependency_versions.json`. ElPiGraph is version 0.3.4. Restore the latter packages into `_deps` if the temporary runtime has been removed. The fitting code loads that directory explicitly. Use one BLAS thread per worker, as the source does.

- `ensemble.py --workers 2` runs the declared 4,800 fresh fits; the source also contains its native-sampling comparison mode. See `ENSEMBLE_DESIGN.md` for seeds, masks and starts.
- `ensemble_summary.py` derives retained-observation frequencies and pair-order counts.
- `ensemble_competition.py` compares recovered objectives only within identical samples and penalties.
- `ensemble_findings.py` writes the compact report tables, start-family outcomes and page-versus-sheet comparisons.
- `ensemble_audit.py` verifies the completed fits and accounting; the independent weighted replay is separate.
- `ensemble_independent_audit.py` makes exactly two additional growth calls and independently checks their saved nodes, positions and energy components.
- `cassidy_h_continuation.py`, then `cassidy_summary.py`, run the event-aware equation-8 continuation and summarize it. `cassidy_independent_audit.py` rechecks derivatives and selected paths independently. `cassidy_return_check.py` is the additional excursion-return diagnostic.
- `prepare_scegret.py` recreates the observed feature matrices and candidate orders without a new curve fit.

## Original scEGRET R functions

The authors' original R sources, license and pinned commit/checksums are under `scegret_source`. Only `egret_flap.R` and `egret_fly.R` are sourced to execute the assessment; they are not edited. The calculation requires methods, matrixStats, matrixTests, reshape2, dplyr and jsonlite. Exact versions are recorded in `scegret_results.json`.

This run used Node 24.15.0 and WebR 0.6.0 (R 4.6.0), because a native R interpreter was unavailable. To restore the temporary interpreter, install `webr@0.6.0` under `_r_runtime` with npm. On Windows, invoke Node with `--import ./webr_windows_register.mjs ./run_scegret.mjs`. The import adapter converts Windows filesystem paths to file URLs; it does not alter R or scEGRET calculations. R packages are downloaded from the WebR package repository into the local interpreter; manuscript data are processed locally.

Running without a segment argument executes all 648 cases. The recorded run split independent feature batches across workers; `merge_scegret.py` checks that the retained partitions cover every case once. `normalize_scegret_output.py` preserves factor codes and corrects original-bin labels in the initial wrapper's output. Its correction changes labels only, not scores or bin membership. The current wrapper already exports original labels correctly. `scegret_independent_audit.py` recomputes every continuity score independently from the CSV inputs and checks the R outputs.

## Figures and report

For the requested actual-page removal follow-up, run `python calibration/omission_test/run_omissions.py`, then `python calibration/omission_test/build_report.py`. It regenerates 7,596 frozen continuous-reference draws and applies seven masks to each draw, keeping all original inputs and calibration outputs unchanged. Its separate `DESIGN.md` defines the targeted and whole-path comparisons, broad-section definitions and rebinning behavior. The original reference means and SDs must reproduce exactly before results are accepted.

The completed scEGRET calibration is under `calibration/`. Read its `DESIGN.md`, `MATCHED_DESIGN.md` and `BYPASS_DESIGN.md`; the latter two explicitly record additions made after preliminary results were inspected. After verifying the original R outputs with `scegret_independent_audit.py`, run `calibration/run_calibration.py`, `calibration/matched_power.py`, `calibration/bypass_calibration.py`, `calibration/summarize.py`, and `calibration/build_report.py` in that order. The 147,888 production synthetic datasets are streamed, not retained. The Monte Carlo score is an exact NumPy implementation validated against the original R evaluations; no R runtime or new research dependency is required for this follow-up. Complete source, input hashes, null reference distributions, false-alarm checks, power conditions and numerical results are included. These calibrations condition on the original Golf order; they do not refit every alternative trajectory.

Run `plot_cassidy_h.py`, `plot_ordering_ensemble.py`, `plot_ordering_ensemble.py sheet_control`, `scegret_summary.py`, and `plot_scegret.py` against the verified results. Plotting does not refit models. `build_gallery.py` converts `RESULTS.md` to a local HTML report, writes the image bundle, and checks local links and PNG integrity. `finalize_manifest.py` records the verified release hashes.

Temporary package directories and duplicate partition checkpoints can be removed after verification; retain the source, final merged metrics, hashes and small feature matrices. Run `python tools/artifact_storage.py --check` from the repository root as the final storage gate.
