# Published-method replication for Golf uncertainty and H influence

Started 16 October 2026. The requested outputs are a combined page/bifolium ordering-sensitivity figure, a two-plane H-position influence figure at lambda 0.01, and a report connecting implementation to the original research.

## Methods and scope

1. Reproduce ElPiGraph's principal-graph ensemble procedure (Albergante et al., 2020, section 2.10, Figure 1C): random 90% subsamples and 100 realizations, with declared multiple starts and parameter comparisons. Use the actual published package for graph growth. Document native observation resampling versus physical-bifolium block adaptation, fixed geometry, weights, initialization, and the frequency interpretation of the display before the main run.
2. Implement Cassidy (2023), section 2.2, equations 7-8: implicit-function sensitivity from the inverse Hessian and mixed derivative, then predictor/corrector continuation for H-coordinate changes. The objective is smooth only while nearest-node assignments remain fixed. Detect crossings, correct and relinearize; do not extend a local derivative across an assignment event without checking it. Validate derivatives, prediction/correction errors, direction dependence, and the resulting preference boundary.
3. Keep historical probability, objective preference, route recovery, and local solution persistence distinct. No claim of chronology correctness, complete basin enumeration, or a temporal Lyapunov exponent is intended.

References:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC7516753/ (ElPiGraph ensemble procedure).
- https://github.com/j-bac/elpigraph-python (published implementation).
- https://link.springer.com/article/10.1007/s11538-023-01200-0 (Cassidy continuation, including Remark 1 on general smooth objectives).

## Resource estimate

The workspace gate passed at approximately 18.22 GiB before this study. Expected temporary dependency footprint is below 300 MB. Compact metrics, projected positions, charts and selected representatives should remain below 100 MB. No full ensemble of high-dimensional node arrays will be retained. Expected peak additional storage is below 0.5 GiB, within the 20 GiB budget. Dependencies will be removed after numerical verification; source, report, compact metrics and hashes will remain. Run `python tools/artifact_storage.py --check` at completion.

Existing studies and figures remain unchanged until replacement outputs are verified. No private or unpublished paper is needed for this analysis.

## Added scEGRET continuity analysis

The user explicitly requested the scEGRET method during this run. Execute its original `egret.flap` and `egret.fly` R functions from pinned public source commit `93517e5ce3f09ccdd3c0d604d0a1b53e70c5cc70`, with text-transition counts or declared per-profile normalized transition frequencies taking the feature-matrix role. Preserve its source-defined binning, detection filtering, mean/log normalization and enthalpy calculation. Document differences between source defaults and prose descriptions, sparse-bin exclusions, multiple bin resolutions and alternative trajectories. Do not substitute JS distances into a renamed EGRET statistic. A temporary local WebR runtime is planned because R is not installed; anticipated additional storage is below 200 MB and will be removed after verification.
