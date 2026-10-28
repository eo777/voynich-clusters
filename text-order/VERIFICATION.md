# Verification of the text entry points

Checked on 28 October 2026 using Python 3.9, NumPy 2.0.2, SciPy 1.13.1,
ElPiGraph 0.3.4 and the pinned scientific dependencies.

Eight fresh raw-text replays (sample 0, start 0 in each of the two
representations × two cohorts × two penalties) converged. Every projected
position, total objective, node MSE, edge energy and star energy matched the
corresponding frozen fit **exactly** in this environment (maximum absolute
error 0). This is a smoke test, not a rerun of all 6,400 fits.

All eight reconstructed distance matrices, MDS coordinate arrays, 100-sample
masks and all 64 initial node configurations also matched the published
design byte-for-byte. A guarded reconstruction rejected reads of saved JSON,
NPZ and CSV scientific inputs/results in the research tree; reconstruction
completed without those reads.

The portable tests check unreadable-word barriers, preservation of line
boundaries, invariance of counts to reporting labels and input row order,
pairwise ties/reversal and missing-pair rejection. They also perform real
ElPiGraph growth on an artificial curved dataset, verify the requested node
count and convergence, and check closest-segment positions under reversal.
The included artificial text CSV was run through the complete text-to-order
command with two starts, eight nodes, both penalties and three full-inclusion
samples. Full inclusion is an installation check, not evidence of sampling
stability or recovery of historical chronology.

The exported repository was tested separately from the source workspace.
Its portable tests passed, and its IVTFF custom-text command fitted all 50 ZL
pooled bifolia with two starts, both penalties, 28 nodes and two full-inclusion
samples. Q13 (`M`) and Q20 (`T`) were present in the post-fit group summaries.
This checks the executable workflow; it is not a full unlabeled ensemble study.

The existing v2 tests still pass: saved-fit replay agrees with the published
inspector, conditional selection precedes best-fit choice, and missingness
and exported scores retain their original meaning. No published score or
classification was changed by adding these entry points.

Readers can run `verify_refit.py` against either a smoke test or a full refit.
That separate verification step reads the frozen reference outcomes; the
raw-text fitting commands do not.
