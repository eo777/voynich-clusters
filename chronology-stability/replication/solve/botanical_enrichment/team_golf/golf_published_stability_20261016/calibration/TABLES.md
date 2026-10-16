# Calibration numerical tables

All tails below are conditional simulation references. “Adjusted” scans boundaries, three bin counts, two scales, and applies a factor of four across channels/cohorts. They are not historical probabilities.

## H boundaries at primary bandwidth 0.12, sheet-resampling model

| Cohort | Channel | Bins | Observed | Continuous mean | Continuous 95% range | Local tail | Adjusted tail |
|---|---|---|---:|---:|---|---:|---:|
| paragraph207 | line | 7 → 8 | 0.279 | 0.217 | 0.162–0.275 | 0.018 | 1.000 |
| paragraph207 | line | 8 → 9 | 0.286 | 0.217 | 0.156–0.274 | 0.010 | 1.000 |
| paragraph207 | word | 7 → 8 | 0.271 | 0.221 | 0.163–0.284 | 0.059 | 1.000 |
| paragraph207 | word | 8 → 9 | 0.270 | 0.222 | 0.163–0.291 | 0.073 | 1.000 |
| pool50 | line | 5 → 6 | 0.258 | 0.202 | 0.150–0.265 | 0.040 | 1.000 |
| pool50 | line | 6 → 7 | 0.180 | 0.170 | 0.121–0.224 | 0.306 | 1.000 |
| pool50 | word | 5 → 6 | 0.251 | 0.219 | 0.156–0.275 | 0.101 | 1.000 |
| pool50 | word | 6 → 7 | 0.170 | 0.176 | 0.123–0.235 | 0.575 | 1.000 |

## Supplementary two-sided excursion diagnostic

B = min(E entering, E leaving) − E bypass. This is our diagnostic, not a published scEGRET test. Adjusted tails also include a factor of two for the additional statistic family.

| Cohort | Channel | Model | Bandwidth | Bins | Scale | Middle bin | Observed B | Local tail | Adjusted tail |
|---|---|---|---:|---:|---|---:|---:|---:|---:|
| paragraph207 | line | dirichlet_multinomial | 0.12 | 10 | per_million | 7 | 0.123 | 0.349 | 1.000 |
| paragraph207 | line | dirichlet_multinomial | 0.12 | 12 | per_million | 8 | 0.102 | 0.467 | 1.000 |
| paragraph207 | line | dirichlet_multinomial | 0.12 | 14 | per_million | 9 | -0.055 | 0.450 | 1.000 |
| paragraph207 | line | sheet_residual | 0.12 | 10 | per_million | 7 | 0.123 | 0.155 | 1.000 |
| paragraph207 | line | sheet_residual | 0.12 | 12 | per_million | 8 | 0.102 | 0.216 | 1.000 |
| paragraph207 | line | sheet_residual | 0.12 | 14 | per_million | 9 | -0.055 | 0.495 | 1.000 |
| paragraph207 | word | dirichlet_multinomial | 0.12 | 10 | per_million | 7 | 0.122 | 0.391 | 1.000 |
| paragraph207 | word | dirichlet_multinomial | 0.12 | 12 | per_million | 8 | 0.101 | 0.545 | 1.000 |
| paragraph207 | word | dirichlet_multinomial | 0.12 | 14 | per_million | 9 | -0.072 | 0.776 | 1.000 |
| paragraph207 | word | sheet_residual | 0.12 | 10 | per_million | 7 | 0.122 | 0.177 | 1.000 |
| paragraph207 | word | sheet_residual | 0.12 | 12 | per_million | 8 | 0.101 | 0.243 | 1.000 |
| paragraph207 | word | sheet_residual | 0.12 | 14 | per_million | 9 | -0.072 | 0.681 | 1.000 |
| pool50 | line | dirichlet_multinomial | 0.12 | 6 | per_million | 4 | -0.103 | 0.852 | 1.000 |
| pool50 | line | dirichlet_multinomial | 0.12 | 8 | per_million | 6 | -0.041 | 0.112 | 1.000 |
| pool50 | line | sheet_residual | 0.12 | 6 | per_million | 4 | -0.103 | 0.771 | 1.000 |
| pool50 | line | sheet_residual | 0.12 | 8 | per_million | 6 | -0.041 | 0.355 | 1.000 |
| pool50 | word | dirichlet_multinomial | 0.12 | 6 | per_million | 4 | -0.114 | 0.892 | 1.000 |
| pool50 | word | dirichlet_multinomial | 0.12 | 8 | per_million | 6 | -0.058 | 0.412 | 1.000 |
| pool50 | word | sheet_residual | 0.12 | 6 | per_million | 4 | -0.114 | 0.768 | 1.000 |
| pool50 | word | sheet_residual | 0.12 | 8 | per_million | 6 | -0.058 | 0.514 | 1.000 |

## Other boundaries flagged under sheet resampling

All qualifying results are retained, including scale and bandwidth sensitivity. These are deviations from a specified smooth reference, not proof of missing text.

| Cohort | Channel | Bandwidth | Bins | Scale | Original bins | Score | Adjusted tail |
|---|---|---:|---:|---|---|---:|---:|
| paragraph207 | line | 0.08 | 10 | per_million | 2 → 3 | 0.171 | 0.032 |
| paragraph207 | line | 0.08 | 12 | raw | 4 → 5 | 0.426 | 0.020 |
| paragraph207 | line | 0.12 | 10 | per_million | 2 → 3 | 0.171 | 0.048 |
| paragraph207 | line | 0.12 | 12 | raw | 4 → 5 | 0.426 | 0.028 |
| paragraph207 | line | 0.20 | 10 | per_million | 2 → 3 | 0.171 | 0.044 |
| paragraph207 | line | 0.20 | 12 | raw | 4 → 5 | 0.426 | 0.024 |
| paragraph207 | word | 0.08 | 10 | per_million | 2 → 3 | 0.183 | 0.024 |
| paragraph207 | word | 0.08 | 10 | raw | 2 → 3 | 0.317 | 0.008 |
| paragraph207 | word | 0.12 | 10 | per_million | 2 → 3 | 0.183 | 0.024 |
| paragraph207 | word | 0.12 | 10 | raw | 2 → 3 | 0.317 | 0.020 |
| paragraph207 | word | 0.12 | 12 | per_million | 10 → 11 | 0.135 | 0.048 |
| paragraph207 | word | 0.20 | 10 | per_million | 2 → 3 | 0.183 | 0.024 |
| paragraph207 | word | 0.20 | 10 | raw | 2 → 3 | 0.317 | 0.024 |
| paragraph207 | word | 0.20 | 10 | raw | 4 → 5 | 0.393 | 0.048 |
| paragraph207 | word | 0.20 | 12 | per_million | 10 → 11 | 0.135 | 0.044 |
| paragraph207 | word | 0.20 | 14 | per_million | 3 → 4 | 0.201 | 0.048 |
| pool50 | line | 0.08 | 8 | raw | 3 → 4 | 0.507 | 0.024 |
| pool50 | line | 0.08 | 10 | raw | 3 → 5 | 0.541 | 0.024 |
| pool50 | line | 0.12 | 6 | raw | 2 → 3 | 0.462 | 0.020 |
| pool50 | line | 0.12 | 8 | raw | 3 → 4 | 0.507 | 0.016 |
| pool50 | line | 0.12 | 10 | raw | 3 → 5 | 0.541 | 0.016 |
| pool50 | line | 0.20 | 6 | raw | 2 → 3 | 0.462 | 0.008 |
| pool50 | line | 0.20 | 8 | raw | 3 → 4 | 0.507 | 0.008 |
| pool50 | line | 0.20 | 10 | raw | 3 → 5 | 0.541 | 0.008 |
| pool50 | word | 0.08 | 10 | raw | 2 → 3 | 0.392 | 0.028 |
| pool50 | word | 0.12 | 8 | raw | 2 → 3 | 0.348 | 0.048 |
| pool50 | word | 0.12 | 10 | raw | 2 → 3 | 0.392 | 0.032 |

## Model-sensitive H-neighborhood flags

These are all H-touching boundaries that pass the adjusted screen under any sensitivity model. They must be read alongside results under empirical sheet resampling and per-million scaling.

| Cohort | Channel | Model | Bandwidth | Bins | Scale | Boundary | Adjusted tail |
|---|---|---|---:|---:|---|---|---:|
| pool50 | line | dirichlet_multinomial | 0.08 | 8 | raw | 6 → 7 | 0.024 |
| pool50 | line | dirichlet_multinomial | 0.12 | 8 | raw | 6 → 7 | 0.016 |
| pool50 | line | dirichlet_multinomial | 0.20 | 8 | raw | 5 → 6 | 0.008 |
| pool50 | line | dirichlet_multinomial | 0.20 | 8 | raw | 6 → 7 | 0.028 |
| pool50 | word | dirichlet_multinomial | 0.08 | 8 | raw | 6 → 7 | 0.012 |
| pool50 | word | dirichlet_multinomial | 0.12 | 8 | raw | 6 → 7 | 0.012 |
| pool50 | word | dirichlet_multinomial | 0.20 | 8 | raw | 6 → 7 | 0.004 |

## Exact decision-rule power: gap join at 0.65

| Cohort | Channel | Noise model | Gap width | Detected / 300 | Detection rate (95% interval) | Continuous false alarms / 500 |
|---|---|---|---:|---:|---|---:|
| paragraph207 | line | dirichlet_multinomial | 10% | 30 | 10.0% (7.1%–13.9%) | 5 |
| paragraph207 | line | dirichlet_multinomial | 20% | 140 | 46.7% (41.1%–52.3%) | 5 |
| paragraph207 | line | dirichlet_multinomial | 30% | 263 | 87.7% (83.5%–90.9%) | 5 |
| paragraph207 | line | sheet_residual | 10% | 1 | 0.3% (0.1%–1.9%) | 3 |
| paragraph207 | line | sheet_residual | 20% | 3 | 1.0% (0.3%–2.9%) | 3 |
| paragraph207 | line | sheet_residual | 30% | 5 | 1.7% (0.7%–3.8%) | 3 |
| paragraph207 | word | dirichlet_multinomial | 10% | 12 | 4.0% (2.3%–6.9%) | 4 |
| paragraph207 | word | dirichlet_multinomial | 20% | 63 | 21.0% (16.8%–26.0%) | 4 |
| paragraph207 | word | dirichlet_multinomial | 30% | 161 | 53.7% (48.0%–59.2%) | 4 |
| paragraph207 | word | sheet_residual | 10% | 4 | 1.3% (0.5%–3.4%) | 5 |
| paragraph207 | word | sheet_residual | 20% | 6 | 2.0% (0.9%–4.3%) | 5 |
| paragraph207 | word | sheet_residual | 30% | 4 | 1.3% (0.5%–3.4%) | 5 |
| pool50 | line | dirichlet_multinomial | 10% | 9 | 3.0% (1.6%–5.6%) | 10 |
| pool50 | line | dirichlet_multinomial | 20% | 77 | 25.7% (21.1%–30.9%) | 10 |
| pool50 | line | dirichlet_multinomial | 30% | 193 | 64.3% (58.8%–69.5%) | 10 |
| pool50 | line | sheet_residual | 10% | 9 | 3.0% (1.6%–5.6%) | 9 |
| pool50 | line | sheet_residual | 20% | 10 | 3.3% (1.8%–6.0%) | 9 |
| pool50 | line | sheet_residual | 30% | 25 | 8.3% (5.7%–12.0%) | 9 |
| pool50 | word | dirichlet_multinomial | 10% | 12 | 4.0% (2.3%–6.9%) | 5 |
| pool50 | word | dirichlet_multinomial | 20% | 49 | 16.3% (12.6%–20.9%) | 5 |
| pool50 | word | dirichlet_multinomial | 30% | 133 | 44.3% (38.8%–50.0%) | 5 |
| pool50 | word | sheet_residual | 10% | 5 | 1.7% (0.7%–3.8%) | 3 |
| pool50 | word | sheet_residual | 20% | 9 | 3.0% (1.6%–5.6%) | 3 |
| pool50 | word | sheet_residual | 30% | 12 | 4.0% (2.3%–6.9%) | 3 |

## Literal removal versus same-size random thinning

This uses the paper-style largest adjacent CPM score at the middle bin count, with a separate same-size random-deletion reference. It does not use the standardized multi-view decision rule above.

| Cohort | Channel | Model | Center | Span | Removed | Gap detection | Random deletion false alarms |
|---|---|---|---:|---:|---:|---:|---:|
| paragraph207 | line | dirichlet_multinomial | 0.35 | 10% | 16 | 6.3% | 5.0% |
| paragraph207 | line | dirichlet_multinomial | 0.35 | 20% | 44 | 15.0% | 6.3% |
| paragraph207 | line | dirichlet_multinomial | 0.35 | 30% | 66 | 15.3% | 3.7% |
| paragraph207 | line | dirichlet_multinomial | 0.65 | 10% | 5 | 5.0% | 4.3% |
| paragraph207 | line | dirichlet_multinomial | 0.65 | 20% | 34 | 12.3% | 4.0% |
| paragraph207 | line | dirichlet_multinomial | 0.65 | 30% | 72 | 35.3% | 6.7% |
| paragraph207 | line | dirichlet_multinomial | 0.85 | 10% | 29 | 9.7% | 3.7% |
| paragraph207 | line | dirichlet_multinomial | 0.85 | 20% | 54 | 10.7% | 4.7% |
| paragraph207 | line | dirichlet_multinomial | 0.85 | 30% | 82 | 21.3% | 7.0% |
| paragraph207 | line | sheet_residual | 0.35 | 10% | 16 | 13.7% | 8.0% |
| paragraph207 | line | sheet_residual | 0.35 | 20% | 44 | 14.0% | 4.3% |
| paragraph207 | line | sheet_residual | 0.35 | 30% | 66 | 18.3% | 4.7% |
| paragraph207 | line | sheet_residual | 0.65 | 10% | 5 | 6.0% | 5.7% |
| paragraph207 | line | sheet_residual | 0.65 | 20% | 34 | 6.7% | 2.3% |
| paragraph207 | line | sheet_residual | 0.65 | 30% | 72 | 26.3% | 5.3% |
| paragraph207 | line | sheet_residual | 0.85 | 10% | 29 | 10.7% | 9.3% |
| paragraph207 | line | sheet_residual | 0.85 | 20% | 54 | 11.0% | 4.3% |
| paragraph207 | line | sheet_residual | 0.85 | 30% | 82 | 22.7% | 8.0% |
| paragraph207 | word | dirichlet_multinomial | 0.35 | 10% | 16 | 10.3% | 8.0% |
| paragraph207 | word | dirichlet_multinomial | 0.35 | 20% | 44 | 7.7% | 3.0% |
| paragraph207 | word | dirichlet_multinomial | 0.35 | 30% | 66 | 13.7% | 5.0% |
| paragraph207 | word | dirichlet_multinomial | 0.65 | 10% | 5 | 4.7% | 4.3% |
| paragraph207 | word | dirichlet_multinomial | 0.65 | 20% | 34 | 13.3% | 5.0% |
| paragraph207 | word | dirichlet_multinomial | 0.65 | 30% | 72 | 36.0% | 6.3% |
| paragraph207 | word | dirichlet_multinomial | 0.85 | 10% | 29 | 4.7% | 3.0% |
| paragraph207 | word | dirichlet_multinomial | 0.85 | 20% | 54 | 8.7% | 3.7% |
| paragraph207 | word | dirichlet_multinomial | 0.85 | 30% | 82 | 25.7% | 9.7% |
| paragraph207 | word | sheet_residual | 0.35 | 10% | 16 | 2.3% | 1.0% |
| paragraph207 | word | sheet_residual | 0.35 | 20% | 44 | 14.3% | 4.7% |
| paragraph207 | word | sheet_residual | 0.35 | 30% | 66 | 14.7% | 4.3% |
| paragraph207 | word | sheet_residual | 0.65 | 10% | 5 | 6.3% | 4.7% |
| paragraph207 | word | sheet_residual | 0.65 | 20% | 34 | 10.0% | 4.7% |
| paragraph207 | word | sheet_residual | 0.65 | 30% | 72 | 16.7% | 3.0% |
| paragraph207 | word | sheet_residual | 0.85 | 10% | 29 | 7.3% | 3.7% |
| paragraph207 | word | sheet_residual | 0.85 | 20% | 54 | 12.0% | 5.0% |
| paragraph207 | word | sheet_residual | 0.85 | 30% | 82 | 14.7% | 4.3% |
| pool50 | line | dirichlet_multinomial | 0.35 | 10% | 3 | 16.3% | 6.0% |
| pool50 | line | dirichlet_multinomial | 0.35 | 20% | 11 | 0.0% | 2.0% |
| pool50 | line | dirichlet_multinomial | 0.35 | 30% | 16 | 0.0% | 6.7% |
| pool50 | line | dirichlet_multinomial | 0.65 | 10% | 2 | 15.7% | 6.0% |
| pool50 | line | dirichlet_multinomial | 0.65 | 20% | 9 | 0.3% | 4.0% |
| pool50 | line | dirichlet_multinomial | 0.65 | 30% | 16 | 1.3% | 3.7% |
| pool50 | line | dirichlet_multinomial | 0.85 | 10% | 5 | 0.3% | 4.7% |
| pool50 | line | dirichlet_multinomial | 0.85 | 20% | 12 | 0.0% | 4.0% |
| pool50 | line | dirichlet_multinomial | 0.85 | 30% | 19 | 0.0% | 4.7% |
| pool50 | line | sheet_residual | 0.35 | 10% | 3 | 17.0% | 4.7% |
| pool50 | line | sheet_residual | 0.35 | 20% | 11 | 0.0% | 8.3% |
| pool50 | line | sheet_residual | 0.35 | 30% | 16 | 0.0% | 5.0% |
| pool50 | line | sheet_residual | 0.65 | 10% | 2 | 4.3% | 5.0% |
| pool50 | line | sheet_residual | 0.65 | 20% | 9 | 2.7% | 5.3% |
| pool50 | line | sheet_residual | 0.65 | 30% | 16 | 8.7% | 5.7% |
| pool50 | line | sheet_residual | 0.85 | 10% | 5 | 1.3% | 5.3% |
| pool50 | line | sheet_residual | 0.85 | 20% | 12 | 0.0% | 4.7% |
| pool50 | line | sheet_residual | 0.85 | 30% | 19 | 0.0% | 6.7% |
| pool50 | word | dirichlet_multinomial | 0.35 | 10% | 3 | 10.7% | 5.7% |
| pool50 | word | dirichlet_multinomial | 0.35 | 20% | 11 | 0.0% | 9.0% |
| pool50 | word | dirichlet_multinomial | 0.35 | 30% | 16 | 0.0% | 4.3% |
| pool50 | word | dirichlet_multinomial | 0.65 | 10% | 2 | 5.7% | 5.0% |
| pool50 | word | dirichlet_multinomial | 0.65 | 20% | 9 | 0.3% | 4.7% |
| pool50 | word | dirichlet_multinomial | 0.65 | 30% | 16 | 0.3% | 4.7% |
| pool50 | word | dirichlet_multinomial | 0.85 | 10% | 5 | 0.0% | 5.0% |
| pool50 | word | dirichlet_multinomial | 0.85 | 20% | 12 | 0.0% | 4.7% |
| pool50 | word | dirichlet_multinomial | 0.85 | 30% | 19 | 0.0% | 8.0% |
| pool50 | word | sheet_residual | 0.35 | 10% | 3 | 8.0% | 3.0% |
| pool50 | word | sheet_residual | 0.35 | 20% | 11 | 0.0% | 3.0% |
| pool50 | word | sheet_residual | 0.35 | 30% | 16 | 0.0% | 5.3% |
| pool50 | word | sheet_residual | 0.65 | 10% | 2 | 2.3% | 5.0% |
| pool50 | word | sheet_residual | 0.65 | 20% | 9 | 4.3% | 4.7% |
| pool50 | word | sheet_residual | 0.65 | 30% | 16 | 14.3% | 3.7% |
| pool50 | word | sheet_residual | 0.85 | 10% | 5 | 2.0% | 3.3% |
| pool50 | word | sheet_residual | 0.85 | 20% | 12 | 0.3% | 6.7% |
| pool50 | word | sheet_residual | 0.85 | 30% | 19 | 0.0% | 9.0% |

## Null-model diagnostics

Zero fractions refer to the channel’s features observed at least once anywhere; all-zero features have no effect on EGRET. Diagnostics are essential: a conditional calibration is only as credible as its generator.

| Cohort | Channel | Model | Bandwidth | Detected features: data / median simulated | Zero fraction: data / median simulated | Residual squared error: data / median simulated | Held-out 5% scan false alarms |
|---|---|---|---:|---|---|---|---:|
| paragraph207 | line | dirichlet_multinomial | 0.08 | 249 / 208 | 0.761 / 0.798 | 0.00480 / 0.00525 | 5.2% |
| paragraph207 | line | dirichlet_multinomial | 0.12 | 249 / 206 | 0.761 / 0.799 | 0.00517 / 0.00551 | 4.4% |
| paragraph207 | line | dirichlet_multinomial | 0.20 | 249 / 204 | 0.761 / 0.799 | 0.00566 / 0.00594 | 4.0% |
| paragraph207 | line | sheet_residual | 0.08 | 249 / 247 | 0.761 / 0.746 | 0.00480 / 0.00455 | 6.6% |
| paragraph207 | line | sheet_residual | 0.12 | 249 / 244 | 0.761 / 0.747 | 0.00517 / 0.00479 | 3.4% |
| paragraph207 | line | sheet_residual | 0.20 | 249 / 240 | 0.761 / 0.749 | 0.00566 / 0.00515 | 4.6% |
| paragraph207 | word | dirichlet_multinomial | 0.08 | 214 / 172 | 0.774 / 0.814 | 0.00442 / 0.00486 | 5.2% |
| paragraph207 | word | dirichlet_multinomial | 0.12 | 214 / 171 | 0.774 / 0.814 | 0.00477 / 0.00509 | 7.0% |
| paragraph207 | word | dirichlet_multinomial | 0.20 | 214 / 169 | 0.774 / 0.814 | 0.00522 / 0.00548 | 5.6% |
| paragraph207 | word | sheet_residual | 0.08 | 214 / 216 | 0.774 / 0.760 | 0.00442 / 0.00432 | 4.8% |
| paragraph207 | word | sheet_residual | 0.12 | 214 / 214 | 0.774 / 0.760 | 0.00477 / 0.00453 | 3.6% |
| paragraph207 | word | sheet_residual | 0.20 | 214 / 210 | 0.774 / 0.763 | 0.00522 / 0.00490 | 5.4% |
| pool50 | line | dirichlet_multinomial | 0.08 | 347 / 266 | 0.601 / 0.701 | 0.00204 / 0.00255 | 6.2% |
| pool50 | line | dirichlet_multinomial | 0.12 | 347 / 262 | 0.601 / 0.704 | 0.00242 / 0.00279 | 5.4% |
| pool50 | line | dirichlet_multinomial | 0.20 | 347 / 258 | 0.601 / 0.706 | 0.00286 / 0.00314 | 4.8% |
| pool50 | line | sheet_residual | 0.08 | 347 / 354 | 0.601 / 0.583 | 0.00204 / 0.00218 | 3.6% |
| pool50 | line | sheet_residual | 0.12 | 347 / 352 | 0.601 / 0.583 | 0.00242 / 0.00245 | 3.0% |
| pool50 | line | sheet_residual | 0.20 | 347 / 346 | 0.601 / 0.587 | 0.00286 / 0.00279 | 4.2% |
| pool50 | word | dirichlet_multinomial | 0.08 | 306 / 222 | 0.625 / 0.731 | 0.00192 / 0.00242 | 6.0% |
| pool50 | word | dirichlet_multinomial | 0.12 | 306 / 218 | 0.625 / 0.734 | 0.00228 / 0.00264 | 4.8% |
| pool50 | word | dirichlet_multinomial | 0.20 | 306 / 214 | 0.625 / 0.737 | 0.00270 / 0.00295 | 2.0% |
| pool50 | word | sheet_residual | 0.08 | 306 / 312 | 0.625 / 0.610 | 0.00192 / 0.00214 | 6.4% |
| pool50 | word | sheet_residual | 0.12 | 306 / 309 | 0.625 / 0.613 | 0.00228 / 0.00236 | 4.0% |
| pool50 | word | sheet_residual | 0.20 | 306 / 303 | 0.625 / 0.615 | 0.00270 / 0.00270 | 6.0% |
