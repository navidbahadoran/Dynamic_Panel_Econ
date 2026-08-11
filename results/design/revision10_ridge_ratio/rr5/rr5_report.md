# RR5 fresh Revision-10 rank-selection preflight

## Scope and integrity

- Approved implementation commit: `24ae1b5817bedc405c8debf12d6b6f9bf8cb9375`
- RR4 implementation commit: `7325c323bd98f6ccf7bf89f4f29e0841c9fec7b9`
- Locked manuscript commit: `47af37560d6e0752519d6a6267a152871d8a5157`
- Master seed: `2026081101`
- Attempts: 4 DGPs × 5 panel shapes × 3 true-rank vectors × 3 replications = 180
- Workers: requested 12, effective 12, one native numerical-library thread per worker
- Targets/inference: none; no Riesz, splits, standard errors, or inferential targets were run
- Frozen calibration SHA-256: `51b4ea76fdd2295a2cfeb840ef7b2b3f468d32c45ad2c123523363d5b968f17e`
- RR5 config SHA-256: `7175c86bf76f2b276bdc62b5d0c491dd2100e18274105050a6ba0578cf524e80`
- Locked manuscript SHA-256: `d1c7224e0962df6a54fb92a550671f72b66258f30afd64d292f3bf60c7d87be9`

## Accounting

The completed evidence contains 180 attempt rows, 180 distinct semantic replication IDs, 180 distinct realization hashes, 60 calibration/rank-design cells, and exactly three pilot-start diagnostics per attempt (540 start rows). All 180 statuses are `rank_selection_numerically_unresolved`. There are zero accepted pilots, zero selected rank vectors, zero cap selections, and zero final post-refits.

The first execution wrapper expired after 7,200 seconds with 39 atomic chunks (117 attempts) safely checkpointed. The resume-safe continuation skipped those chunks and completed the remaining 21 chunks (63 attempts) without changing any scientific setting. This operational interruption is recorded in `rr5_interruption_log.json`; it does not create duplicate semantic attempts.

## Results by true rank

| true_rank_vector | attempted_replications | numerically_acceptable_spectral_pilots | exact_recovery | underselection_only | overselection_only | mixed_rank_errors | unresolved_selections | final_postrefit_valid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| (1,0,2) | 60 | 0 | 0 | 0 | 0 | 0 | 60 | 0 |
| (1,1,1) | 60 | 0 | 0 | 0 | 0 | 0 | 60 | 0 |
| (2,1,1) | 60 | 0 | 0 | 0 | 0 | 0 | 60 | 0 |

## Results by DGP

| dgp | attempted_replications | numerically_acceptable_spectral_pilots | exact_recovery | underselection_only | overselection_only | mixed_rank_errors | unresolved_selections | final_postrefit_valid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 45 | 0 | 0 | 0 | 0 | 0 | 45 | 0 |
| 2 | 45 | 0 | 0 | 0 | 0 | 0 | 45 | 0 |
| 3 | 45 | 0 | 0 | 0 | 0 | 0 | 45 | 0 |
| 4 | 45 | 0 | 0 | 0 | 0 | 0 | 45 | 0 |

## Results by panel dimensions

| N | T | attempted_replications | numerically_acceptable_spectral_pilots | exact_recovery | underselection_only | overselection_only | mixed_rank_errors | unresolved_selections | final_postrefit_valid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50 | 50 | 36 | 0 | 0 | 0 | 0 | 0 | 36 | 0 |
| 100 | 100 | 36 | 0 | 0 | 0 | 0 | 0 | 36 | 0 |
| 100 | 200 | 36 | 0 | 0 | 0 | 0 | 0 | 36 | 0 |
| 200 | 100 | 36 | 0 | 0 | 0 | 0 | 0 | 36 | 0 |
| 200 | 200 | 36 | 0 | 0 | 0 | 0 | 0 | 36 | 0 |

## Numerical diagnostics

- Converged/stationarity-passing starts: 53 / 540
- Non-passing starts: 487 / 540
- Attempts with at least one stationarity-passing start: 50 / 180
- Attempts with at least two stationarity-passing starts: 3 / 180
- Attempts accepted under the complete frozen cap+1 diagnostics: 0 / 180
- Median stationarity residual: `1.32622364e-05`
- 95th-percentile stationarity residual: `0.00356981726`
- Maximum stationarity residual: `0.02679843`
- Starts invoking constrained fallback: 201 / 540
- Starts marked boundary-active: 180 / 540
- Maximum recorded constraint violation: `3.55271368e-15`
- Best-pilot stationarity failures: 166 / 180
- Objective-stability failures: 180 / 180
- Feasibility failures: 0 / 180
- Nonfinite failures: 0 / 180
- Boundary-active pilot selections: 79 / 180
- Median across-start objective spread: `0.00161073475`
- 95th-percentile across-start objective spread: `0.0117848741`

Each dedicated truth-block audit contains 60 replication rows. Because all pilots were unresolved, zero selected B ranks are observed for truth `(1,0,2)`, zero selected A ranks are observed for truth `(2,1,1)`, and zero selected H ranks are observed for truth `(1,0,2)`; all 180 block-audit outcomes are marked not reached. The cap-selection audit contains zero instances. These are numerical-unresolved counts, not evidence of statistical rank-zero or rank-two selection failure.

The rectangular summaries remain separate: `(100,200)` and `(200,100)` each contain 36 attempts, zero valid pilots, and 36 unresolved selections. Block-recovery quantities are unavailable because no selector was reached.

Normalized pilot spectra and ridge ratios are deliberately unavailable: the frozen implementation computes them only after accepting the joint pilot. The audit tables mark these fields `not_reached_due_to_pilot_failure`; they do not fabricate singular values, ratios, selected ranks, or recovery classifications.

## Manuscript-report compatibility

The evidence schema supplies the quantities requested by the locked Revision-10 Monte Carlo design: attempted replications, acceptable pilots, exact/under/over/mixed/unresolved counts, cap selections, selected-rank distributions, block diagnostics, normalized spectra, ridge ratios, ratio margins, and rectangular-panel summaries. For this NO-GO run, only attempted and unresolved counts plus pilot-start diagnostics are observed. Selector-dependent columns remain explicitly unavailable. This is compatible with the manuscript's numerical-unresolved rule and its requirement to count failure reasons.

## Post-run validation

- `pytest -q`: 167 passed in 58.10 seconds
- `ruff check src scripts tests`: passed
- `git diff --check`: passed
- Scientific source diff from RR4 under `src/` and `scripts/`: empty
- Only RR5 evidence is uncommitted

## Decision

**RR5-NO-GO.** The fresh evidence does not support proceeding to larger rank-selection experiments. No scientific source, manuscript, calibration value, selector definition, or numerical tolerance was changed, and no post-hoc rerun or tuning was performed.
