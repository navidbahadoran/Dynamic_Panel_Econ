# RR5e fresh post-engineering scientific preflight

## Integrity, seed, and design

- Pre-run Git state: `HEAD == main == origin/main == fec018fb2d995007d192de7060ab0640b133b4a9`; worktree clean.
- Scientific implementation: RR5d commit `caa8572556a8d892ae07936fc4e6460c73819285`.
- Manuscript: commit `47af37560d6e0752519d6a6267a152871d8a5157`, SHA-256 `d1c7224e0962df6a54fb92a550671f72b66258f30afd64d292f3bf60c7d87be9`.
- Calibration: commit `24ae1b5817bedc405c8debf12d6b6f9bf8cb9375`, SHA-256 `51b4ea76fdd2295a2cfeb840ef7b2b3f468d32c45ad2c123523363d5b968f17e`.
- Pre-run validation: 182 tests passed; Ruff passed; `git diff --check` passed.
- Master seed: `2026081102`. Exact working-tree and Git-history searches found no pre-launch occurrence or prior scientific use.
- Authoritative engine configuration hash: `f6c655e22d7c8930`; fingerprint hash: `1da871755ec1acaf314fddb18c8a9bdba1074e4ecdde1f3cbfdb5752f80f0009`.
- Design: 4 DGPs × 3 truths × 5 panel cells × 3 replications = 180 semantic tasks.
- Execution: requested/effective `n_jobs=8`, one BLAS/OpenMP thread per worker, no nested pool.
- Scope: rank-only. The inference diagnostics dataset has zero rows; no Riesz equation, split fit, standard error, interval, coverage, or rejection calculation ran.

The prelaunch fingerprint also records the dry-run configuration hash `9b7c3db446b5383f`; the engine's authoritative hash differs because it is computed after frozen group-gap resolution (`resolved_group_gap=0.2`). Both identify the same immutable input file, whose SHA-256 is `2da353b511e487a501d3d31d4c03726d0792c788fb65dba22e80437dfdc93023`.

## First question: numerical acceptance

| Quantity | RR5e | Frozen RR5 comparator |
|---|---:|---:|
| Attempted pilots | 180 | 180 |
| Accepted pilots | 0 | 0 |
| Unresolved pilots | 180 | 180 |
| Total starts | 540 | 540 |
| Individually valid starts | 35 | 53 |
| Interior-path successes | 13 | not separately frozen |
| Constrained-path successes | 22 | not separately frozen |
| Stationarity/KKT failures | 505 | 487 non-passing starts |
| Boundary-active starts | 173 | 180 |

Start-path accounting is complete: 352 starts used the interior fast path, of which 13 passed; 188 invoked constrained fallback, of which 22 returned `success` and 166 returned `constrained_optimality_failure`. Maximum constraint violation was numerical noise (`8.17e-14`). Median, p95, and maximum stationarity/KKT residuals were `1.1647e-5`, `0.0031914`, and `0.0126333`.

Pilot valid-start counts were: 150 with zero, 25 with one, and 5 with two valid starts. None of the five valid pairs met the normalized objective-agreement tolerance `1e-6`; consequently all 180 pilots failed objective stability. Across all three starts, the median and p95 objective spreads were `0.00179179` and `0.0142077`.

The comparison answers the first question adversely. RR5d's resumability and performance engineering worked, but the improved solver did not fix universal pilot rejection under the frozen gate.

## Second question: rank selection

The selector was not reached. Exact, under-only, over-only, and mixed counts are all zero because these categories are conditional on an accepted pilot; they are not statistical failures.

### By truth

| Truth | Attempted | Accepted | Unresolved | Exact | Under | Over | Mixed | Valid post-refit |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `(1,0,2)` | 60 | 0 | 60 | 0 | 0 | 0 | 0 | 0 |
| `(1,1,1)` | 60 | 0 | 60 | 0 | 0 | 0 | 0 | 0 |
| `(2,1,1)` | 60 | 0 | 60 | 0 | 0 | 0 | 0 | 0 |

### By DGP

| DGP | Attempted | Accepted | Unresolved | Exact | Under | Over | Mixed | Valid post-refit |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 45 | 0 | 45 | 0 | 0 | 0 | 0 | 0 |
| 2 | 45 | 0 | 45 | 0 | 0 | 0 | 0 | 0 |
| 3 | 45 | 0 | 45 | 0 | 0 | 0 | 0 | 0 |
| 4 | 45 | 0 | 45 | 0 | 0 | 0 | 0 | 0 |

### By panel dimensions

| N×T | Attempted | Accepted | Unresolved | Exact/under/over/mixed | Median pilot runtime (s) |
|---|---:|---:|---:|---:|---:|
| 50×50 | 36 | 0 | 36 | 0/0/0/0 | 43.375 |
| 100×100 | 36 | 0 | 36 | 0/0/0/0 | 48.517 |
| 100×200 | 36 | 0 | 36 | 0/0/0/0 | 52.199 |
| 200×100 | 36 | 0 | 36 | 0/0/0/0 | 52.640 |
| 200×200 | 36 | 0 | 36 | 0/0/0/0 | 68.913 |

## Block, cap, post-refit, and rectangular audits

- Zero-rank B: 60 truth-`(1,0,2)` tasks, all marked `not_reached_due_to_pilot_failure`; selected-B frequencies 0/1/2/3 are all zero.
- Rank-two A: 60 truth-`(2,1,1)` tasks, all not reached; selected-A frequencies 0/1/2/3 are all zero.
- Rank-two H: 60 truth-`(1,0,2)` tasks, all not reached; selected-H frequencies 0/1/2/3 are all zero.
- Cap selections: zero observed because no rank vector was selected.
- Final post-refits: zero attempted, zero valid, and 180 not reached before post-refit.
- Rectangular cells: both 100×200 and 200×100 had 36 attempts, zero accepted pilots, 36 unresolved tasks, and no selector-dependent block or post-refit result. Their median runtimes were 52.199 and 52.640 seconds. No aspect-ratio inference is drawn.

The block audit files retain one explicit not-reached row per relevant task and use missing values for spectra, ratios, ratio gaps, and selected ranks. No numerical quantity was fabricated.

## Performance and resumability

- Total execution wall time: 1,582.012 seconds (26.37 minutes).
- Task time median/p95/max: 54.051 / 181.211 / 299.094 seconds.
- Summed pilot time: 11,031.644 seconds; 89.46% of summed task wall time.
- Summed final-refit time: zero because no pilot was accepted.
- CPU utilization: 37.39% of 20-logical-CPU capacity; worker utilization 97.44%.
- Peak worker RSS: 132,042,752 bytes (125.93 MiB).
- Serialization/write overhead: 0.405 / 0.762 seconds.
- Process-spawn setup and idle tail: 26.862 / 50.322 seconds.
- Maximum pending queue depth recorded: 15.
- Interruptions/resume events/skipped-on-resume/corrupt bundles: 0/0/0/0.
- All 180 semantic task IDs and 180 realization hashes are unique. Fingerprint validation passed.

Compared with RR5's approximately 2.760-hour wall time and approximately 98.3% pilot share, RR5e was operationally faster and had a lower pilot share. Runtime is computational evidence only, not statistical evidence.

## Frozen-scope confirmation and decision

No tuning, fallback, additional replication, inference, scientific-source edit, manuscript edit, or calibration edit occurred after results appeared. The frozen ridge, caps, starts, tolerances, box, selector, and rank-zero rule were unchanged.

Post-run validation passed: ordinary `pytest -q` reported 182 passed in 24.75 seconds, Ruff reported all checks passed, and `git diff --check` passed. Git showed no tracked change under `src/`, `scripts/`, `tests/`, `configs/`, or the manuscript tree; only RR5e evidence was new.

Decision: **RR5e-NUMERICAL-NO-GO**. The selector was not meaningfully evaluated because pilot acceptance remained 0/180. The correct scientific interpretation is numerical non-evaluation, not a finding that the ridge-ratio rank rule statistically fails.
