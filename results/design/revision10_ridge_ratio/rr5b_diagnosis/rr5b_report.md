# RR5b diagnosis and simulation-engineering audit

## Finding

RR5 failed before rank selection. The exact universal mechanism was `objective_stability_pass=False` in 180/180 pilots: 177 had fewer than two individually valid starts, and three had two valid starts but objective gaps above `1e-6`. Start-level termination and stationarity failures coincide in 487/540 starts. All 540 objectives were finite, all were box feasible, every numerical rank was `[4,4,4]`, and no exact-rank condition was imposed.

The insufficient-valid/gap-failure split was DGP1 `45/0`, DGP2 `45/0`, DGP3 `42/3`, and DGP4 `45/0`; by truth it was `(1,0,2)` `60/0`, `(1,1,1)` `57/3`, and `(2,1,1)` `60/0`. By panel it was `50x50` `36/0`, `100x100` `35/1`, `200x200` `36/0`, `100x200` `36/0`, and `200x100` `34/2`. Every group therefore had 100% objective-agreement failure. The complete counts, rates, and 60 chunk-level summaries are in `rr5_failure_summary.csv`; worker-specific attribution is unavailable because RR5 did not record worker IDs.

## Stationarity/KKT

Across the heterogeneous stored residual column, min/p05/median/p95/max were `2.03501e-07`, `1.56966e-06`, `1.32622e-05`, `0.00356982`, `0.0267984`. Relative to each path's correct threshold, the corresponding ratios were `0.204`, `0.784`, `6.68`, `42`, `268`.

The interior projected-gradient path used `1e-6`: 14/339 passed, median residual/threshold `6.55`. The constrained factor-space KKT path used `1e-4`: 39/201 passed, median ratio `7.55`. Residuals range from narrowly above to orders of magnitude above threshold; this is not a single near-threshold artifact.

Using each start's correct threshold, 53 passed; 42 failures were between 1x and 2x threshold, 248 were between 2x and 10x, 192 were between 10x and 100x, and 5 exceeded 100x. Thus most failures are modest-to-large rather than orders of magnitude, while a small tail is extreme.

Historical Revision-9 rank-one full fixed-rank diagnostics had median residual `6.47e-7` and 58/75 passing starts. This descriptive comparison uses different designs and does not justify changing tolerances; it only shows that the width-four joint pilot is numerically harder.

## Multistart objectives

All starts had finite objectives. The relative all-start spread min/p05/median/p95/max was `7.31017e-05`, `0.000247209`, `0.00141641`, `0.0117849`, `0.0464768`. Even the minimum spread was about `73.1` times the frozen tolerance. Best starts were distributed across starts 1/2/3, so there is no single consistently dominant initialization and no repeated all-start solution under the maintained rule.

## Feasibility, termination, spectra

Every final start was finite and box feasible; maximum violation was `3.55271e-15` and maximum envelope ratio `1.0000000000000004`. Boundary activity occurred in 180/540 starts and in 90/180 replications, but boundary activity is not a rejection rule. Constrained fallback was used 201 times: 39 successes and 162 `constrained_optimality_failure`. Their unconstrained maximum coefficients had median 15.79, p95 90.87, and maximum 447.84; constrained subproblem iterations had median 514,983 and maximum 1,865,590. The interior path had 14 converged and 325 nonconverged starts. No solver encountered a nonfinite or final-feasibility failure.

Finite coefficient matrices and collapsed SVD diagnostics existed for all rejected starts. Full blockwise spectra were not serialized, and normalization occurs only after pilot acceptance, so no normalized spectra or ridge ratios exist. They were not reconstructed and no ranks were computed.

## Runtime and engineering

Wall time was 2.760 hours including interruption/resume. Task runtime min/median/p95/max was `28.1/485.3/2070.2/3471.8` seconds. Pilot fits accounted for 98.3% of summed task runtime. The tail is heavy, but exact worker-idle tail cannot be reconstructed because task timestamps and worker mappings were not stored.

Current resumability is useful but partial: deterministic IDs/seeds, terminal failure preservation, chunk skipping and deterministic aggregation are present; task-atomic durability, fingerprint checks, a live state manifest, corruption validation, explicit fsync and automatic interruption logging are absent or partial. The proposed design adds fingerprinted per-task atomic bundles and canonical aggregation without changing random-number semantics.

Prior deterministic benchmarks provisionally support `n_jobs=12`; RR5 itself does not contain enough CPU/memory telemetry to validate 12 robustly for cap+1. A future non-DGP engineering benchmark is required before any larger authorized run.

## Classification and scope

**RR5b-NUMERICAL-SOLVER.** No implementation bug or acceptance-scale inconsistency is established. No new Monte Carlo draw, DGP generation, estimator call, rank computation, tuning, source edit, calibration edit, or manuscript edit occurred in RR5b.

## Validation

- `pytest -q`: 167 passed in 17.23 seconds.
- `ruff check src scripts tests`: passed.
- `git diff --check`: passed.
- Tracked source/manuscript/calibration changes: none; only `rr5b_diagnosis/` is new.
