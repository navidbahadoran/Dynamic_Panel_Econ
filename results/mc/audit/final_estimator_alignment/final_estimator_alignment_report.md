# Final estimator-alignment report

## Decision

The approved Revision-8 DGP calibration and literal entrywise box-constrained fixed-rank estimator are activated. This work performed deterministic checks and small estimator integration tests only. No medium diagnostic, rank-stress production, power, N=400, or 1,000-replication job was started.

## Frozen calibration activation

The maintained official Monte Carlo configurations use `configs/mc/frozen_dgp_calibration.toml`. Production configuration without a frozen table is rejected. Startup validation verifies the requested cell, analytical `c_H`, intended pooled-R2 target, and stored deterministic coefficient envelopes against the active DGP parameters.

The activated table contains 52 cells and has SHA-256:

```text
e8983cadc4fbca990feeba6363420542a99ee056cf445e867532e2a6ea0e7d62
```

`activated_calibration_summary.csv` records cell-level `c_H`, `c_xi`, intended and achieved R2, `pi_H`, and coefficient envelopes.

## Coefficient envelope

The active common coefficient bound is `B = 10` with simulation margin `c_B = 1`. Across all 52 approved cells, the maximum deterministic true-coefficient envelope is:

```text
max C_Theta = 8.410761115894578
B - max C_Theta = 1.5892388841054217
(B - c_B) - max C_Theta = 0.5892388841054217
```

Every cell satisfies `C_Theta <= B - c_B`. The cell-level verification is in `coefficient_envelope_verification.csv`.

## Estimator algorithm

All maintained fixed-rank routes use the shared `estimation.fit_fixed_rank` implementation: supplied-rank full-panel fits, retained candidate post-refits, rank-at-most-cap pilot refits, and all four split fits.

The algorithm first runs the existing unconstrained alternating least-squares fit. When its reconstructed matrices are strictly inside the box by the numerical interior tolerance, that solution is retained as the fast path. Otherwise, the fallback starts from a common deterministic rescaling of the fitted matrices into the box and alternates:

1. row-factor convex quadratic subproblems with exact linear inequalities on reconstructed entries;
2. time-factor convex quadratic subproblems with the same exact entrywise inequalities;
3. joint objective and feasibility checks after each sweep.

QR factor renormalization preserves the reconstructed matrices and ranks. No entrywise clipping is used. A factor-coordinate active-set KKT residual is computed using nonnegative least squares on active box normals. Successful boundary solutions are valid constrained estimates; infeasibility, solver failure, or excessive KKT residual receive explicit failure statuses.

## Deterministic validation

`constrained_solver_tests.csv` records:

- an interior case retained by the fast path;
- an outside-box unconstrained case routed to the constrained fallback;
- exact box feasibility within tolerance;
- active-boundary KKT diagnostics and rank preservation;
- an explicit forced solver-failure test reference;
- successful evaluation of the existing inference computation at a valid boundary estimate.

`estimator_path_alignment.csv` records the shared optimizer route for each maintained estimator path. The separate boundary audit explains why the implementation does not use the interior normal equation as a zero identity and why theorem interpretation nevertheless remains tied to interiority.

## Run authorization

The implementation is ready for bounded validation and independent review. Starting any paused medium or production experiment remains outside this audit and requires a separate explicit authorization after review of these artifacts.

## Validation result

The final local validation completed successfully:

```text
pytest -q: 121 passed, 1 non-test-failure cache warning
ruff check src scripts tests: All checks passed
git diff --check: passed
```
