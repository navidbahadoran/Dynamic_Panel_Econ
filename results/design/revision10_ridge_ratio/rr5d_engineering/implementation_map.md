# RR5d implementation map

## Cap+1 solver

- `src/dynamic_panel_econ/cap_plus_one.py` adds the isolated Revision-10 cap+1 kernel. Ordinary `fit_fixed_rank` is unchanged.
- `rank_selection.fit_revision10_spectral_pilot` alone calls `fit_cap_plus_one`; the final selected-rank post-refit and every supplied-rank/Revision-9 path continue to call the legacy fixed-rank estimator.
- Interior updates use reversible column equilibration, SVD/QR least squares, balanced product-preserving block factorizations, monotone rollback, condition/stationarity/objective traces and matrix-free Gauss-Newton refinement.
- Constrained updates solve the identical row/time convex quadratic programs by a deterministic active-set algorithm under reversible internal equilibration. The literal product box and existing factor-space KKT evaluator are retained.
- The three maintained start seeds, cap widths, `B`, objective, stationarity/KKT thresholds and objective-stability gate are unchanged.

## Resumability and orchestration

- `src/dynamic_panel_econ/resumability.py` implements canonical task specifications, scientific fingerprints, one atomic JSON bundle per task, payload/file hashes, a sole-writer atomic live manifest, corruption quarantine, terminal failure preservation and resume reconciliation.
- `monte_carlo.run_monte_carlo` integrates the store before execution, checkpoints each returned task before writing legacy chunk outputs, skips validated terminal bundles, requeues only missing/stale work and records interruption events.
- Existing semantic RNG derivation is recorded, not replaced. Canonical legacy chunk aggregation remains sorted by semantic record keys.

## Performance engineering

- Workers still receive immutable resolved configuration/calibration once through the pool initializer, use one native thread and generate task-owned arrays locally.
- Per-task metrics now include wall/CPU time, pilot/final-fit phase time, RSS, native threadpool state, completion time, queue depth, serialization and atomic-write time.
- The run manifest summarizes utilization, memory, spawn/setup, phase time, idle tail and queue depth.
- `scripts/benchmark_rr5d_engineering.py` runs the non-scientific 24-task, five-worker-count benchmark.
- `scripts/validate_rr5d_engineering.py` writes deterministic solver/equivalence evidence without calling any scientific DGP.

## Tests

- `tests/test_cap_plus_one_engineering.py`: rank 0/1/2 inside width four, gauge/product/objective invariance, bad scaling, near-collinearity, active box/KKT, monotone rollback, frozen three-start gate, supplied-rank regression and Windows spawn determinism.
- `tests/test_resumability.py`: uninterrupted versus interrupted/resumed equality, stale-worker recovery, terminal failure preservation, fingerprint/task-set refusal, corruption/incomplete quarantine and manifest states.
- `tests/test_revision10_rank_selection.py` now mocks the isolated pilot solver at the correct boundary; selector formulas and gate assertions are unchanged.
