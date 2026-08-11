# RR5d final engineering report

## Outcome

**RR5d-PASS.** The isolated cap+1 solver, task-atomic resume engine and performance instrumentation are implemented and validated on deterministic non-scientific fixtures.

## Solver

The new path combines balanced product-preserving SVD gauges, reversibly equilibrated SVD/QR least squares, monotone rollback and matrix-free Gauss-Newton refinement for interior fits. Literal box fallback uses deterministic active-set QPs and the frozen factor-KKT residual. Rank 0/1/2 products fit within width four, active box constraints pass, and the real three-start fixture passes the unchanged `1e-6` objective gate with all three starts valid.

The solver is isolated to the Revision-10 pilot. Ordinary supplied-rank fitting, split fitting, final selected-rank post-refitting and Revision-9 IC selection retain `estimation.fit_fixed_rank` unchanged. Deterministic old/new coefficient products and objectives agree where the legacy solver already converges.

## Resumability

Each semantic task has a fingerprinted, checksummed, fsynced and atomically renamed bundle. The atomic live manifest records pending/running/completed/failed/unresolved/corrupt states. Resume validates fingerprints and bundles, quarantines corruption, skips all terminal work, preserves failures and aggregates in semantic order. Uninterrupted and interrupted/resumed scientific payloads match exactly in the deterministic suite.

## Performance

The 24-task benchmark completed three repetitions at 1/4/8/12/14 workers with exact cross-worker hashes, one native thread and zero worker failures. Median walls were `20.330/3.077/2.745/3.015/3.284` seconds. Eight workers is the sole count within 5% of best and is recommended. Four deliberately short-budget task pilots were nonconverged identically at every count; dedicated acceptance fixtures converge and this invariant was not used to tune the solver.

## Scope

No scientific Monte Carlo, DGP generation, RR5 seed, rank-recovery evaluation, inference, calibration/config/manuscript change or outcome-driven tuning occurred. The exact file list, validation results, commit and push status are reported at handoff.

Final validation: Pytest `182 passed in 42.08s`; Ruff passed; `git diff --check` passed.
