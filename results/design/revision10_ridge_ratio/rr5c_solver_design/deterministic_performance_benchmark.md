# Deterministic performance benchmark design

## Workload

Use at least 24 independent tasks built from fixed synthetic matrices and literal seeds dedicated to engineering tests. Do not call scientific DGP generators, use RR5 seeds/calibration cells, or score rank recovery. Tasks exercise the same NumPy/SciPy/threadpool stack and the proposed numerical kernels:

- 8 interior width-`(4,4,4)` joint fits, including balanced and ill-conditioned factor gauges;
- 8 active-box constrained fits with repeated dense row/time QPs;
- 4 lower-rank-boundary/full-width fits;
- 4 fixed selected-rank-style post-refits with widths representative of `(1,1,1)`, `(2,1,1)` and `(1,0,2)` only as computational dimensions, not truths.

Use a balanced mix of `50x50`, `100x100`, `100x200`, `200x100`, and a limited number of `200x200` problems. Calibrate repetition of the *fixed kernel operations* in a separate dry engineering step so the serial suite lasts minutes, not hours; do not tune solver rules or data from rank outcomes. Freeze the literal matrices and workload repetition count before comparing workers.

## Execution matrix

Run requested worker counts `1,4,8,12,14`; 14 is justified by the machine's reported 14 physical cores. Ensure at least 24 tasks so each requested count has a distinct effective count. Perform one unreported warm-up for import/cache stabilization, followed by at least three measured repetitions per worker count in rotated deterministic order to reduce time drift.

For each repetition record all measurements required by `parallel_performance_freeze_spec.md`, plus effective worker count, machine CPU topology, available/total memory, Python/NumPy/SciPy/BLAS versions, scientific commit, benchmark-input hash and output hash.

## Scientific-equivalence check

Before timing interpretation, compare worker count 1 against every other count for semantic task IDs, seeds, status, products, objectives, spectra, diagnostics and canonical ordering. Require exact equality for discrete/serialized fields and only existing accepted tolerances for floating-point fields. A mismatch invalidates that worker count and triggers diagnosis; timing cannot excuse it.

## Reporting

For each count report median and range across repetitions for wall time, tasks/minute, CPU utilization, per-worker/peak memory, spawn/serialization/write overhead and idle tail. Include p50/p95 task and phase time, queue-depth timeline, failures/retries and one-thread verification. Select workers using the frozen 5%-of-best rule and memory/stability gates, never rank outcomes.

No benchmark is run in RR5c; this document freezes the future protocol only.
