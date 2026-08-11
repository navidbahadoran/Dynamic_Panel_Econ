# Frozen performance and parallelization audit

## Mandatory measurements

Before every larger authorized run, execute the non-scientific benchmark in `deterministic_performance_benchmark.md`. For every worker count and task record:

- total wall time and completed tasks/minute;
- per-task elapsed and process CPU time;
- pilot, final-fit and other phase times;
- coordinator and aggregate worker CPU utilization sampled through the run;
- worker RSS distribution, memory per worker and system peak memory;
- Windows process-spawn/pool-initialization overhead;
- parent-to-worker serialization/deserialization time;
- atomic bundle write/fsync time and aggregate write time;
- idle-tail duration from the last time all workers were busy until final task completion;
- ready/running/completed queue depth over time;
- worker PID/task timeline and abnormal exits.

Definitions, sampling interval, logical/physical CPU denominator and memory units must be recorded in a machine-readable manifest. CPU utilization must not be replaced by `sum(task time)/(jobs*wall)`; that ratio may be reported only as an accounting proxy.

## Concurrency invariants

1. Set and verify one BLAS/OpenMP thread per worker both by environment and `threadpoolctl` inspection inside every worker.
2. Use one Windows-spawn-safe outer process pool and prohibit nested process pools.
3. Derive seeds from semantic task identity, never queue position or worker identity.
4. Submit/aggregate in deterministic semantic order while permitting dynamic work-conserving dispatch.
5. Initialize immutable config/calibration/fingerprint objects once per worker, verify their hashes, and avoid reparsing per task.
6. Do not repeatedly pickle large immutable arrays. Generate task-owned arrays in the worker or use verified read-only shared storage when scientifically identical.
7. Bound in-flight futures and total memory; record the bound.
8. Keep runtime telemetry observational; it may not control scientific branches.

## Worker-count selection rule

Eligible counts are `1,4,8,12` and optionally the machine's 14 physical cores. Reject a count if it violates determinism, one-thread compliance, memory headroom, worker stability, or resumability. Among remaining counts, choose the smallest count within 5% of the best median throughput across at least three repeated benchmark runs; break ties by lower peak memory, then shorter idle tail. Do not use statistical outcomes.

The selected count and full benchmark evidence are frozen into the next run manifest. A machine, dependency stack, solver architecture or representative-dimension change requires a new benchmark.

## Larger-run gate

No larger run may start unless:

- the solver correctness suite passes;
- interrupted/resumed equivalence passes;
- performance benchmark completes without unresolved engineering failures;
- estimated peak memory plus 25% safety margin is below available memory;
- projected duration and checkpoint volume are recorded;
- Git commit/config/calibration fingerprints are frozen.
