# Frozen safe performance improvements

Every item below is permitted only after implementation tests show unchanged semantic task IDs, RNG streams, coefficient outputs, statuses and deterministic aggregates.

| Improvement | Exactness guard |
|---|---|
| Initialize and hash immutable config/calibration once per worker | Same canonical objects and selected rows as per-task parsing; no mutable cache |
| Avoid repeated TOML parsing | Parse once before/pool initialization; verify frozen hashes in task bundle |
| Reuse one spawn-safe process pool | No nested pools; task seeds independent of worker/order |
| Dynamic work-conserving scheduling | Canonical task list/output ordering; completion order never enters RNG or summaries |
| Balance dispatch by frozen computational fixture class | Scheduling metadata only; no statistical outcome or adaptive scientific branching |
| Send semantic keys rather than large panels | Generate task-owned inputs using unchanged seed semantics in worker |
| Avoid unnecessary array copies and reuse workspaces | Alias tests prevent mutable cross-task state; numeric kernel equivalence required |
| Cache row/time layout and half-sweep constraint structures | Cache only objects invariant for that scope; invalidate deterministically when factors change |
| Reuse exact QR/SVD/QP factorizations within one unchanged subproblem or half-sweep | Key cache by complete numeric inputs; no approximate reuse after inputs change |
| Warm-start constrained active sets | Same convex QP and tolerances; cold/warm reference equality tests |
| Compact lossless serialization | Round-trip exact arrays and canonical fields; no float downcast/compression loss |
| Write atomic task bundles as tasks finish | One bundle per semantic ID; fsync/rename protocol; canonical later aggregation |
| Read-only shared/memory-mapped immutable data | Byte/hash equality under Windows spawn; never share mutable estimator state |
| Explicit native-thread limit | Verify one BLAS/OpenMP thread inside every worker |
| Resource telemetry and queue tracing | Observational namespace excluded from scientific hashes and decisions |

The preferred QP-specific solver, internal reversible preconditioning, interior product-preserving gauge, and monotone safeguards are semantically exact numerical improvements but must pass the dedicated solver-equivalence tests before they are treated as performance improvements.

Not authorized as engineering: fewer/changed starts, looser thresholds, altered `B`/caps/ridge/anchor, singular-value selection thresholds, approximate ratio spectra, early acceptance, changed DGP/calibration, or skipped difficult tasks.
