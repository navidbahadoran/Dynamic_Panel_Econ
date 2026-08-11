# Safe performance and engineering improvements

## SAFE / SEMANTICALLY EXACT

| Improvement | Exactness requirement |
|---|---|
| Per-task atomic checkpoint and fingerprinted manifest | Store existing records without numeric transformation; preserve seed keys. |
| Dynamic scheduling | Keep semantic task list and output ordering canonical; completion order must not affect seeds or aggregates. |
| Deterministic batching, if useful | Batch transport only; retain one semantic seed, fingerprint and atomic terminal record per task. |
| Process-pool reuse | One spawn-safe outer pool, no nested pools. |
| Cache immutable config/calibration per worker | Verify the same canonical hashes; never cache mutable scientific state. |
| Precompute immutable index/layout objects per `(N,T,rank caps)` in each worker | Confirm byte/numeric equivalence of constructed arrays. |
| Avoid unnecessary array copies and preallocate workspaces | Require exact/tolerant kernel equivalence tests. |
| Read-only shared or memory-mapped immutable inputs | Use only where Windows spawn semantics and byte identity are verified; never share mutable estimator state. |
| More compact task serialization | Send semantic keys only; generate panels in the owning worker as now. |
| Deterministic final aggregation | Sort by semantic key and validate one terminal record per task. |
| Explicit one-thread native-library control | Retain `threadpool_limits` and environment settings. |
| Resource telemetry | Record worker PID, CPU time, peak RSS and task timestamps without affecting computation. |

## REQUIRES SCIENTIFIC DECISION

Fewer starts, changed stationarity/KKT tolerance, altered objective-stability rule, different optimizer or stopping rule, different caps or ridge, approximate SVD that changes spectra, and skipping difficult cells all change the frozen numerical/statistical procedure. None is recommended or implemented in RR5b.
