# Memory and serialization audit

- Scientific arrays remain float64; no downcast, lossy compression or clipping is used.
- Task bundles use compact canonical JSON with exact Python float round trips and tagged nonfinite values. Payload and complete-file SHA-256 values detect corruption.
- Benchmark serialization time was below one millisecond in median aggregate per 24-task leg; atomic output writes were approximately 0.08–0.10 seconds per leg.
- Peak worker RSS was approximately 122–124 MB across counts. Summed per-worker peaks rose from 123.6 MB at one worker to 1.70 GB at 14, well below installed memory but with no throughput benefit beyond eight workers.
- Resolved configuration and calibration objects are pickled once per spawned worker, not per task. Panels remain worker-owned. Only result records and compact metrics return to the coordinator.
- Reusable half-sweep arrays and constraint layouts remain local to a fit. Mutable solver state is never cached across semantic replications.
- Atomic per-task writes eliminate the former requirement to hold an entire three-task chunk before durable output. Legacy chunks are reconstructed from validated task rows for compatibility.
