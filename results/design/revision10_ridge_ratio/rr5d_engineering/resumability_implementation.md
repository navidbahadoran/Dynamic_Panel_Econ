# Resumability implementation

## Identity and fingerprint

Each task specification contains DGP, `N`, `T`, true-rank vector, replication, master seed and selector. Its canonical SHA-256 is the task ID. The scientific fingerprint contains Git commit, current source-tree hash, resolved config hash, calibration hash, selector, master seed, relevant DGP/estimation/inference configuration and bundle schema. Resume refuses any fingerprint or expected-task-set mismatch before computation.

## Atomic bundle protocol

One compact lossless JSON bundle is written per semantic task. NumPy values are converted without precision reduction; nonfinite floats use explicit tagged encodings. The payload carries its own canonical SHA-256. The write uses a unique same-directory temporary file, flush, file `fsync`, close, atomic `os.replace`, then validation and manifest update. On platforms supporting directory descriptors, the directory is also fsynced.

Temporary, corrupt, wrong-schema, wrong-task, wrong-fingerprint or hash-mismatched bundles are atomically moved to `task_quarantine` and recorded. They are never counted as terminal work or silently deleted.

## Manifest and restart

The sole coordinator writes `task_manifest.json` by the same atomic protocol. It exposes expected IDs plus pending, running, completed, failed, unresolved and corrupt lists. Terminal `failed` and `unresolved` bundles are durable and skipped. On restart, a valid bundle is authoritative even if the preceding process died before updating the manifest; a stale running task without a bundle returns to pending.

Ctrl-C and parent/worker exceptions are logged before re-raising. A worker-process death leaves already returned bundles durable and missing/running work recoverable. No task seed is derived from worker or queue order.

## Aggregation

Validated bundle rows are placed into the existing deterministic chunk writer, which sorts by semantic replication and record keys. Completion order is irrelevant. Existing parquet/CSV aggregate layouts remain compatible.
