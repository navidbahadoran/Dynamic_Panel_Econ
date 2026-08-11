# Resume-equivalence results

The mandatory deterministic non-scientific collection passed.

- Uninterrupted run: six semantic tasks, including completed, failed and unresolved terminal states.
- Interrupted run: first three terminal bundles written; a Ctrl-C event recorded; a fourth task left in simulated stale-running worker state.
- Resume: three terminal tasks skipped, stale running task requeued, remaining tasks completed.
- Comparison: exact task specifications/IDs, seeds, terminal statuses, discrete ranks, objectives and canonical task ordering matched the uninterrupted run.
- Parent restart was represented by constructing a fresh store from disk.
- Worker failure/restart was represented by stale-running reconciliation with no bundle.
- Terminal failure and numerical-unresolved states were preserved rather than recomputed.
- Fingerprint and expected-task-set mismatches were rejected.
- Corrupt final and incomplete temporary bundles were detected, quarantined and recorded.

The tests compare scientific payloads only. Volatile PIDs, timestamps and resource measurements live in the metrics namespace and are intentionally outside scientific equality.
