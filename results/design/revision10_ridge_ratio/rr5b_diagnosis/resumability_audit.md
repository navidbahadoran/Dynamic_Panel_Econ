# Resumability audit

| Capability | Status | Evidence |
|---|---|---|
| Deterministic semantic task IDs | PRESENT | ID is a function of DGP, N, T, replication and truth. |
| Execution-order-independent seeds | PRESENT | `seed_sequence` uses semantic keys; worker order is absent. |
| One durable record per completed task | PARTIAL | Results remain in parent memory until all three tasks in a chunk finish; durability is five files per three-task chunk. |
| Atomic write | PARTIAL | Each parquet uses `.partial` then `os.replace`; file close is implicit, but no explicit file or directory `fsync`, and the five-file chunk is not transactional. |
| Expected/running/completed/failed/unresolved manifest | ABSENT | Final manifest is written only after task execution; there is no live state ledger. |
| Restart detection | PARTIAL | `--resume` checks existence of all five chunk files only. It does not validate schema, hashes or row identities. |
| Exact config verification | PARTIAL | Output directory contains a resolved-config hash, but resume rewrites config files and does not compare an existing manifest. |
| Exact code and calibration verification | ABSENT | Commit and calibration hash are recorded only in the final manifest and are not checked before skipping/mixing chunks. |
| Skip completed | PARTIAL | Complete five-file chunks are skipped; individual completed tasks in an incomplete chunk are rerun. |
| Protection from mixed commits/config/calibration | PARTIAL | Config-hash directory helps, but code changes and calibration-content changes can share it. |
| Failed/unresolved preservation | PRESENT | Failure and replication rows are checkpointed like successes. |
| Deterministic aggregation | PRESENT | Chunk names and rows are sorted; aggregation is independent of completion order. |
| Automatic interruption log | ABSENT | RR5's log was manually created after the wrapper timeout. |
| Resume after user interruption | PARTIAL | Completed chunks survive; in-flight/incomplete chunks do not. |
| Resume after reboot | PARTIAL | Replaced chunk files usually survive, but explicit durability and transaction validation are absent. |
| Resume after worker crash | PARTIAL | Completed chunks survive; `future.result()` aborts the run and no task-state recovery occurs. |

RR5 empirically demonstrated useful chunk-level resume: 117 attempts in 39 complete chunks were skipped and 63 remaining attempts completed. That success does not make the driver task-atomic or fingerprint-safe.
