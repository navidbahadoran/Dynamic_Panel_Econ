# Frozen resumability specification

## Equality contract

For a fixed scientific fingerprint and task collection,

```text
uninterrupted run == interrupted + resumed run
```

for semantic task identities, seeds, scientific arrays/records, terminal statuses, diagnostics and canonical summaries. Checkpointing may not transform a numeric output or consume random numbers.

## Semantic identity and scientific fingerprint

The canonical task specification contains, in fixed key order:

- master seed and the existing execution-order-independent task seed key/state;
- DGP identifier, `N`, `T`, true-rank vector and replication index;
- selector method and task type;
- scientific code commit;
- canonical resolved-config SHA-256;
- complete frozen-calibration SHA-256 and selected-row SHA-256;
- output schema version.

Serialize this specification as UTF-8 canonical JSON (sorted keys, fixed separators, no NaN) and define `semantic_task_id = sha256(canonical_json)`. The run fingerprint hashes the immutable collection-level fields and the canonically sorted expected task IDs. Existing RNG derivation is retained exactly; the semantic ID records and verifies it but does not replace it.

Launch is permitted only from a clean scientific commit. A resume requires exact equality of every fingerprint field. A mismatch is a hard stop, never an invitation to mix results.

## One atomic bundle per task

Each task produces one immutable self-contained bundle at `tasks/<first-two-hash>/<semantic_task_id>.bundle`. The bundle contains:

1. canonical task specification and fingerprint;
2. complete seed metadata/state needed for audit;
3. scientific outputs and coefficient arrays;
4. spectra, ratios, ranks, status and diagnostics;
5. schema/member hashes and terminal classification;
6. no volatile timestamps inside the scientific payload.

The bundle is a deterministic single-file container with fixed member order, fixed metadata and lossless binary arrays. Write protocol:

1. write a uniquely named temporary file in the destination directory;
2. flush every stream and `fsync` the file;
3. close and reopen it; validate schema, member hashes, semantic ID and bundle hash;
4. atomically rename/replace it to the final task path;
5. `fsync` the directory when supported;
6. only then atomically update the live manifest to a terminal state.

The parent/coordinator is the sole manifest writer. If workers write bundles directly, each writes only its unique task path and returns the verified hash; no two workers may own the same task concurrently.

## Durable live manifest

The manifest contains the immutable run fingerprint, expected task IDs, aggregate version, interruption log and for each task exactly one state:

- `expected`: eligible, never durably started or reset after stale-running recovery;
- `running`: attempt ID, worker PID, start time and heartbeat; nonterminal;
- `completed`: valid successful scientific bundle hash;
- `failed`: valid terminal computational/scientific failure bundle hash;
- `unresolved`: valid terminal numerical-unresolved bundle hash.

`failed` and `unresolved` are preserved terminal results and are skipped on resume just like `completed`. Every manifest mutation is written to a new temporary manifest, flushed, fsynced, validated and atomically replaced. An append-only interruption/event log is stored in the manifest snapshots or a separately checksummed atomic journal; it records Ctrl-C, worker death, parent exception, stale-running recovery, quarantines and resume time without entering scientific equality comparisons.

## Resume algorithm

1. Acquire an exclusive run lock; refuse two coordinators.
2. Verify code cleanliness and exact run fingerprint.
3. Scan all final task bundles. Validate container integrity, schema, semantic ID, fingerprint, bundle/member hashes and one-to-one expected membership.
4. Quarantine incomplete, corrupt, duplicate or mismatched files by atomic rename to `quarantine/`; record the reason. Never silently overwrite or count them.
5. Reconcile manifest with bundles. A valid terminal bundle is authoritative even if a crash occurred before the manifest update. A `running` task with no valid bundle becomes `expected` after confirming its worker/attempt is stale.
6. Skip every valid `completed`, `failed`, or `unresolved` task without invoking its RNG or estimator.
7. Queue only `expected` tasks, in canonical semantic-task order; dynamic completion order is allowed.
8. After all tasks are terminal, aggregate validated bundles in canonical semantic-task order. Write each aggregate by the same temp/fsync/atomic-rename protocol and record its hash.

## Crash guarantees

- **Ctrl-C:** coordinator stops submission, records interruption, allows a short bounded drain, atomically persists received bundles/states, then terminates workers.
- **Worker crash:** mark attempt abandoned, validate whether a bundle exists, and requeue only if none is valid.
- **Python parent crash:** on restart, recover valid bundles written after the last manifest snapshot and reset stale `running` tasks.
- **Machine reboot/power loss:** fsynced final bundles survive; incomplete temp files are ignored/quarantined; reconciliation reconstructs state.

No guarantee depends on a three-task chunk completing. Durability granularity is one semantic task.

## Deterministic aggregation

Never aggregate in completion order. Sort by the canonical semantic task specification and validate exactly one terminal bundle per expected ID. Summary code consumes those ordered records and retains terminal failures. Volatile runtime telemetry is reported separately or in explicitly excluded fields so it cannot make scientific outputs differ.

This specification is frozen for every future long run and must be implemented and pass `resume_equivalence_test.md` before any larger Monte Carlo is authorized.
