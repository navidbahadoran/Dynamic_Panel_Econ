# Minimum exact resumability design

## Identity and fingerprint

Define a canonical semantic task key `(master_seed, DGP, N, T, true_rank, replication_index)`. Hash a canonical task specification containing that key plus scientific commit, resolved configuration hash, frozen calibration file hash, selected calibration-row hash, selector method and output schema version. Random-number keys remain exactly unchanged.

## Durable task protocol

1. At launch, create a run manifest atomically with the immutable run fingerprint and every expected task key in `expected` state.
2. Workers compute one task and return only its existing scientific records.
3. The parent validates semantic IDs, seed metadata, schema, finiteness/accounting and task fingerprint.
4. Write one self-contained task bundle to a unique temporary file, flush and `fsync` it, close it, atomically rename it to the content-addressed final task path, then `fsync` the directory where supported.
5. Atomically update the manifest state to `completed`, `failed`, or `unresolved`, preserving all three as terminal completed computations. `running` includes attempt time and worker PID but is never treated as scientific completion.
6. On resume, require exact run fingerprint equality before any skip. Validate every terminal task file's hash, schema and semantic key. Skip valid terminal tasks without regeneration. Quarantine incomplete/corrupt/mismatched files and require an explicit recovery action; never silently mix them.
7. Requeue only expected tasks lacking a valid terminal bundle. A stale `running` state after interruption becomes `expected` only after confirming no valid task bundle exists.
8. Aggregate validated bundles in canonical semantic-key order. Write aggregates atomically and record their hashes. Aggregation is reproducible and can be rerun without scientific computation.

## Equality contract

An uninterrupted run and interrupted-plus-resumed run must have identical task keys, seed states, realization hashes, selected ranks, statuses, stored diagnostics and canonical aggregates. Numeric comparisons use only already-approved tolerances; checkpointing itself performs no numeric transformation.

The minimum implementation unit is per task, not the current three-task/five-file chunk. A compatibility layer can still emit the existing aggregate/chunk layout after all task bundles validate.
