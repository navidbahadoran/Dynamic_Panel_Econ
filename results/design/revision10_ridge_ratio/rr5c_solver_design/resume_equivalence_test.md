# Mandatory interrupted/resumed equivalence test

## Collection

Create at least 24 deterministic non-scientific tasks spanning interior, active-box, lower-rank, deliberately ill-conditioned, terminal failure and numerical-unresolved fixtures from the solver test plan. Use synthetic fixed arrays, not DGP1-DGP4, RR5 seeds, calibrations or statistical rank-recovery targets.

Run the collection under one immutable fingerprint in two separate output roots:

- **A — uninterrupted:** complete all tasks once.
- **B — interrupted/resumed:** start the identical collection, deliberately interrupt after a deterministic set of task bundles is durable, then resume to terminal completion.

Repeat B for four failure modes: coordinated Ctrl-C, forced worker death, abrupt parent-process death, and reboot-equivalent restart after all processes are gone. A real reboot is optional in CI; the test must at least close all handles/processes and restart from disk in a fresh Python process.

## Required comparisons

After canonical aggregation, compare A and every B run:

- expected and terminal semantic task IDs;
- recorded task seed keys/states;
- realization/input hashes;
- status codes (`completed`, `failed`, `unresolved`);
- coefficient outputs;
- singular spectra and normalized spectra;
- ridge ratios and selected ranks where the fixture legitimately reaches them;
- objectives, stationarity/KKT, feasibility and multistart diagnostics;
- canonical aggregate rows and summaries.

Discrete values, canonical JSON, IDs, hashes, seeds and serialized status fields require exact equality. Floating-point arrays/fields require bitwise equality when serialization and platform are identical; otherwise use only the already accepted production numerical tolerance, documented field by field. No new tolerance may be selected from observed results.

Volatile attempt IDs, PIDs, wall-clock timestamps and resource telemetry are expected to differ and must live outside the scientific comparison namespace.

## State-machine assertions

- No valid terminal task is recomputed; instrument estimator/RNG entry counts.
- A stale `running` task with no bundle is requeued exactly once.
- A bundle committed before a manifest crash is recovered and skipped.
- Truncated temp and final files are detected; corrupt finals are quarantined, not counted.
- Fingerprint mismatch prevents resume before task execution.
- Terminal `failed` and `unresolved` bundles are preserved and skipped.
- Aggregation order is identical despite deliberately reversed completion order.
- Interruption log records every induced event without changing scientific hashes.

## Pass condition

Every scientific comparison and state-machine assertion must pass for every failure mode. Any silent recomputation, mixed fingerprint, lost terminal failure, changed seed stream, or completion-order-dependent summary is a failure. This gate is mandatory before any larger Monte Carlo.
