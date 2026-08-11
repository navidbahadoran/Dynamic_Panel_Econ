# RR5 decision

## Decision: RR5-NO-GO

The fresh Revision-10 rank-selection preflight is a NO-GO. All 180 prespecified attempts ended as `rank_selection_numerically_unresolved`; none produced a numerically acceptable cap+1 spectral pilot. Consequently, the frozen ridge-ratio selector, selected-rank distribution, cap audit, and final literal post-refit were not reached. No rank result is inferred from a failed pilot.

This decision uses the frozen implementation at `24ae1b5817bedc405c8debf12d6b6f9bf8cb9375`, master seed `2026081101`, the frozen calibration table, and the exact 180-attempt design. It does not use Revision-9 candidate/IC evidence and it does not authorize tuning or RR5 expansion.

Numerically, 53 of 540 deterministic pilot starts passed the frozen stationarity threshold, while 487 did not. Fifty attempts had at least one stationarity-passing start and three had at least two, but zero attempts passed the complete maintained cap+1 pilot acceptance diagnostics. The stationarity residual median was `1.32622364e-05` (maximum `0.02679843`) against the frozen `1e-6` tolerance.

The initial process was stopped by an external two-hour execution-wrapper timeout after 117 attempts had been atomically checkpointed. The repository's resume-safe continuation used the identical configuration, seed, and 12-worker setting, skipped those completed chunks, and generated only the remaining 63 attempts. Accounting verifies 180 unique semantic IDs and zero checkpointed replication reruns.
