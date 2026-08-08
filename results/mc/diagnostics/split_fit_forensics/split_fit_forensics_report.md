# Split-fit forensics

This diagnostic preserves the accepted estimator, B=9, split formula, Riesz construction, IC rate, DGPs, and all numerical tolerances. No observations were trimmed or removed.

## Existing N=50 preflight

The table contains 56 actual fits from the four-fit bundles of every replication with a required split failure; 27 fits are invalid. Every invalid fit is classified as `coefficient_bound_hit`.
Fixed rank has 15 invalid fits (DGP counts {1: 3, 2: 2, 3: 7, 4: 3}); selected rank has 12 (DGP counts {1: 3, 2: 2, 3: 6, 4: 1}). No other requested primary-cause category occurs.
Fixed-rank time/unit failure rates are 0.273/0.409; selected-rank rates are 0.278/0.389.
All invalid N=50 split fits converged below the iteration cap, passed stationarity, and retained their supplied computational rank. Split objective stability is not assessed because the prescribed split stage uses one fit per half; blank stability fields are therefore not evidence of instability.
Per-block singular diagnostics show no A/B/H rank collapse. The historical fit records store only the maximum coefficient envelope across A, B, and H, not the block attaining it, so matrix-specific bound attribution cannot be reconstructed without replaying these N=50 fits; no such replay was authorized. This limitation is stated rather than guessed.

Local/plugin targets do not require split correction; broad targets do. Retention is:
- fixed_rank, broad_split_corrected: attempted 96, point retained 88, inference retained 24.
- fixed_rank, local_plugin: attempted 120, point retained 110, inference retained 110.
- selected_rank, broad_split_corrected: attempted 96, point retained 72, inference retained 24.
- selected_rank, local_plugin: attempted 120, point retained 90, inference retained 90.

## Full versus split conditioning

Across the failed bundles, the median full-panel envelope ratio is 0.881, versus 1.124 for invalid split fits. Median stationarity residuals are 2.36e-05 versus 4.62e-05; median sigma-r ratio is 1.000 versus 1.000; and median runtimes are 0.027s versus 0.040s. Split objective-stability gaps are unavailable by design, while full-panel multi-start gaps are saved in the conditioning table.
The failed bundles therefore pair well-behaved 50x50 full fits with smaller 50x25 or 25x50 fits whose coefficient envelopes cross B. Stationarity and singular-rank diagnostics remain satisfactory, making finite-small-panel envelope activity the descriptive mechanism. This is not evidence of asymptotic failure.

## DGP 4 rank-cap replay

- replication 3: cap 2 selected [1, 1, 1]; cap 3 selected [1, 1, 1].
- replication 4: cap 2 selected [1, 1, 2]; cap 3 selected [1, 1, 2].
- replication 5: cap 2 selected [1, 1, 2]; cap 3 selected [1, 1, 2].
Rank 3 is selected in 0/3 diagnostic replays.
The two former (1,1,2) cap hits remain (1,1,2), so cap=(2,2,2) does not appear to truncate these DGP-4 solutions. This does not automatically change or freeze the main cap.

## Fixed-rank size escalation

N=50 four-split success is 0.250; N=100 is 0.917. Time-split failure changes from 0.273 to 0.000; unit-split failure changes from 0.409 to 0.042.
At N=100, full-panel success is 12/12 with 0 bound hits. Point retention is 1.000; local/broad/total inference retention is 1.000/0.917/0.963. Gram and Riesz failure rates are 0.000/0.000. Total runtime is 80.956s (mean 6.746s per replication).
The sole N=100 split failure is the DGP-3 replication-7 unit half 2 coefficient-bound hit; it converged, passed stationarity, and preserved rank. No other failure mechanism occurs.
These 12 new evaluations are a numerical gate, not publication Monte Carlo evidence.

## Gate decision

The evidence supports proceeding to the medium diagnostic, but it was not launched.

## Reconciliation

Fit-level rows: 56; invalid rows: 27. Summary and retention denominators are generated directly from the consolidated run records.
