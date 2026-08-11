# RR5 evidence map

This map covers only evidence committed in RR5. No panel or estimator was reconstructed.

| Diagnostic | Committed source | Granularity | Exact content / limitation |
|---|---|---|---|
| Pilot objectives by maintained start | `rr5/beb5729f4f89ab13/fit_diagnostics.parquet` and `.csv`; duplicated as JSON arrays in `raw/*.parquet` and `rr5_replication_records.*` | 540 starts; 180 replications | Initial/final objective, iterations, and start number. |
| Stationarity/KKT residual and pass flag | `fit_diagnostics.parquet` | Start | `stationarity_residual`, `stationarity_pass`, fallback flag and solver status. The serialized field combines the interior projected-gradient metric with constrained factor-space KKT; `stationarity_type` itself was not serialized. |
| Feasibility and literal box | `fit_diagnostics.parquet`; `raw/*.parquet`; `rr5_feasibility_audit.csv` | Start and replication | Envelope ratio, maximum coefficient, violation, fallback and boundary flags. Constraint-message text was not serialized. |
| Termination | `fit_diagnostics.parquet` | Start | Convergence flag, iteration-cap flag, constrained solver status, iterations. Detailed constrained subproblem messages are missing. |
| Numerical ranks | `fit_diagnostics.parquet` | Start | One three-block vector per start; all are `[4,4,4]`. |
| Singular values | `fit_diagnostics.parquet` | Start | Only collapsed `sigma_1`, `sigma_r`, and their ratio (minimum summaries across positive-rank blocks). Full blockwise singular spectra and factor matrices were not stored. |
| Normalized spectra / ridge ratios | `rr5_block_diagnostics.csv`, block audits | Replication/block | Explicitly marked not reached. No values were recorded because the pilot gate raised before normalization. |
| Task elapsed runtime | `rr5_replication_records.*`, `attempted_replications.*`, `raw/*.parquet` | Task | End-to-end replication runtime. No start/end timestamps. |
| Pilot start runtime | `fit_diagnostics.parquet` | Start | Per-fit runtime; constrained runtime and subproblem iteration count also present. |
| Wall runtime | `rr5_runtime_summary.csv`, `rr5_run_manifest.json` | Run | Includes the wrapper interruption and resume. |
| Interruption / restart | `rr5_interruption_log.json` | Run | 117 preserved attempts, 63 completed after resume, zero checkpointed attempts rerun. This log was created for RR5, not automatically by the driver. |
| Task/checkpoint files | `raw/`, `rank/`, `fit/`, `inference/`, `replications/` | Three-replication chunk | Five parquet files per chunk. Empty rank/inference files are intentional. |
| Worker/process information | `run_manifest.json`, `resolved_config.*`, `rr5_environment.txt` | Run | Requested/effective workers and native-thread policy. Worker IDs, task-to-worker map, CPU samples and RR5 memory metrics were not recorded. |
| Configuration and command | `rr5_config.toml`, `resolved_config.*`, `command.txt`, manifests | Run | Frozen caps/tolerances, seed, jobs and command. |
| Output hashes | `rr5_output_hashes.csv` | File | SHA-256 inventory of RR5 evidence. |

Missing diagnostics are not inferred: per-block rejected-pilot spectra, full coefficient/factor arrays, ALS design condition numbers, zero-column norms, detailed QP messages, worker PIDs, CPU/memory traces, task timestamps, and automatic progress-state manifests do not exist in RR5.
