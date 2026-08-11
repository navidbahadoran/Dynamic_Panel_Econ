# Parallel-efficiency audit

## RR5 execution

- Requested/effective jobs: 12/12; one native BLAS/OpenMP thread per worker; no nested pool.
- The driver creates one Windows-spawn `ProcessPoolExecutor`, installs config/calibration and thread limits once per worker, and keeps at most `2*n_jobs` futures pending. `FIRST_COMPLETED` provides dynamic scheduling.
- Panels are generated inside workers; large panel arrays are not pickled from the parent. Returned diagnostics are small records. Calibration/config payloads are pickled once per worker.
- Checkpoint parquet writes occur in the parent only when a three-task chunk completes, avoiding worker disk contention but delaying durability.
- Frozen calibration is parsed/calibrated once in the parent and the resulting 60-cell payload is installed once per worker; it is not reparsed per task. Each task does rebuild its panel, regressors, factor arrays and repeated row/time least-squares workspaces.

RR5 recorded 2.760 wall hours including interruption/resume. Summed task runtime was 106270.2 seconds, and summed pilot-start fit time was 104455.6 seconds (98.3% of task time). The task median was 485.3 seconds, p95 2070.2, and maximum 3471.8; max/median was 7.15. This heavy tail creates load-imbalance risk.

The accounting utilization proxy `sum(task runtime)/(12*wall)` is 89.1%. It is not CPU utilization and is depressed by killed in-flight work at the wrapper timeout. RR5 stored neither worker IDs/completion timestamps nor CPU/memory samples, so end-of-run idle time and peak memory cannot be measured exactly.

## Bottlenecks

The dominant cost is numerical: three width-(4,4,4) joint ALS/box fits per task, usually reaching 200 sweeps. Constrained fallback ran in 201 starts; its median row/time QP count was 514,983 and maximum 1,865,590. Before fallback, maximum coefficients reached 447.84, so the constrained path had materially difficult starting points even though it restored exact box feasibility. Repeated SVD summaries are small relative to those fits. DGP construction and serialization together account for roughly the remaining 1.7% of recorded task time.

Potential engineering overheads are Windows spawn startup, one-time config/calibration pickling to 12 workers, repeated construction of regressors and invariant index structures inside tasks, and five sequential checkpoint writes per chunk. None appears dominant relative to pilot fitting. Full RR5 memory per worker was not recorded.

## Worker-count conclusion

The prior deterministic kernel benchmark on this 14-physical-core/20-logical-processor/~40-GB machine found 12 workers at 25.646 seconds versus 25.552 for 14, with one-thread compliance. Its 12-worker mean CPU metric was 25.04% of 20 logical processors (about five logical cores on average), individual worker peak RSS was roughly 119-122 MB, and the conservative summed peak-RSS bound was 1.35 GiB. The benchmark was short and lighter than cap+1, so these are not RR5 resource measurements. They support 12 as a provisional engineering choice but do not robustly validate the heavier cap+1 workload because RR5 lacks CPU/memory/worker traces and suffered an interruption. Before a larger authorized phase, use a deterministic, non-DGP cap+1 engineering benchmark with resource telemetry; do not choose workers from rank outcomes.
