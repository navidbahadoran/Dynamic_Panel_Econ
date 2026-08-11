# Performance instrumentation

Future long-run task bundles record worker wall/process CPU time, pilot time, final-fit time, worker PID, current/peak RSS, native threadpool state, coordinator completion time, pending count, serialization time and atomic-write time.

The final run manifest records total execution wall time, worker count, one-thread setting, logical-capacity CPU utilization, worker-time utilization, peak worker RSS, spawn/setup proxy, summed pilot/final-fit time, serialization/write time, idle-tail duration and maximum pending count. Full per-task metrics are written separately from scientific aggregates.

The outer architecture remains one Windows-spawn process pool, no nested pool, bounded dynamic scheduling, execution-order-independent semantic seeds and one native BLAS/OpenMP thread per worker. Immutable resolved config and frozen calibration payloads are installed once per worker. Task-owned panels are generated in workers, avoiding large parent-to-worker panel transfers.

Instrumentation is observational. No timestamp, PID, CPU, memory or queue measurement controls scientific computation or acceptance.
