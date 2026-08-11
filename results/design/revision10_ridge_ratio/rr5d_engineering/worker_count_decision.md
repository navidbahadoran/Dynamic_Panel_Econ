# Worker-count decision

## Recommendation

**Recommended default: `n_jobs=8`.** Keep `n_jobs=1` for debugging.

The benchmark used 24 identical deterministic non-scientific tasks, three repetitions per worker count, fixed 12-sweep kernels and one native thread. All 360 worker-count/task comparisons had identical scientific hashes. There were zero worker/process failures. Four short-budget pilot tasks were numerically nonconverged at every count; this invariant workload property is reported separately and was not used as a statistical criterion.

| Workers | Median wall (s) | Tasks/min | CPU % logical capacity | Worker utilization % | Peak worker RSS | Summed worker peaks | Idle tail (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 20.330 | 70.83 | 4.92 | 99.34 | 123.6 MB | 123.6 MB | 0.000 |
| 4 | 3.077 | 468.03 | 11.81 | 59.37 | 122.8 MB | 490.1 MB | 0.371 |
| 8 | 2.745 | 524.59 | 17.05 | 42.73 | 122.7 MB | 977.0 MB | 0.414 |
| 12 | 3.015 | 477.66 | 20.21 | 33.93 | 122.6 MB | 1.46 GB | 0.463 |
| 14 | 3.284 | 438.53 | 20.25 | 29.24 | 122.5 MB | 1.70 GB | 0.488 |

Eight workers had the best median wall time. Its 5% boundary was 2.882 seconds; no other count qualified. Counts 12 and 14 spent more on spawn/dispatch, memory and idle tail without improving throughput. The decision uses computation only—no DGP or rank-recovery outcome.
