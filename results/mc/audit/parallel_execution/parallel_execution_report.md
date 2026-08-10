# Parallel execution audit

The scientific section replays four locked Revision-9 realizations only. The scaling section is a deterministic computational workload using the same fixed-rank estimation kernel, including box-constrained estimation. Neither section is new Monte Carlo evidence.

## Architecture

Previously, pools were repeatedly created inside cell/chunk loops and only a small chunk could run concurrently. The revised executor uses one bounded Windows-spawn outer pool across DGP x N x T x semantic-replication tasks. Nuclear paths, cap pilots, candidate refits, local completion, and split fits remain sequential inside each worker; there is no nested process pool.

Configuration and frozen calibrations are initialized once per worker. Panel and design arrays are generated inside their owning worker instead of copied from a parent cache.

## A. Scientific equivalence

| Requested | Effective | Tasks | Wall seconds | Tasks/s |
|--:|--:|--:|--:|--:|
| 1 | 1 | 4 | 9763.337 | 0.000409696 |
| 4 | 4 | 4 | 3710.367 | 0.00107806 |

Scientific equivalence passed: `True`. The comparison covers semantic IDs, realization hashes, candidate sets, selected ranks, objectives, Q_hat, IC values, target estimates, statuses, and failure classifications.

The previously completed requested-8 replay also had effective worker count 4. Its checkpoint is preserved as `.superseded_scientific_njobs_8_effective4.json`, but its redundant timing is excluded from this redesigned benchmark and is not interpreted.

## B. Parallel scaling

| Requested | Effective | Tasks | Wall seconds | Tasks/s | Mean total CPU % | Peak RSS upper bound GiB |
|--:|--:|--:|--:|--:|--:|--:|
| 1 | 1 | 20 | 128.979 | 0.155064 | 4.92 | 0.118 |
| 4 | 4 | 20 | 34.469 | 0.580232 | 16.72 | 0.455 |
| 8 | 8 | 20 | 26.674 | 0.749783 | 22.85 | 0.903 |
| 12 | 12 | 20 | 25.646 | 0.779854 | 25.04 | 1.350 |
| 14 | 14 | 20 | 25.552 | 0.782732 | 25.84 | 1.571 |

Fastest scaling setting: requested/effective `14`/`14`, with `5.05x` wall-time speedup over serial.
All scaling outputs deterministic across worker settings: `True`.
One-native-thread policy passed for every detected BLAS/OpenMP library: `True`.
Per-process resident-memory measurement available: `True`.

Peak CPU was not sampled reliably; aggregate worker CPU time divided by wall time and 20 logical processors is reported as mean total CPU utilization. Peak RSS is a conservative sum of per-worker lifetime peaks, not a synchronized system-wide sample.

## Recommendation

For this 14-physical-core/20-logical-processor/approximately-40-GB machine, use `--n-jobs 12` for workloads with enough outer tasks. It is the smallest setting within 2% of the fastest wall time and therefore avoids extra worker memory for a negligible timing difference. Retain `--n-jobs 1` for deterministic debugging. Recheck Task Manager during the first separately authorized production launch because full selected-rank tasks are heavier than this scaling kernel.
