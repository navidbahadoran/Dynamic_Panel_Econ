# FP-DGP1-N50 supplied-rank production report

## Scope and integrity

This is the single authorized DGP1, `N=T=50`, `R=1000` production cell. It used master seed `2026081103`, semantic replication indices 1 through 1000, supplied ranks `(1,1,1)`, eight workers, and one BLAS/OpenMP thread per worker. No DGP1 cell at a larger panel size and no DGP2-DGP4 cell was generated.

The pre-run gate passed at commit `2e66212f132ece3f96987f2f7c0398d81a9dc2f5`: the worktree was clean; local `main`, `HEAD`, and `origin/main` agreed; the manuscript and calibration hashes matched their frozen values; 182 tests passed; Ruff passed; and `git diff --check` passed. Exact searches of the worktree and all Git history found no earlier occurrence of master seed `2026081103`, so it was classified as fresh before any task seed was derived.

The launch command was:

```text
.\.venv\Scripts\python.exe scripts\run_mc.py --config results\production\fixed_rank_revision10\dgp1_n50_r1000\config.toml --n-jobs 8
```

The directly usable resume command is:

```text
.\.venv\Scripts\python.exe scripts\run_mc.py --config results\production\fixed_rank_revision10\dgp1_n50_r1000\config.toml --n-jobs 8 --resume
```

The engine fingerprint is `a570a03ad413263b80b526bcb3228801f23e544914d49260c7cdfd32a33280a9`; the resolved configuration hash is `e5824158072b2d7c`. There were no interruptions, resume events, corrupt bundles, or tasks skipped on resume.

## Estimator and DGP confirmation

The implementation's `fixed` execution path passed the true supplied rank `(1,1,1)` directly to the existing fixed-rank estimator. The scientific run is recorded semantically as `rank_mode=supplied_true_rank`. It did not call the Revision-10 ridge-ratio selector, cap+1 pilot, Revision-9 IC, threshold selection, candidate enumeration, or any other rank selector. Rank-sensitivity paths were disabled.

The run loaded the frozen DGP1 calibration file unchanged, used the diagonal DGP1 conditional-variance estimator, retained `B=10` and `c_B=1`, and used the frozen two-way correction for broad targets. The deterministic coefficient envelope was 6.149704, leaving 3.850296 to the box boundary and satisfying the required simulation margin without clipping or recalibration.

## Numerical and retention results

All 1,000 semantic tasks reached a terminal state. The DGP and target truths were generated in all 1,000 tasks. The full supplied-rank fit was accepted in 889 tasks and failed the frozen numerical gate in 111. Among the 889 accepted full fits, 21 were boundary-active for inference, no supplied-rank loss occurred, and the frozen point-valid rule retained 889 point estimates for every target.

For each of the ten local/fixed-time targets, inference was valid in 868 replications: 111 were excluded by full-fit failure and 21 by boundary/interiority. All Riesz solves on fitted tasks converged, no target was rejected for Riesz or empirical tangent-Gram conditioning, and all recorded variance estimates were finite.

For every one of the eight broad corrected targets, point estimates were retained in 889 replications but inference was retained in zero. Their 889 fitted replications were classified as 648 split-fit failures and 241 boundary/interiority failures under the frozen rules. The split coefficient computations themselves recorded only eight replications with a nonconverged split and no split rank loss; the much larger inference suppression is the maintained split numerical/interiority acceptance gate. Accordingly, coverage, mean SE, interval length, and SE/SD are undefined for all broad corrected targets in this `N=50` cell.

The complete target table is in `summaries/by_target.csv`. Selected local-target results are:

| Target | R point | R inference | Bias | RMSE | MC SD (inference set) | Mean SE | SE/SD | Coverage | MC SE coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A entry | 889 | 868 | -0.0407 | 0.1767 | 0.1721 | 0.1333 | 0.774 | 0.820 | 0.0130 |
| A fixed-time mean | 889 | 868 | -0.0345 | 0.1392 | 0.1342 | 0.0988 | 0.736 | 0.828 | 0.0128 |
| A group-1 fixed time | 889 | 868 | -0.0351 | 0.1396 | 0.1347 | 0.0996 | 0.739 | 0.831 | 0.0127 |
| A group-2 fixed time | 889 | 868 | -0.0338 | 0.1405 | 0.1356 | 0.0998 | 0.736 | 0.834 | 0.0126 |
| B entry | 889 | 868 | 0.0341 | 0.4502 | 0.4504 | 0.3034 | 0.674 | 0.809 | 0.0133 |
| B fixed-time mean | 889 | 868 | 0.0103 | 0.2655 | 0.2642 | 0.1879 | 0.711 | 0.840 | 0.0124 |
| B group-1 fixed time | 889 | 868 | 0.0091 | 0.2658 | 0.2643 | 0.1907 | 0.721 | 0.843 | 0.0123 |
| B group-2 fixed time | 889 | 868 | 0.0116 | 0.2730 | 0.2720 | 0.1911 | 0.703 | 0.843 | 0.0123 |

The A and B fixed-time group contrasts are reported only as weak-target stress diagnostics outside headline Assumption-10 coverage. Their coverage rates were 0.910 and 0.893, respectively, and they must not be treated as headline uniform-theorem evidence.

Rejection against zero (which is also the nonzero-truth power statistic under the replication-specific truths) is stored target by target. For the local/fixed-time levels it ranges from 0.521 for B entry to 0.932 for the A fixed-time means; for the weak contrasts it is 0.091 for A and 0.145 for B.

## Riesz, Gram, variance, and extreme records

Across 16,002 fitted target diagnostics, every Riesz solve converged. Riesz residuals ranged from `4.76e-9` to `1.00e-8`. The empirical tangent-Gram smallest eigenvalue ranged from 0.0231 to 0.3757, and its condition number ranged from 51.5 to 2588.7; no Gram failure was recorded. All 16,002 variance estimates were finite, so no variance failure was recorded.

Six inference-valid finite records satisfy `|error / SE| > 8`; none was trimmed or winsorized. They are B-entry replications 145, 638, and 901, and B fixed-time group-contrast replications 140, 419, and 786. Their values and standardized errors are retained in `summaries/extreme_records.csv` and in every reported statistic.

## Performance

The engine execution wall time was 1,399.86 seconds (23m20s), and end-to-end command time including final aggregation was approximately 1,540 seconds. Median, p95, and maximum task runtimes were 5.72, 35.02, and 188.66 seconds. Throughput was 0.714 tasks/second. Worker utilization was 99.19%; reported process CPU use was 30.84% of logical machine capacity. Peak worker RSS was 135,438,336 bytes. Serialization took 18.39 seconds, output writing 5.65 seconds, and the idle tail 13.66 seconds. The maximum pending queue was 15 tasks.

Full-fit diagnostic runtime summed to 1,891.84 worker-seconds and split-fit runtime to 4,884.45 worker-seconds. The current frozen instrumentation does not expose Riesz-solve and variance-construction time as separate clocks; total inference runtime was 309.35 worker-seconds. No estimator code was changed to add timing.

## Interpretation

This small-sample cell has three distinct conclusions. First, the supplied-rank full estimator is point-valid in 88.9% of attempted replications. Second, local/fixed-time inference is retained in 86.8%, but its conditional coverage is materially below 95% and mean SE understates the replication-specific error SD. Third, broad corrected inference is not operational at `N=T=50` under the frozen split/interiority gates: its inference retention is zero, so no conditional coverage claim is available.

These are evidence, not tuning criteria. No method, tolerance, calibration, seed, target, or failure rule was altered after results were observed. The prespecified DGP1 `N=100,200,400` sequence remains planned, but none of those cells was launched.
