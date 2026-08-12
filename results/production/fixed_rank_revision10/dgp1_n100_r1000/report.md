# FP-DGP1-N100 supplied-rank production report

## Integrity and execution

The N=50 lightweight commit `5acf224947372b9afcf3a29105c75d9969623264` was pushed successfully before this run; `HEAD`, local `main`, and `origin/main` then agreed. The only pre-existing untracked content was the authorized N=50 raw evidence. The N=100 preflight matched the frozen manuscript, calibration, and source hashes; 182 tests passed; Ruff and `git diff --check` passed.

The balanced-sequence master seed remained `2026081103`. The 1,000 N=100 semantic task IDs and 1,000 derived DGP seed states were internally unique and had zero overlap with their N=50 counterparts.

This cell ran only DGP1, `N=T=100`, replications 1-1000, supplied rank `(1,1,1)`, with eight workers and one native thread per worker. The fixed supplied-rank path ran directly; no automatic rank selector, cap+1 pilot, IC, threshold, or path-derived rank was invoked. The engine fingerprint is `d5bf1eff50bc6db46bedee9ace9c7511fb7ff8365c9f189f80103c39b5c20eab` and config hash is `32bb33c3cfe3047d`.

Launch:

```text
.\.venv\Scripts\python.exe scripts\run_mc.py --config results\production\fixed_rank_revision10\dgp1_n100_r1000\config.toml --n-jobs 8
```

Resume:

```text
.\.venv\Scripts\python.exe scripts\run_mc.py --config results\production\fixed_rank_revision10\dgp1_n100_r1000\config.toml --n-jobs 8 --resume
```

There were no interruptions, resumes, corrupt bundles, or terminal task failures. The deterministic coefficient envelope was 6.509671, leaving 3.490329 inside `B=10` and exceeding the required `c_B=1` margin.

## Retention and inference

All 1,000 full fits were accepted and every target has 1,000 point-valid estimates. Local/fixed-time targets have 999 inference-valid records; broad corrected targets have 972. Across target records, frozen inference accounting produced 17,766 successes, 218 boundary/interiority suppressions, and 16 split-fit suppressions. Full fits had no rank loss. Twenty-six replications used constrained fallback/became boundary-active in at least one maintained fit.

Every one of 18,000 Riesz solves converged; residuals ranged from `5.35e-9` to `1.00e-8`. Tangent-Gram smallest eigenvalues ranged from 0.1024 to 0.3391 and condition numbers from 64.2 to 666.5. No Riesz, Gram, or variance failure occurred. All split coefficient fits converged and retained rank; two broad-target replications were suppressed because a split stationarity residual narrowly exceeded `1e-6`.

Selected target results follow; the full table is `summaries/by_target.csv`.

| Target | R inference | Bias | RMSE | Mean SE | SE/SD | Coverage (MC SE) |
|---|---:|---:|---:|---:|---:|---:|
| A entry | 999 | -0.0174 | 0.1134 | 0.0990 | 0.883 | 0.894 (0.0097) |
| A fixed-time mean | 999 | -0.0148 | 0.0867 | 0.0759 | 0.887 | 0.904 (0.0093) |
| A full mean | 972 | -0.0162 | 0.0187 | 0.0097 | 1.033 | 0.611 (0.0156) |
| B entry | 999 | 0.0057 | 0.2816 | 0.2231 | 0.792 | 0.886 (0.0101) |
| B fixed-time mean | 999 | -0.0072 | 0.1765 | 0.1491 | 0.844 | 0.898 (0.0096) |
| B full mean | 972 | 0.0107 | 0.0225 | 0.0217 | 1.100 | 0.944 (0.0073) |

The A and B fixed-time group contrasts remain weak-target stress diagnostics outside headline Assumption-10 coverage. No finite record was trimmed. Two inference-valid extremes were retained: B entry replication 388 (`error/SE=-8.26`) and B fixed-time group-contrast replication 639 (`error/SE=-55.17`).

## N=50 comparison

Relative to frozen N=50 evidence, N=100 full-fit/point retention rose from 88.9% to 100%; local inference retention rose from 86.8% to 99.9%; broad inference retention rose from 0% to 97.2%. Local bias and RMSE generally decreased, and local coverage improved but remains below 95%. Broad B-target coverage is near 95%, whereas broad A-level coverage remains poor despite SE/SD near one because its persistent bias is large relative to its shrinking sampling dispersion. This comparison is descriptive and caused no tuning.

## Performance and reproducibility

Engine wall time was 1,663.74 seconds and end-to-end command time about 1,802.7 seconds. Median/p95/max task times were 11.99/15.02/186.00 seconds; throughput was 0.601 tasks/second. Worker utilization was 98.42%, CPU use 32.91% of logical capacity, peak worker RSS 179,695,616 bytes, serialization 19.95 seconds, writes 6.97 seconds, and idle tail 33.10 seconds. Full-fit, split-fit, and aggregate inference totals were 859.78, 1,299.92, and 808.49 worker-seconds. Separate Riesz and variance clocks are unavailable in the frozen instrumentation, as documented in the field audit.

All planned table/figure scientific fields are present. No tuning, scientific-source modification, manuscript modification, calibration change, N=200/N=400 task, or other DGP run occurred.
