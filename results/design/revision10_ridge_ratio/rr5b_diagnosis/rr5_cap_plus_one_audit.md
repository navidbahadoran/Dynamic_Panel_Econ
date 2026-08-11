# Cap+1 factorization audit

The cap+1 pilot calls the shared `fit_fixed_rank` kernel three times at widths `(4,4,4)`. A width-four factorization parameterizes matrices of rank at most four; exact numerical rank is deliberately not required by the acceptance gate. Every stored RR5 start nevertheless had numerical-rank vector `[4,4,4]`, so there is no recorded rank collapse.

Static implementation findings:

- Each block has an `N x 4` loading and `T x 4` factor (`estimation.py`, `_initial_blocks`, lines 100-130).
- Random starts use independent normal entries scaled by 0.1. There are no intentionally zero factor columns.
- Every loading update solves a joint 12-column least-squares problem, followed by QR renormalization of each block; every time update solves the corresponding joint 12-column problem.
- QR moves scale from loadings into factors without changing coefficient matrices. No recorded evidence indicates an invariance/normalization mismatch.
- The literal box fallback solves alternating linear box-QP subproblems. It was invoked for 201/540 starts; 162 of those failed the frozen `1e-4` KKT requirement.
- Before fallback, those 201 unconstrained solutions exceeded `B=10`: their median maximum coefficient was 15.79, p95 was 90.87, and the maximum was 447.84. The constrained solver restored literal feasibility in every case, but used a median 514,983 row/time subproblem iterations (maximum 1,865,590).
- 339 starts remained on the interior ALS path; 325 failed its `1e-6` projected-gradient criterion.
- All 540 objectives, final coefficient envelopes, numerical ranks, and collapsed singular summaries were finite. The maximum final absolute coefficient was `10.000000000000004` and maximum recorded box violation was `3.55271e-15`.
- Collapsed `sigma_r/sigma_1` across starts ranged from `0.192596` to `0.839314` (median `0.529876`), inconsistent with wholesale numerical rank collapse.

Not recorded: factor-column norms, ALS design condition numbers, per-block full singular spectra, intermediate factor ranks, or detailed constrained-subproblem messages. Those pathologies therefore cannot be affirmatively excluded, but there is no stored symptom of a cap/dimension mismatch. The evidence supports the narrower conclusion that the three full-width joint fits are substantially harder to bring to the maintained stationarity and common-objective gate than historical rank-one fixed fits.
