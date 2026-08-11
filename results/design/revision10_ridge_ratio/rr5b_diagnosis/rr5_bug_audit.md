# RR5 implementation-bug audit

No implementation bug is established from code and committed RR5 evidence.

| Question | Expected logic | Actual logic / evidence | Finding |
|---|---|---|---|
| Gate applied to wrong object? | Validate each cap+1 start, then compare valid starts | `fit_revision10_spectral_pilot` validates each returned `FitResult`; all diagnostics map one-to-one to 540 observer rows | No mismatch |
| Wrong caps? | Reporting `(3,3,3)`, pilot `(4,4,4)` | All 540 `requested_rank` and `numerical_rank` values are `[4,4,4]` | Correct |
| Exact-rank condition leaked in? | Rank-at-most pilot must allow collapse | Pilot calls `fit_invalid_reasons(... require_exact_numerical_rank=False)` | Correct |
| Residual on wrong scale? | Use path-appropriate frozen diagnostic | Interior fits use normalized tangent projected gradient with `1e-6`; constrained fits use factor-space box KKT with `1e-4`; `fit_invalid_reasons` consumes the precomputed constrained pass flag | No rejection-scale bug found |
| Objective normalization mismatch? | Same normalized LS objective across starts | Every path uses SSE divided by `2NT`; stability denominator is `max(1, abs(best))` | Correct |
| Final-fit rule applied to pilot? | Final multistart only after selection | Pilot raises at line 1494; final post-refit is below spectra/ratio code | No leakage |
| Revision-9 IC/status leakage? | None | Revision-10 dispatch calls `select_ranks`; no IC/candidate code is reached | None |
| NaN/equality handling? | Fewer than two valid starts must fail | Second objective/gap are NaN by construction, while `stable = len(valid)>=2 and ...` safely evaluates false | Intended |
| Worker serialization altered booleans? | Preserve diagnostics | Raw records, consolidated records and 540 fit rows agree on 180 failures and start diagnostics | No evidence |

One reporting limitation matters for interpretation but did not change decisions: `stationarity_residual` serializes two different metrics under one column and omits the explicit `stationarity_type`. The fallback flag and solver status make the correct thresholds reconstructible. Treating every residual against `1e-6` after the fact would be wrong, but the live gate did not do that.

The pilot's factor width is fixed at four rather than rank-adaptively reduced. This still parameterizes the mathematical rank-at-most-four set and exact rank is not required. It is a numerical architecture choice, not evidence that the frozen estimator was coded as a different statistical selector.
