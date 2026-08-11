# RR5 cap+1 solver failure mechanism

## Dominant established mechanism

RR5 failed at the numerical credibility gate, before ridge-ratio selection. All 180 pilots failed multistart objective stability. Of these, 177 lacked two individually valid starts; only three had two valid starts and then missed the `1e-6` agreement requirement. Across the 540 starts, 487 failed termination/stationarity. Therefore loosening objective agreement would not address the dominant failure.

The width-four joint pilot is also the dominant computational cost: 98.3% of summed task runtime. The constrained fallback was invoked 201 times; although it restored product feasibility in every case, 162 ended in constrained optimality failure. Its median and maximum accumulated row/time subproblem iteration counts were 514,983 and 1,865,590. On the interior path, 325/339 starts failed the projected-gradient rule.

## Mechanism classification

| Suspected mechanism | Classification | Evidence and limit |
|---|---|---|
| Three full-width joint factorizations are much harder than historical rank-one fits | **ESTABLISHED** | 487/540 cap+1 starts failed termination/stationarity; historical R9 rank-one full fits had 58/75 residual passes. Designs differ, so the comparison is descriptive. |
| Factor rotation/scaling nonidentifiability | **ESTABLISHED** mathematically | The invariant transformation `(U,V)->(UG,V G^-T)` is present for every block. Existing QR fixes only one gauge after loading updates. |
| Redundant directions and rank-deficient factor Grams can make factor coordinates ill-conditioned | **PLAUSIBLE** | The pilot width exceeds the data-generating ranks, and the parameterization admits lower-rank strata. RR5 did not store factor Gram spectra. |
| RR5 solutions numerically collapsed to lower rank | **NOT SUPPORTED** | Every stored numerical-rank vector was `[4,4,4]`; collapsed `sigma_r/sigma_1` was at least 0.1926. |
| Near-collinearity in the 12-column joint row/time systems | **PLAUSIBLE** | Joint designs combine all blocks and can be ill-conditioned; no condition numbers were recorded. |
| Existing QR can transfer extreme scale into `V` and worsen the next loading system | **PLAUSIBLE** | Algebra and update order support the channel, but factor norms were not recorded. |
| Poor scaling propagates across `A`, `B`, and `H` through joint updates | **PLAUSIBLE** | All blocks enter each joint solve. Cross-block design condition numbers were not stored. |
| Pure factor scaling by itself caused coefficient products of 15, 90, or 400+ | **NOT SUPPORTED** | Gauge scaling leaves the product unchanged. Product blow-up must arise from the fitted direction itself; ill-conditioned/underdetermined joint solves are a plausible contributor. |
| Unconstrained product excursions make fallback initialization and active constraints difficult | **ESTABLISHED** for the excursions and workload; **PLAUSIBLE** as causation | Fallback starts had median/p95/max unconstrained envelopes 15.79/90.87/447.84 and extremely high subproblem iteration counts. No controlled solver comparison isolates causality. |
| Generic SLSQP is inefficient for the repeated convex quadratic/linear-box subproblems | **ESTABLISHED** as architectural mismatch; performance gain remains **PLAUSIBLE** | The objective is quadratic and constraints linear, yet every row/time call uses a general nonlinear optimizer. High iteration totals are recorded. |
| Convergence slows specifically near lower-rank boundaries | **PLAUSIBLE**, not established | Bilinear geometry is singular at rank-deficient points, but RR5 final products were numerically rank four and no trajectory rank data were stored. |
| Box feasibility caused pilot rejection | **NOT SUPPORTED** | All 540 final products were feasible; maximum violation was `3.55e-15`. |
| Nonfinite arithmetic or objective failure | **NOT SUPPORTED** | All objectives and final coefficient products were finite. |
| Acceptance thresholds were applied on inconsistent scales | **NOT SUPPORTED** | The live gate correctly used `1e-6` for interior projected gradients and the stored `1e-4` constrained KKT pass. |
| Wrong caps, exact-rank leakage, IC leakage, or a selector bug | **NOT SUPPORTED** | Caps were `(4,4,4)`, exact-rank validation was disabled, and the exception preceded spectra, ratios, ranks, and final fit. |
| Interior objective-increase branch returns post-sweep factors with the prior stored objective | **ESTABLISHED** as a latent code-path defect; RR5 role **NOT SUPPORTED** | Source exits without rollback or appending the increased objective. RR5 did not store the needed history/branch flag, so it cannot be attributed to the 180 failures. |

## Conditioning diagnosis

The numerical representation has more degrees of freedom than the coefficient products because every block has an invertible gauge. Near lower-rank strata, those gauges become singular and ordinary alternating least squares can alternate through flat or badly scaled coordinates. The current QR orthogonalizes only `U`; the singular scale is accumulated in `V`, so the next row-design can remain poorly scaled. The three blocks are then concatenated into one solve, allowing within- and cross-block near-dependence.

For fallback, every one of (N+T) row/time updates per sweep calls SLSQP on a small convex QP. Constraint matrices change with the opposite factors and may inherit their scaling. The code warm-starts from the previous coordinate value, but it does not use a QP-specific active-set method, explicit variable equilibration, or cached linear algebra. This explains the *possibility* of slow progress and large iteration counts without claiming an unobserved condition number.

## Conclusion

RR5b supports a numerical-solver diagnosis, not a statistical-selector diagnosis. The newly documented rollback/bookkeeping hazard must be repaired as ordinary semantically exact solver engineering, but its RR5 causal role is unknown. The evidence is sufficient to motivate a redesign, but not to claim any one unrecorded conditioning pathology as proven. Deterministic, non-scientific tests must isolate those mechanisms before implementation is accepted.
