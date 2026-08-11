# Deterministic non-scientific solver test plan

## Test data policy

Use hand-constructed fixed matrices and deterministic algebraic regressors. Do not call DGP1-DGP4, use RR5 seeds, calibrations, or compare rank-recovery outcomes. Small sizes (`N,T` between 6 and 30) cover exact unit tests; representative `50x50`, `100x100`, and rectangular shapes are engineering integration tests.

## Mandatory correctness cases

| Case | Construction | Required assertions |
|---|---|---|
| Rank-one inside width four | Exact outer-product coefficient plus fixed regressors | finite solution; rank-at-most feasibility; objective no worse than start; frozen residual evaluated |
| Rank-two inside width four | Sum of two deterministic outer products | same, with both directions represented |
| Rank-zero block | One block identically zero and other blocks nonzero | zero product can be represented; dormant-direction setup preserves it; no rank threshold invoked |
| Badly scaled equivalent factors | `(U*1e8,V/1e8)` and reverse | products/objectives identical before and after allowed gauge; solver does not overflow |
| Nearly collinear columns | columns differing by deterministic `1e-10` perturbation | stable least-squares/QP result, finite telemetry, no normal-equation dependence |
| Lower-rank boundary | exact deficient product in full width | no forced full-rank assertion; zero padding preserves product; available dormant direction can be activated by a deterministic descent problem |
| Active box | optimum deliberately touches `+B` and `-B` | feasibility `<=1e-8`; frozen factor KKT `<=1e-4`; correct active signs |
| Interior optimum | strict envelope below `B-1e-8` | interior path retained; product projected gradient `<=1e-6` |
| Cross-block collinearity | deterministic regressors with near-dependent joint columns | stable minimum-norm solves and deterministic output |

## Transformation tests

1. For random-looking but fixed literal arrays, verify QR, orthogonal rotation, and balanced-SVD transformations preserve every product and objective.
2. Verify zero-padding and one-sided dormant directions add exactly zero product.
3. Apply transformations repeatedly and bound accumulated product/objective error.
4. Confirm internal diagonal equilibration is bijective: original and transformed QP solutions map to the same `x`, objective, constraints and KKT decision.
5. Demonstrate that non-orthogonal global scaling can change the finite factor-KKT norm; this guards against accidentally applying balanced scaling on the constrained path.
6. Verify identical coefficient products immediately before/after every accepted gauge normalization using an instrumented solver hook.

## Solver-path tests

- Assert loading rows, gauge step, and time rows execute in the frozen order.
- Compare the specialized QP against an independent high-accuracy reference on deterministic convex subproblems, including degenerate active sets.
- Verify warm and cold active-set starts reach the same objective/KKT solution.
- Require accepted coordinate updates and safeguarded full sweeps to be monotone under the existing allowance.
- Force an objective-increase rejection and verify that factors, reconstructed products and the stored objective all roll back to the same pre-sweep state.
- Verify a rejected accelerated step returns to the unaccelerated accepted iterate exactly.
- Evaluate stationarity with the existing production residual functions, not a test-only surrogate.
- Check the unchanged stop logic at values just below, equal to, and above `1e-6`, `1e-4`, and `1e-6` objective agreement.
- Verify all three maintained starts retain exactly their existing IDs and RNG-generated initial arrays.
- Verify exact-rank validation remains disabled only for the cap+1 pilot.
- Verify an unacceptable pilot returns `rank_selection_numerically_unresolved` and no spectra, ratios, ranks, final fit, or inference.

## Determinism and platform tests

- Repeat each case in one process and in a Windows-spawn worker with one native thread; compare serialized scientific fields.
- Repeat with task completion order deliberately reversed; canonical outputs must agree.
- Assert no nested process pool can be created from a worker.
- Capture factor norms, Gram/design condition proxies, QP iterations, active-set sizes and monotonicity only as diagnostics; none may alter branches except the specified deterministic solver safeguards.

## Acceptance gate for implementation

All tests above must pass, followed by the frozen resume-equivalence test and deterministic performance benchmark. Only then may a separately authorized tiny non-scientific integration run be considered. RR5 or any scientific Monte Carlo remains prohibited until a later explicit decision.
