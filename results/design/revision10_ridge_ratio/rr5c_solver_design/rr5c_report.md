# RR5c solver and simulation-engineering design report

## 1. Dominant numerical failure mechanism

The universal failure was the frozen multistart gate, driven primarily by per-start numerical invalidity: 177/180 pilots had fewer than two valid starts and 487/540 starts failed termination/stationarity. The width-four joint factorization combines gauge nonidentifiability, potentially ill-conditioned 12-column row/time systems and, on fallback, hundreds of thousands of generic SLSQP iterations over repeated convex QPs. Factor/Gram conditioning itself was not recorded, so specific conditioning channels are classified as plausible rather than asserted as observed fact.

## 2. Why objective-agreement relaxation is not the solution

Only three pilots reached two valid starts. Changing `1e-6` objective agreement would leave the other 177 unresolved and would violate the freeze. The baseline must first reliably produce individually valid solutions under `1e-6` interior projected gradient and `1e-4` constrained KKT.

## 3. Current solver conditioning diagnosis

Each RR5 block uses `N x 4` and `T x 4` factors; every row/time update jointly solves 12 variables. Factor rotations/scalings are nonidentified. Existing QR makes loadings orthonormal after a loading half-sweep but transfers triangular scale to time factors. Redundant/lower-rank directions and cross-block collinearity can make the alternating coordinates poorly conditioned. Product envelopes of 15.79 median, 90.87 p95 and 447.84 maximum before fallback are established; factor scaling alone cannot cause product blow-up, and the exact cause was not measured. Source inspection also found a latent interior rollback/bookkeeping defect: an increasing sweep can leave new factors paired with the prior objective. RR5 did not record whether that branch fired, so its role is unestablished.

## 4. Semantically exact options

Options assessed are product-preserving gauges, rank-revealing QR/SVD coordinate solves, dormant full-width direction handling, reversible equilibration, QP-specific active-set solves, warm active sets, cached invariant factorizations/layouts, monotone safeguards, deterministic line searches, reusable workspaces and richer condition telemetry. Each retains the coefficient products, objective and rank-at-most feasible set. Arbitrary non-orthogonal constrained-path gauge changes are excluded from the baseline because the finite factor-KKT norm is coordinate dependent.

## 5. Preferred architecture

Use deterministic gauge-aware safeguarded ALS plus quotient-aware Gauss-Newton refinement on the interior, and the same alternating product-box subproblems solved by a deterministic dense active-set QP method under reversible internal equilibration on fallback. Preserve the existing update order, maintained seeds/starts, widths, maximum-iteration and acceptance rules. Evaluate final residuals with the frozen production functions.

## 6. Semantic-equivalence argument

`(U,V)->(UG,V G^-T)` preserves `UV'`; balanced SVD and zero padding reconstruct the same product. The normalized loss and product box therefore remain identical. `x=Dz` is a bijection between original and preconditioned QPs. Interior tangent stationarity is product invariant. Exact constrained KKT zero is coordinate invariant, but its finite norm is not; mapping back to maintained coordinates before evaluation preserves the existing finite diagnostic. No penalty, cutoff, approximate ratio or feasible-direction deletion is introduced.

## 7. Changes ruled out without paper decision

Changing `B`, pilot/reporting caps, number or identity/distribution of starts, stationarity/KKT tolerances, objective agreement, one-start acceptance, singular-value rank thresholds, ridge, rank-zero anchor, exact ratio SVD, DGP/calibration, or the mathematical pilot estimator is prohibited.

## 8. Deterministic solver tests

The frozen suite covers rank-one and rank-two products inside width four, a zero block, badly scaled equivalent factors, nearly collinear columns, lower-rank boundaries, active boxes, interior optima, cross-block collinearity, transformation/product invariance, QP reference equality, monotonicity, frozen residual boundaries, maintained start/RNG identity, unresolved behavior, Windows spawn and order determinism. It uses no scientific DGP or RR5 outcomes.

## 9. Frozen resumability architecture

Every task has a canonical semantic ID and scientific fingerprint containing the code commit, config/calibration hashes, selector and master seed while retaining existing RNG semantics. Each task writes one validated, fsynced, atomically renamed bundle. A sole-writer atomic manifest tracks expected/running/completed/failed/unresolved. Resume verifies fingerprints and bundle hashes, quarantines corruption, skips all valid terminal states, requeues only missing work and aggregates in canonical semantic order.

## 10. Interrupted/resumed equivalence

At least 24 deterministic synthetic tasks are run uninterrupted and under deliberate Ctrl-C, worker death, parent death and reboot-equivalent restart. Exact discrete/serialized equality and existing accepted floating tolerances are required for IDs, seeds, statuses, products, spectra, ratios, ranks, diagnostics and summaries. Terminal failure preservation, stale-running recovery, corruption handling and no-recomputation are explicitly asserted.

## 11. Performance/parallel protocol

Every larger run is gated by measurements of wall/task/phase time, worker and system CPU, memory per worker/system peak, spawn, serialization, bundle writes, idle tail, queue depth and worker timeline. One native thread, one outer spawn pool, no nested pools, execution-order-independent seeds and deterministic aggregation are mandatory.

## 12. Safe performance improvements

Permitted after equivalence tests: worker-local immutable config/calibration caches, no repeated TOML parsing, pool reuse, dynamic semantic scheduling, smaller transport payloads, copy/workspace reduction, invariant half-sweep caches, exact factorization reuse, warm QP active sets, compact lossless bundles, direct atomic completion writes, verified read-only sharing, native-thread enforcement and observational telemetry.

## 13. Worker benchmark plan

Use at least 24 fixed synthetic tasks across representative square/rectangular dimensions and interior, constrained, lower-rank and post-refit kernels. Measure `1,4,8,12,14` workers with a warm-up and at least three rotated repetitions. Require scientific equivalence first. Select the smallest stable count within 5% of best median throughput, then lower memory and idle tail. No benchmark was run in RR5c.

## 14. Decision

**RR5c-PASS-SEMANTIC-ENGINEERING.**

## 15. Implication for ridge-ratio route

The ridge-ratio route remains viable but untested by RR5 because no pilot passed and no normalized spectrum or rank was produced. A later implementation phase must pass the frozen deterministic solver, resume-equivalence and performance gates before any separately authorized scientific rerun.

## 16–18. Scope confirmations

- No Monte Carlo, DGP generation or RR5 seed was run.
- No tuning occurred.
- Scientific source, optimizer, scripts, tests, configs, calibration and manuscript are unchanged.

## 19. Files created

Thirteen design/audit documents were created in `results/design/revision10_ridge_ratio/rr5c_solver_design/`:

1. `current_solver_mathematical_map.md`
2. `rr5_solver_failure_mechanism.md`
3. `semantically_exact_solver_options.md`
4. `preferred_cap_plus_one_solver.md`
5. `solver_semantic_equivalence.md`
6. `deterministic_solver_test_plan.md`
7. `resumability_freeze_spec.md`
8. `resume_equivalence_test.md`
9. `parallel_performance_freeze_spec.md`
10. `deterministic_performance_benchmark.md`
11. `safe_performance_improvements.md`
12. `rr5c_decision.md`
13. `rr5c_report.md`

## 20–22. Validation and Git

RR5c requires `git diff --check`; its result and the exact RR5c commit/local-remote equality are reported in the final handoff. Tests and Ruff are intentionally not rerun because no code changed. The frozen parent commit is `2536ce0a40510f9977697eb54083d5f7bd5cf5f0`; the intended commit message is `Design robust cap+1 solver and simulation engineering`.
