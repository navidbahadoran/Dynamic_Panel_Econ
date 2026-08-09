# Paper-to-code coefficient-estimator alignment

## Maintained mathematical problem

For supplied rank vector `r`, the paper estimator is

`argmin_{Theta in M(r), ||Theta||_max <= B} (1/(2NT)) sum_it [y_it-Y_it(Theta)]^2`.

The entrywise box is part of the finite-sample feasible set.  An unconstrained
optimum outside the box is therefore a trigger for constrained optimization,
not an estimator failure.  The current implementation does not solve that
boundary problem.

## Audit table

| Path | Code path | Current feasible set and optimizer | Current treatment of B | Classification | Required correction |
|---|---|---|---|---|---|
| Fixed-rank full panel | `estimation.fit_fixed_rank`; `rank_selection.fit_fixed_rank_multistart`; fixed branch in `monte_carlo.run_replication` | Exact supplied-rank factorization; unconstrained alternating row/time least squares | Computes `max_envelope_ratio` after fitting; full-panel result at or outside B becomes `coefficient_bound_hit` | `INTERIOR_EQUIVALENT_ONLY` | Keep a feasible, accepted unconstrained optimum as the fast path; invoke a literal entrywise-box-constrained fixed-rank fallback otherwise. |
| Selected-rank candidate post-refit | `rank_selection.fit_fixed_rank_multistart`; `rank_selection.fit_invalid_reasons`; `rank_selection.select_ranks` | Exact candidate rank; unconstrained multistart ALS | `coefficient_bound_active` makes a post-refit invalid and prevents it minimizing IC | `INTERIOR_EQUIVALENT_ONLY` | Use the same constrained fallback for each retained candidate; boundary activity must be diagnostic rather than automatic invalidation. |
| Rank-at-most-cap pilot | rank-adaptive routes in `rank_selection`; every joint route refit ultimately calls `fit_fixed_rank` | Rank-adaptive outer search, but each route is fitted by unconstrained ALS; starting matrices are rescaled into the box | Starts are made interior, but iterates may leave; outside fits are rejected as invalid | `INTERIOR_EQUIVALENT_ONLY` for an interior returned route; not the literal cap-bounded objective on boundary events | Each route/refit used by the cap-bounded pilot must solve its fixed-rank box-constrained subproblem, with boundary-feasible routes retained under the existing best-valid-route policy. |
| Time-half splits | `inference.prepare_split_fits.fit_one` -> `estimation.fit_fixed_rank` | Supplied full-panel rank on the time subset; one unconstrained ALS fit | Records envelope ratio; downstream split validity treats a hit as failure | `INTERIOR_EQUIVALENT_ONLY` | Apply the identical fast-path/fallback rule to both time-half coefficient fits. |
| Unit-half splits | `inference.prepare_split_fits.fit_one` -> `estimation.fit_fixed_rank` | Supplied full-panel rank on the unit subset; one unconstrained ALS fit | Records envelope ratio; downstream split validity treats a hit as failure | `INTERIOR_EQUIVALENT_ONLY` | Apply the identical fast-path/fallback rule to both unit-half coefficient fits. |

The nuclear path is screening only.  Its Dykstra proximal step does impose an
entrywise box together with nuclear shrinkage, but it is not the unpenalized
fixed-rank paper estimator and does not make any final/post-refit/split path
exact.

No audited final coefficient-estimation path is `EXACT_PAPER_ESTIMATOR` when a
box boundary can bind.  On a genuinely interior unconstrained optimum the same
solution is feasible and the implementations are numerically equivalent to the
constrained problem, which is the limited meaning of `INTERIOR_EQUIVALENT_ONLY`.

## Inference on a boundary-active estimate

No manuscript source is present in this Git working tree; the repository README
only states that the numerical routine solves the interior problem and treats
an envelope hit as failed interiority.  That is insufficient to establish an
explicit paper rule for inference after a valid finite-sample constrained
boundary solution.  This remains a paper/code clarification needed from the
author.  No alternative Riesz procedure or suppression rule is invented here.

## STOP-gate consequence

The requested constrained fallback is the necessary code correction for the
literal estimator.  It was not implemented in this audit because the task first
requires a finite deterministic DGP envelope, and expressly commands STOP if
calibration does not furnish one.  The `c_h` calculation fails that prerequisite.

