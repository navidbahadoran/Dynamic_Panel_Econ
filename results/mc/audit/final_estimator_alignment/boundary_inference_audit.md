# Boundary inference audit

## Scope and conclusion

This audit concerns the implemented computation when a valid fixed-rank estimate lies on the entrywise coefficient box. It does not extend the Revision-8 asymptotic theorem to boundary parameters.

The existing inference calculation can be evaluated at a valid boundary estimate. The implementation does **not** replace the weighted residual term by zero using an interior first-order normal equation. Boundary activity is now retained as an explicit diagnostic. The maintained theorem interpretation continues to require the paper's interiority condition.

## Code-path audit

The empirical tangent-space Riesz problem is constructed from the target direction, design arrays, and fitted tangent spaces. Its right-hand side is not formed from the regression residual. After the Riesz weights have been solved, the weighted-residual normal-equation quantity is computed and stored only as a diagnostic (`weighted_residual_identity`).

The point estimate remains the plug-in target `target_value(direction, fit.theta)`. The residual enters the influence and variance calculation through the Riesz weights; no line of the implementation algebraically substitutes a zero residual identity. The returned inference diagnostics therefore set:

```text
normal_equation_used_as_identity = false
```

## Boundary handling

A constrained solution is admissible for computation when the solver reports success, the box violation is within the numerical feasibility tolerance, and the constrained KKT residual is within its tolerance. Such a fit records `boundary_active = true` and `constrained_fallback_used = true`; boundary activity alone is not classified as a numerical failure.

The deterministic boundary fixture in `constrained_solver_tests.csv` confirms that inference is computed at a successful boundary solution while `normal_equation_used_as_identity` remains false.

## Interpretation

This is an implementation audit, not a new theorem. A reported result whose fitted coefficient is on the box boundary must be flagged for interpretation because Revision 8 obtains its stated asymptotic result under an interior coefficient-bound condition. The finite-sample software calculation remains defined, but theorem-based claims should continue to rely on the verified deterministic DGP envelope and its positive interior margin.
