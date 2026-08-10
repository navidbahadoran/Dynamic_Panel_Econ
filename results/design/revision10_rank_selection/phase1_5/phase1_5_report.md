# Phase-1.5 theory report

## Conclusion

The Phase-1 ST2 gap is formalized but not closed. The correct classification is **ST3 — THEORY
COST TOO HIGH** for the proposed self-tuning route under the current manuscript theory.

## Gain representation

For block M, RSS improvement along a rank-one normal direction is exactly
`2<S_M,D>-Gamma_MM(D,D)`. After residualizing existing tangent coordinates, effective curvature
`mu_M` yields the local bound `Delta_loc<=||S_M^eff||_op^2/mu_M`. A joint post-refit gain is at
least this local gain, not automatically at most it; a reverse bound needs profiling and
localization. Endpoint RSS errors perturb a gain by at most twice the uniform endpoint error.

The block scores are `e_it y_i,t-ell`, `e_it x_it,k`, and `e_it` for A, B, and H. Conditional
heteroskedasticity and regressor-specific moments rule out one rate-sharp common scalar.

## Existing theory and envelope

Revision 9 establishes conditional uniform score order `sqrt(NT zeta_NT)`, empirical prediction
curvature with an additive `NT zeta_NT` term, and overfit RSS improvement `O_p(NT zeta_NT)` on a
fixed lattice. These are sufficient for a deterministic asymptotic penalty comparison but do not
give an observable probability-indexed operator-score quantile.

The best candidate uses a cap-residual block score, a spatially weighted estimate of the row
predictable covariance, a column variance proxy, empirical effective curvature, and a rectangular
matrix Freedman boundary with `x=log((N+T)L/alpha_NT)`. Exact recovery requires
`alpha_NT->0`; a fixed nominal level is insufficient. No numerical level or sequence is selected
here.

## Mixed states and order

In the original greedy state, an omitted L-block signal can enter the M-block score even when M
is correctly ranked. Current tangent residualization does not remove that omitted normal
direction. Greedy largest-only, add-all, and joint-maximum updates therefore lack coordinatewise
no-overfit control.

Blockwise nuisance-at-cap profiling resolves this definitionally: compare rank s and s+1 for M
while all other blocks may range up to their fixed caps. Decisions then become scalar by block and
order invariant. Dynamic profile completion avoids dependence on nuclear-screen omissions and
uses only a linear number of profile endpoints.

## Separation, zero ranks, and caps

Strong singular values and identification make every underfit profile class order NT away from
truth. Successive order-NT marginal gains still require profiled hereditary detectability. Once a
block contains truth, its gain is `O_p(NT zeta_NT)`. A true zero block is protected by the same
null comparison at s=0. Truth-at-cap is selected after all lower increments pass; truth outside
the cap remains outside the theorem and does not trigger automatic expansion.

## Assumption and inference consequences

Exactly two candidate new conditions remain: uniform observable block-envelope validity and
profiled marginal detectability. The first requires substantial new operator-valued spatial
variance/concentration theory and is not established. If both were proved, exact rank recovery
would transfer every supplied-rank recovery and inference theorem unchanged. Until then, the
conditional pseudocode in this package is research scaffolding only and must not be implemented.

No Monte Carlo, fit, source-code change, manuscript change, or Revision-9 outcome-based tuning was
performed.
