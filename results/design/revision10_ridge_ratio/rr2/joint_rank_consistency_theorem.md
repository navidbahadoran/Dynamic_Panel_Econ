# Manuscript-ready rank theorem

## Theorem (blockwise ridge-ratio rank consistency)

Suppose Assumptions `a:stab`, `a:exog`, `a:geometry`, `a:ned`, `a:moments`, `a:signal`,
`a:identification`, `a:gram`, and `a:growth` hold, the fixed true vector `r_0` belongs to the fixed reporting
set `R_max`, and a joint pilot satisfies

`max_M ||M_hat^pil-M_0||_op=O_p(sqrt(NT zeta_NT))`

conditionally on the common-shock sigma-field. Form the dimensionless spectra with the existing fitted-value
weights,

`lambda_hat_M,j=(w_M/w_A,1)^2 sigma_j(M_hat^pil)^2/(NT)`, `j>=1`,

set `lambda_hat_M,0=1` and `a_NT=1/log(NT)`, and define

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`,

`r_hat_M=min argmin_{0<=j<=bar_r_M} R_M,j`.

Then

`P_0(r_hat=r_0 | C) ->_p 1`.

Consequently, Theorems `thm:recovery`, `thm:target_expansion`, and `thm:twoway` remain valid when supplied
ranks are replaced by ridge-ratio selected ranks and all inferential objects are computed from the final
literal fixed-rank post-refit.

No Revision-9 IC penalty, threshold, candidate-coverage, candidate-refit, or IC numerical-error condition is
an assumption of this theorem.

## Proof

The existing weights divided by `w_A,1` are jointly bounded above and away from zero with probability tending
to one under the maintained moment, innovation, and prediction-identification conditions. Weyl's inequality,
the generic operator rate, and the strong signal condition therefore give order-one normalized positive
singular values and `O_p(zeta_NT)` values after truth.

The maintained growth condition gives `a_NT->0` and `zeta_NT/a_NT->0` on every admitted rectangular sequence.
For each positive-rank block, the ratio at truth tends to zero, pre-truth ratios (including the zero anchor)
are bounded below by a fixed positive constant, and post-truth ratios tend to one. For a zero-rank block, the
anchored ratio tends to zero and every positive-index ratio tends to one. The cap+1 spectrum supplies a
genuine stochastic post-truth value when truth is at the reporting cap. Thus each true block rank is the
unique asymptotic minimizer.

Let `E_M={r_hat_M!=r_M,0}`. The number `P+K+1` is fixed, so

`P_0(r_hat!=r_0|C)<=sum_M P_0(E_M|C)->_p0`.

This finite union proves the joint conclusion without growing-model-selection machinery.

## Attempted corollary (joint cap+1 implementation)

The desired corollary does not follow from Revision 9 as written. Its prediction-identification lower bound
only covers candidates within `R_max`. A cap+1 candidate can generate a difference of rank
`2bar_r_M+1`, outside the cap-truth reporting difference class. The basic inequality would give the desired
rate if curvature were extended to that enlarged class, but doing so is a substantive strengthening. The
generic theorem above remains valid; the requested implementation corollary remains unresolved.
