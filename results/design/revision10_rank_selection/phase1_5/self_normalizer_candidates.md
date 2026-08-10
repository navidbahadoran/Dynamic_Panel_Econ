# Observable self-normalizer candidates

Let `Theta_hat^cap` be the rank-at-most-cap profile fit and
`uhat_it^cap=y_it-Y_it(Theta_hat^cap)`. Using this common residual avoids inflating a scale merely
because the current tested state is underfit. Put `ghat_M,t=(z_it,M uhat_it^cap)_(i<=N)`.

## A. Residual MSE

- Sample: `sigmahat_cap^2=RSS(Theta_hat^cap)/(NT)`.
- Population: `(NT)^-1 sum_it E[u_it^2|C]`.
- Order: `O_p(1)` and plausibly consistent from existing recovery plus a scalar LLN.
- Dependence: valid for average variance under maintained mixing.
- Path validity: common cap residual is usable along the path if truth is inside the caps.
- Verdict: rejected as a complete normalizer. It ignores z, heteroskedasticity, spatial
  covariance, and the operator-norm maximum.

## B. Diagonal block quadratic variation

- Sample:
  `vhat_M^diag=max{max_i sum_t ghat_M,it^2, max_t sum_i ghat_M,it^2}`.
- Population: the same row/column sums of conditional marginal variances of `z_it,M u_it`.
- Order: typically `N max T` times a block variance, with logarithmic tail inflation.
- Consistency: needs uniform row/column LLNs and cap-residual replacement.
- Path validity: using cap residual avoids underfit contamination.
- Verdict: invalid under the maintained spatial dependence because covariance terms can govern
  the row-side matrix variance.

## C. Spatially weighted block variance

- Sample row proxy:
  `Vhat_M^row=sum_t [W_N elementwise (ghat_M,t ghat_M,t')]`;
  sample scalar proxy `vhat_M^row=||Vhat_M^row||_op`.
- Sample column proxy:
  `vhat_M^col=max_t sum_i ghat_M,it^2` (or a proved conservative predictable upper estimator).
- Population:
  `V_M^row=sum_t E[g_M,t g_M,t'|G_t]` and
  `v_M^col=max_t E[||g_M,t||_2^2|G_t]`, with
  `G_t=sigma(x_t,F_(t-1),C)`.
- Desired consistency: simultaneous conservative ratio control of the two predictable variance
  quantities over the fixed block/rank index set.
- Dependence: the spatial weights are compatible in spirit with the manuscript, but the existing
  HAC theorem treats fixed scalar Riesz directions, not this growing operator-norm covariance.
- Path validity: plausible with cap residual and truth inside caps; not proved.
- Verdict: best candidate, but unresolved.

## D. Matrix martingale self-normalization

Write `X_M,t=g_M,t e_t'`, so `S_M=sum_t X_M,t`. Sequential exogeneity makes the true-score
increments matrix martingale differences even with contemporaneous spatial dependence. A
rectangular matrix Freedman boundary would use

`v_M=max{||sum_t E[X_t X_t'|G_t]||_op,
          ||sum_t E[X_t' X_t|G_t]||_op}`.

The second term equals `max_t E||g_M,t||^2`. The realized row quadratic variation is
`sum_t g_M,t g_M,t'=S_M S_M'`; using its operator norm directly is tautological because it equals
`||S_M||_op^2`. A non-tautological empirical-Bernstein result therefore needs a predictable or
independently estimated spatial variance proxy, returning to candidate C.

## Candidate analytical boundary

If candidate C were proved conservative and score columns were truncated so
`||g_M,t||_2<=Rhat_M`, a rectangular Freedman form would be

`qhat_M(alpha)=[sqrt(2 vhat_M x)+(Rhat_M x)/3]^2/muhat_M`,

with `vhat_M=max(vhat_M^row,vhat_M^col)`, observable effective curvature `muhat_M`, and
`x=log((N+T)L/alpha)`, where `L=sum_M cap_M` is the finite comparison count.

This formula contains no fitted rank-recovery cutoff. It does contain a nominal familywise level
alpha. Rank consistency requires `alpha_NT -> 0`; a fixed alpha leaves nonvanishing false-addition
probability. A theorem could prescribe, for example, a polynomially vanishing alpha, but no
sequence is frozen in Phase 1.5 because the required conservative variance theorem is absent.
