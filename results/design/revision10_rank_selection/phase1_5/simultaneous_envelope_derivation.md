# Simultaneous envelope derivation

There are

`L=sum_M cap_M`

profile increments. L is fixed because the number of matrices and all caps are fixed. No
dimension-growing model-count theory is needed.

For the true block score let `g_M,t=(z_it,M u_it)_i` and
`X_M,t=g_M,t e_t'`. Conditional sequential exogeneity yields
`E[X_M,t|G_t]=0`, where `G_t=sigma(x_t,F_(t-1),C)`. Contemporaneous spatial dependence is inside
the vector `g_M,t`, not excluded. Define

`v_M=max{||sum_t E[g_M,t g_M,t'|G_t]||_op,
          max_t E[||g_M,t||_2^2|G_t]}`.

For bounded increments with norm at most R, rectangular matrix Freedman suggests

`P(||S_M||_op > sqrt(2 v_M x)+(R x)/3 | C) <= (N+T) exp(-x)`.

Taking `x=log((N+T)L/alpha_NT)` and a union bound gives familywise probability at most
`alpha_NT`. Combining this with effective curvature gives the candidate gain envelope

`Ehat_M(s)=[sqrt(2 vhat_M x)+(Rhat_M x)/3]^2/muhat_M(s)
            +2 epsilon_RSS`.

This is an exact symbolic candidate, not a proved observable boundary. Three replacements need
proof under the maintained assumptions:

1. cap-pilot residual scores must replace unobserved innovation scores uniformly;
2. `vhat_M` must conservatively estimate the growing row-operator and column predictable
   variances under spatial mixing;
3. truncation/localization and observable `Rhat_M` must preserve the displayed tail probability.

The manuscript proves a blocked matrix `O_p` rate, not these replacements. Its scalar spatial-HAC
proof cannot be promoted to growing operator covariance by a finite union bound: the difficulty is
within each N-dimensional variance operator, not the fixed number L of comparisons.

A fixed nominal alpha is insufficient for exact rank consistency. One needs
`alpha_NT -> 0` and `log(1/alpha_NT)` small enough that `max Ehat_M=o_p(NT)`. A theorem-derived
deterministic function of N,T would be tuning-free under the prompt's definition. Bootstrap or
multiplier quantiles are not adopted: their validity for the conditional spatial/mixing matrix
maximum would require a separate high-dimensional bootstrap theory.

Thus finite multiplicity is fully resolved, but the observable single-comparison operator-tail
calibration is not.
