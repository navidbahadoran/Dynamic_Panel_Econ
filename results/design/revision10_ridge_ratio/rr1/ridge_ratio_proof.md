# Ridge-ratio consistency proof

For every block M, use the enlarged cap+1 pilot, set `a_NT=1/log(NT)`, define the anchor
`lambda_hat_M,0=1`, and compute

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT),
        0<=j<=bar_r_M`.

Select the smallest minimizer:

`r_hat_M=min argmin_(0<=j<=bar_r_M) R_M,j`.

## Positive truth

Let `0<r=r_M,0<=bar_r_M`.

- If `j<r`, both adjacent normalized singular values are bounded above and below by positive
  constants. Uniformly over the finite set, `R_M,j>=c_R>0` with probability tending to one.
- At `j=r`, the denominator is `Theta_p(1)` and the numerator is
  `O_p(zeta_NT)+a_NT=a_NT{1+o_p(1)}`. Hence `R_M,r=O_p(a_NT)->0`.
- If `j>r`, both spectral terms are `O_p(zeta_NT)=o_p(a_NT)`, so `R_M,j->p1`.
- At j=0, positive truth gives `R_M,0` bounded away from zero.

The result includes `r=bar_r_M` because the ratio at the reporting cap uses the pilot's
`bar_r_M+1` singular value, which is stochastic noise rather than mechanical zero.

## Zero truth

The anchor argument in `rank_zero_analysis.md` gives `R_M,0->p0` and every `R_M,j->p1` for j>=1.

## Joint rank vector

For each M, the true index is the unique asymptotic minimizer. The number of blocks and every cap
are fixed, so a finite union gives

`P_0(r_hat_M=r_M,0 for every M | C)->_p1`.

Equivalently, `P_0(r_hat=r_0|C)->_p1`. No growing-model-space, iid, Gaussian, homoskedastic,
independent-unit, independent-time, or relative N/T argument is used. All dependence enters only
through the already-proved pilot operator rate.
