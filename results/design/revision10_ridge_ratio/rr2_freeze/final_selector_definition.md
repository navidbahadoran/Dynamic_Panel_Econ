# Final Revision-10 ridge-ratio selector

## Fixed rank sets

For the fixed coefficient-block collection

`J={A^(1),...,A^(P),B^(1),...,B^(K),H}`,

let `bar_r_M` be the prespecified reporting cap and `0<=r_M,0<=bar_r_M`. The admissible reported-rank set is

`R_max=product_(M in J){0,1,...,bar_r_M}`.

For the spectral pilot only, define

`R_max^+=product_(M in J){0,1,...,bar_r_M+1}`.

The number of blocks and all caps are fixed. There are no growing-rank asymptotics.

## One joint cap+1 spectral pilot

Compute one feasible approximate-global joint least-squares pilot `Theta_hat^pil` subject to
`rank(M)<=bar_r_M+1` for every block and `||Theta||_max<=B`, where

`L_NT(Theta)=(2NT)^(-1)sum_it{y_it-Y_it(Theta)}^2`.

Its normalized objective gap from the infimum over this feasible class must be `o_p(zeta_NT)`.

The extra rank exposes the `(bar_r_M+1)`st singular value. It is never reportable. Inputs are interpreted in
the stored canonical scientific units in which the unchanged coefficient box is defined.

## Frozen block weights and dimensionless spectrum

Use exactly

`w_A,l=[(NT)^(-1)sum_it y_i,t-l^2]^(1/2)`,

`w_B,k=[(NT)^(-1)sum_it x_it,k^2]^(1/2)`,

and `w_H=1`.

With common reference `s_NT=w_A,1`, define for `j>=1`

`lambda_hat_M,j=(w_M/w_A,1)^2 sigma_j(M_hat^pil)^2/(NT)`.

This is dimensionless. If `x_k^*=c_k x_k`, then `w_B,k^*=|c_k|w_B,k` and the equivalent coefficient is
`B_k^*=B_k/c_k`, leaving every normalized singular value unchanged. Under a coherent outcome-unit change
`y^*=d y`, the transformations of `A`, `B`, `H`, the weights, and `w_A,1` likewise cancel up to signs, which
do not affect singular values. Unit changes must be mapped coherently to the canonical box; holding the same
numerical box in different coefficient units would define a different constrained problem.

## Frozen ridge, anchor, ratios, and ranks

Set, without a multiplier,

`a_NT=1/log(NT)` and `lambda_hat_M,0=1`.

Define

`R_M,0=(lambda_hat_M,1+a_NT)/(1+a_NT)`

and, for `j=1,...,bar_r_M`,

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`.

Then

`r_hat_M=min argmin_(0<=j<=bar_r_M) R_M,j`,

where exact numerical ties go to the smaller rank, and

`r_hat=(r_hat_A,1,...,r_hat_A,P,r_hat_B,1,...,r_hat_B,K,r_hat_H)`.

Ranks are selected separately. No common latent rank, zero-rank test, cutoff, or sensitivity multiplier is
part of the definition.

## Final statistical estimator

After rank selection, compute exactly one literal constrained unpenalized joint fixed-rank post-refit
`Theta_hat(r_hat)`. Every reported coefficient, empirical Riesz object, target, four split fits, two-way
correction, variance, interval, and test is constructed from that post-refit. The spectral pilot is never a
reported coefficient estimator.

The nuclear path is not part of the statistical definition. It may only be optional deterministic
initialization for the cap+1 pilot or final post-refit. Its grid has no mathematical rank-selection role.
