# Exact proposed Revision-10 ridge-ratio selector

## Source and scope

This definition uses the RR1 record and the externally supplied author manuscript
`E:\OneDrive\Desktop\ver7_revision9_Montecarlo_appendix_design.tex`, SHA-256
`2525d019ec5d7b28457585ff57c44670d7fafc6c199c60acc753a76e85070138`. The requested filename
with suffix `(1)` was not present. No repository manuscript is substituted.

Let the fixed block collection be

`M in {A^(1),...,A^(P),B^(1),...,B^(K),H}`,

and let each fixed reporting cap be `bar_r_M`, with `r_M,0 in {0,...,bar_r_M}`. The Monte Carlo
reporting caps already fixed by Revision 9 are `(3,3,3)` for its three-block design; in general the
paper retains its prespecified matrix-specific fixed caps. Reported ranks never exceed these caps.

## One joint spectral pilot (not yet theory-frozen)

Define `Theta_hat^pil` as an approximate global solution of

`inf L_NT(Theta)` subject to `rank(M)<=bar_r_M+1` for every block and
`||Theta||_max<=B`,

where `L_NT=(2NT)^(-1) sum_it {y_it-Y_it(Theta)}^2` is exactly the paper's joint loss. The extra
rank is spectral only. It is not an admissible reported rank.

## Dimensionless spectrum

Use the manuscript's existing weights

`w_A,l={ (NT)^(-1) sum_it y_i,t-l^2 }^(1/2)`,
`w_B,k={ (NT)^(-1) sum_it x_it,k^2 }^(1/2)`, and `w_H=1`.

Let `s_NT=w_A,1` and `q_M=w_M/s_NT`. For `j>=1`, define

`lambda_hat_M,j = q_M^2 sigma_j(M_hat^pil)^2/(NT)`.

This is the squared singular spectrum of the fitted-value-scaled coefficient block `q_M M_hat^pil`,
normalized by panel size. It is dimensionless. The reference uses an already-maintained weight and is
not outcome-dependent rule calibration.

## Ratios and rank

Set

`a_NT=1/log(NT)` and `lambda_hat_M,0=1`.

For `j=0,...,bar_r_M`, define

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`.

The block estimate and joint vector are

`r_hat_M=min argmin_{0<=j<=bar_r_M} R_M,j`,

`r_hat=(r_hat_A,1,...,r_hat_A,P,r_hat_B,1,...,r_hat_B,K,r_hat_H)`.

The `min` freezes exact numerical ties in favor of the smaller rank. The asymptotic minimum is unique,
so the convention has no first-order role.

## Final estimator and failure rule

After selecting `r_hat`, compute exactly one literal unpenalized joint constrained fixed-rank post-refit
`Theta_hat(r_hat)`. All reported coefficients, Riesz objects, targets, four split fits, split correction,
and variance objects are built from that post-refit. The spectral pilot is never reported as a coefficient
estimate.

If a future resolution validates the cap+1 pilot, then a pilot with no valid feasible numerical solution, a
failure of the frozen numerical tolerances, or unresolved frozen multistart diagnostics gives rank status
`pilot_numerically_unresolved`; no rank is
selected and selected-rank inference is not reported. Observable diagnostics do not certify the unknown
global objective gap.

The nuclear path is absent from the proposed statistical definition. It may be an optional deterministic warm start.
Grid changes cannot change the mathematical selector when they reach the same valid pilot solution.

The scale-normalized generic selector is exact and its consistency theorem is complete conditional on the
displayed pilot operator rate. The joint cap+1 construction has not been proved to supply that rate under the
unchanged Revision-9 identification domain; see `cap_plus_one_pilot_theorem.md`. Therefore this definition is
not authorized for implementation by RR2.
