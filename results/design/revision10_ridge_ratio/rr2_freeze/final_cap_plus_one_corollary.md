# Final cap+1 spectral-pilot corollary

## Pilot-only identification domain

Define

`D_max^+={Theta-Theta_0: Theta has a rank vector in R_max^+, ||Theta||_max<=B}`.

For the spectral-pilot result only, extend part (i) of conditional prediction identification to require a
fixed `c_+>0` such that

`E_0[sum_it Y_it(Delta)^2|C]>=c_+||Delta||_F^2`

for every nonzero `Delta in D_max^+`. The constant may be smaller than the original reporting-class constant.
This is a localized strengthening of the domain of an existing identification assumption. It is not claimed
to follow logically from Revision 9's restricted statement, and it is not a new stochastic, dependence,
moment, signal, or growth assumption.

For every enlarged difference,

`rank(Delta_M)<=bar_r_M+1+r_M,0<=2bar_r_M+1`.

All such ranks are fixed. Consequently, the uniform score and empirical prediction bounds extend with only
fixed nuclear-norm and covering constants:

`sup_(Delta in D_max^+\{0}) |sum_it u_it Y_it(Delta)|/||Delta||_F
 =O_p(sqrt(NT zeta_NT))`,

`sum_it Y_it(Delta)^2
 >=(c_+/2)||Delta||_F^2-CNT zeta_NT`.

## Corollary

Let the feasible joint cap+1 pilot satisfy

`delta_NT=L_NT(Theta_hat^pil)-inf_(Theta in M_cap+1)L_NT(Theta)=o_p(zeta_NT)`.

Truth is feasible. With `Delta_hat=Theta_hat^pil-Theta_0`, `n=NT`,
`Gamma_hat(Delta,Delta)=sum_it Y_it(Delta)^2`, and
`S(Delta)=sum_it u_it Y_it(Delta)`, approximate optimality and the exact squared-loss expansion give

`(1/2)Gamma_hat(Delta_hat,Delta_hat)<=S(Delta_hat)+n delta_NT`.

Substitution of the enlarged-class concentration bounds yields

`(c_+/4)||Delta_hat||_F^2
 <=O_p(sqrt(n zeta_NT))||Delta_hat||_F
   +O_p(n zeta_NT)+n delta_NT`.

Because `delta_NT=o_p(zeta_NT)`, solving this quadratic gives

`||Theta_hat^pil-Theta_0||_F=O_p(sqrt(NT zeta_NT))`.

Therefore

`max_M||M_hat^pil-M_0||_op=O_p(sqrt(NT zeta_NT))`.

This supplies the generic pilot rate. Feasibility, multistart agreement, objective stability, and KKT
residuals are observable numerical diagnostics; they do not certify the unknown global infimum appearing in
`delta_NT`.
