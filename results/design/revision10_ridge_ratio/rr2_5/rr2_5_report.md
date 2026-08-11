# RR2.5 cap+1 pilot identification report

## Result

**RR2.5-PASS-MINIMAL.** The missing cap+1 pilot rate is available with one explicit, localized extension of
the existing prediction-identification domain. All accepted RR2 selector objects remain closed.

## Enlarged class

The pilot-only rank set is `R_max^+=product_M{0,...,bar_r_M+1}` and

`D_max^+={Theta-Theta_0: Theta has ranks in R_max^+, ||Theta||_max<=B}`.

Every difference block has rank at most

`bar_r_M+1+r_M,0<=2bar_r_M+1`.

The reported set remains `R_max`.

## Identification conclusion

A formally maintained global lower eigenvalue for the joint conditional regressor/prediction Gram would
derive cap+1 identification for all ranks. Revision 9 mentions such nonsingularity only as a sufficient
interpretation; its displayed assumption is restricted to `R_max`. Therefore no automatic derivation is
claimed.

The minimal repair is to require the same population prediction lower bound on `D_max^+` for the spectral
pilot. This is a mathematically stronger but localized domain extension of the existing identification
assumption. It is not a new stochastic or DGP-rate assumption.

## Concentration and pilot proof

On the enlarged class, nuclear duality introduces only the fixed factor
`sqrt(sum_M(bar_r_M+1+r_M,0))`. The covering entropy remains `C_+(N+T)log(Cb_NT)`, where `C_+` depends only on
fixed enlarged ranks. The same truncation, blocking, Bernstein, peeling, and Lipschitz arguments yield

`sup_(Delta in D_max^+\{0}) |sum uY(Delta)|/||Delta||
 =O_p(sqrt(NT zeta_NT))`

and, for a possibly smaller fixed `c_+>0`,

`sum Y(Delta)^2>=(c_+/2)||Delta||^2-CNT zeta_NT`.

For a feasible cap+1 pilot with normalized global objective gap `delta_NT=o_p(zeta_NT)`, loss expansion gives

`(1/2)Gamma_hat(Delta_hat,Delta_hat)
 <=S(Delta_hat)+NT delta_NT`.

Substituting the two enlarged-class bounds and solving the resulting quadratic yields

`||Delta_hat||_F=O_p(sqrt(NT zeta_NT))`,

hence `max_M||Delta_hat_M||_op` has the same order.

## Cost and unchanged objects

The manuscript cost is one pilot-only clause in `a:identification`, definitions of `R_max^+` and `D_max^+`,
and corresponding proof cross-references. Supplied-rank theorems retain the original domain. `B`, `c_B`,
coefficient support, `zeta_NT`, and rectangular growth are unchanged. No alternative extra-spectrum method is
needed.

## Source and scope

The analysis uses RR2 and the externally supplied manuscript
`E:\OneDrive\Desktop\ver7_revision9_Montecarlo_appendix_design.tex`, SHA-256
`2525d019ec5d7b28457585ff57c44670d7fafc6c199c60acc753a76e85070138`. The requested `(1)` filename was not
present. This phase contains theory documents only: no simulation, fit, implementation, tuning, scientific
source change, or manuscript edit.
