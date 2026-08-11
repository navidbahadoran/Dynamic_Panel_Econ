# RR1 ridge-ratio feasibility report

## Result

**RR1-PASS.** A blockwise ridge singular-value-ratio selector is provably feasible under the existing
Revision-9 rectangular-panel assumptions, with no simulation-selected constant and no new substantive DGP
assumption.

## Construction

Compute one joint spectral pilot under fixed caps `bar_r_M+1`. Keep reported ranks in
`{0,...,bar_r_M}`. Define

`lambda_hat_M,j=sigma_j(M_hat^pil)^2/(NT)`, `lambda_hat_M,0=1`, and
`a_NT=1/log(NT)`.

For `j=0,...,bar_r_M`, set

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`

and select the smallest minimizer block by block.

## Proof summary

The enlarged-cap pilot inherits the current `O_p(sqrt(NT zeta_NT))` operator rate. Weyl implies
`lambda_hat_M,j=Theta_p(1)` at and below each positive truth rank and `O_p(zeta_NT)` above it. The exact
rectangular growth algebra gives `zeta_NT/a_NT->0` without balancing N and T.

For positive rank r, ratios below r stay bounded away from zero, the ratio at r tends to zero, and ratios
above r tend to one. If r equals the reporting cap, the cap+1 pilot supplies the required noise singular
value. For true rank zero, the anchored ratio at zero tends to zero and all positive-index ratios tend to one.
The true rank is therefore the unique asymptotic blockwise minimizer. Fixed block count and caps give exact
joint rank-vector recovery by a finite union.

## Assumptions and inference

Sequential exogeneity, predetermined regressors, conditional spatial dependence, mixing/NED,
heteroskedasticity, common shocks, and unequal rectangular growth all remain allowed because the ratio proof
uses only the established pilot rate and signal condition. The only localized changes are cap+1 pilot
feasibility and the same `o_p(zeta_NT)` normalized objective-gap requirement.

On `r_hat=r_0`, the final unpenalized selected-rank post-refit is identical to the supplied-true-rank object.
The current recovery, Riesz, target expansion, two-way correction, spatial variance, and feasible normality
theorems therefore transfer exactly as in the existing selected-rank theorem.

## Scope

The design substantially reduces the conceptual selection workload to one enlarged pilot, spectra, ratios,
and one final selected-rank post-refit. Nuclear screening may remain for initialization but no longer decides
ranks. This is theory feasibility only: no selector, fit, Monte Carlo, scientific code, or manuscript was
changed.
