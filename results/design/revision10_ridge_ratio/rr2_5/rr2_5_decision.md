# RR2.5 decision

## RR2.5-PASS-MINIMAL

The cap+1 pilot operator rate is proved after a localized extension of part (i) of prediction identification
from the reporting-cap difference class `D_max` to the fixed pilot-only class `D_max^+`.

The extension is mathematically stronger and is disclosed as such. It is a domain extension of the existing
identification inequality, not a new stability, exogeneity, dependence, moment, signal, coefficient-bound,
target, Gram, or panel-growth assumption. Uniform concentration extends with only fixed constants depending on
`bar_r_M+1`; `zeta_NT` and rectangular asymptotics are unchanged.

The complete approximate-optimization basic inequality gives

`||Theta_hat^pil-Theta_0||_F=O_p(sqrt(NT zeta_NT))`

and hence the required uniform block operator rate.

Decision: **RR2.5-PASS-MINIMAL**.

RR2 may now proceed from `RR2-PARTIAL` to `RR2-FREEZE`, conditional on adopting and clearly stating this
pilot-only identification-domain extension. No selector implementation or manuscript edit occurs in RR2.5.
