# Target-regularity audit

## Analytic support calculation

With `rho_g=0.5` and innovations in `[-sqrt(3),sqrt(3)]`, the infinite-horizon deterministic AR support is

\[
 |g_{b,t}|\leq \frac{\sqrt{1-0.5^2}\sqrt 3}{1-0.5}=3.
\]

Starting at zero, after `k` transitions the sharper support is `[-3(1-0.5^k), 3(1-0.5^k)]`. Consequently

\[
f_{b,t}=0.6+0.2g_{b,t}\in[0,1.2]
\]

in the uniform-in-horizon closure. At a finite `k` its lower bound is `0.6*0.5^k`, which converges to zero and therefore is not a deterministic positive floor uniform in `T`.

This calculation, rather than realized projection ratios, determines the B-entry result. The entry target needs both unit and time leverage. The bounded positive `lambda_b` supplies unit leverage, but `f_b,t` does not supply a uniform positive time-leverage floor at the fixed paper date. Thus B-entry is not theorem-covered at `kappa_f_b=0.20`.

For time-averaged targets, the relevant time direction is an average rather than one fixed coordinate. The stationary mean of `f_b` is 0.6, so its sample average is separated from zero with probability approaching one under the bounded geometrically mixing AR process. This supports the time-averaged classifications below; it does not rescue the fixed B-entry.

## Classification

| Paper target | Classification | Analytic reason |
|---|---|---|
| A entry | THEOREM_COVERED | `f_a=0.5+0.1g_a` has support `[0.2,0.8]`; the loading also has a positive bounded floor. |
| B entry | NOT_UNIFORMLY_COVERED | `f_b` has support touching zero, so fixed-date time leverage lacks a uniform positive floor. |
| Overall fixed-time A mean | THEOREM_COVERED | The all-unit direction has nonzero loading alignment; A also has a positive time-factor floor. |
| Overall fixed-time B mean | THEOREM_COVERED | The all-unit direction has loading alignment bounded away from zero even when the fixed-date factor is small. |
| DGP1-3 fixed-time group means | THEOREM_COVERED | Each positive group indicator has nonzero alignment with the bounded positive loading vector. |
| DGP1-3 fixed-time group contrasts | NOT_UNIFORMLY_COVERED | Artificial halves have the same loading law; their loading-mean contrast is only stochastic order `N^-1/2`, and fixed-date time alignment cannot uniformly replace it. |
| DGP4 fixed-time group means | THEOREM_COVERED | Each group has positive loading alignment. |
| DGP4 fixed-time group contrasts | THEOREM_COVERED | The loading-mean gaps converge to the nonzero prespecified gaps (0.2 for A and 0.4 for B). |
| Full-panel A mean | THEOREM_COVERED | Both the loading mean and time-factor mean are nonzero. |
| Full-panel B mean | THEOREM_COVERED | Both the loading mean and the ergodic time-factor mean (0.6) are nonzero. |
| Group-specific time averages | THEOREM_COVERED | Positive group loading alignment and nonzero factor time means provide regularity. |
| Time-averaged group contrasts | THEOREM_COVERED | In DGP1-3 the time-average direction aligns with the nonzero factor mean; in DGP4 the nonzero loading gap additionally supplies unit alignment. |

## Code-label comparison

The code already excludes `B_entry` and DGP1-3 fixed-time group contrasts from headline theorem validation. It labels all other listed targets theorem-covered, consistent with this audit. The paper should explicitly state that the DGP1-3 fixed-time contrasts are stress targets outside the target-regularity assumption and that B-entry is outside the uniform theorem at the stated `kappa_f_b`.

Changing `kappa_f_b` alone is not the cleanest repair because the stationary support of `mu_f_b+kappa_f_b g_b` has lower bound `mu_f_b-3|kappa_f_b|`. A strict floor requires `mu_f_b>3|kappa_f_b|` (or another explicit factor restriction). At `mu_f_b=0.6`, `kappa_f_b=0.20` gives equality.
