# Conditioning sigma-field audit

## Conclusion

The paper's field generated only by the four common paths is too small if the theorem treats the true coefficient matrices as fixed after conditioning. A suitable field for a cell `(N,T)` is

\[
\mathcal C_{NT}=\sigma\!\left(
 \{g_{a,t},g_{b,t},g_{h,t},f_{x,t}:t=-50,\ldots,T\},
 \{\lambda_{a,i},\lambda_{b,i},\lambda_{h,i},\lambda_{x,i},
    \sigma_i^2,\sigma_{e,i}^2:i=1,\ldots,N\},
 \{G_i:i=1,\ldots,N\},c_a,c_H,c_\xi
\right).
\]

Here `lambda_a` and `lambda_b` include the DGP-4 group means and their bounded coefficient-loading perturbations. The group labels are deterministic in the code, so including them is harmless and makes the definition portable. The constants `c_H` and `c_xi` are frozen cell constants. The realized `c_a` is measurable with respect to `lambda_a` and `g_a`; listing it explicitly makes the fixed-coefficient claim transparent.

For the code's indexing, the common paths in the field should cover the entire burn-in plus observed horizon used to construct the coefficient matrices and covariates, not just the displayed sample. Equivalently, index them by the code columns `0,...,T+49`, mapped to calendar periods `-49,...,T` after the stated initial values at `-50`.

## Why each enlargement is needed

- `lambda_a`, `lambda_b`, and the DGP-4 loading perturbations are required to make `A_0` and `B_0` fixed conditional on the field.
- `lambda_h` is required to make `H_0` fixed conditional on the field.
- `lambda_x`, `sigma_i^2`, and `sigma_e,i^2` are not needed merely to fix `Theta_0`, but conditioning on them gives a clean conditional mixing/NED statement with fixed unit-specific scales and bounded covariate loadings.
- The deterministic support calculations already give unconditional coefficient envelopes. Conditioning on the unit draws is nevertheless the clean way to align the simulated triangular array with assumptions written conditionally on `C`.

## What should not be included

Do not condition on the time-varying Gaussian disturbance innovations `epsilon_it` or covariate innovations `e_it`. They provide the remaining conditional randomness used for the conditional alpha-mixing and NED arguments. Conditioning on them would trivialize or destroy the intended stochastic structure.

## Findings

| Claim | Four-path field only | Enlarged field |
|---|---:|---:|
| `A_0`, `B_0`, `H_0` fixed conditional on `C` | No | Yes |
| Deterministic coefficient envelopes | Yes, from bounded support, but matrices remain random | Yes |
| Clean conditional alpha-mixing formulation | Incomplete | Yes, retaining time-varying innovations |
| Clean conditional covariate NED formulation | Incomplete | Yes |

Recommended paper action: enlarge the definition of `C_NT` before changing any implementation. No code change is required for this issue.
