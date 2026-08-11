# RR4 selector formula audit

| Frozen item | Literal implementation |
|---|---|
| Pilot class | One joint `fit_fixed_rank` call per numerical start with bounds `bar_r_M+1`; represented matrix rank is at most the factor-column count. |
| Box | Existing literal `||Theta||max <= B` interior ALS / alternating exact linear-box-QP solver; no completed-matrix clipping. |
| `w_A,l` | `sqrt(mean(y_lag_l**2))` on the full estimation sample. |
| `w_B,k` | `sqrt(mean(x_k**2))` on the full estimation sample. |
| `w_H` | Literal `1.0`. |
| Reference | First lag weight `w_A1`; nonfinite or nonpositive is unresolved and is never floored. |
| Spectrum | `(w_M/w_A1)^2 * sigma_j(M_pil)^2/(N*T)` from fitted pilot matrices through cap+1. |
| Ridge | Exactly `1/log(N*T)`. |
| Rank zero | `R_0=(lambda_1+a)/(1+a)` with anchor one. |
| Positive rank | `R_j=(lambda_{j+1}+a)/(lambda_j+a)`. |
| Rank at cap | Uses the actual stored cap+1 pilot singular value. |
| Tie | `numpy.argmin`, hence the first/smallest exact minimizer; no near-tie tolerance. |
| Assembly | Independent ordered A blocks, B blocks, and H. |
| Reported estimator | One selected-vector call to the maintained literal fixed-rank multistart wrapper. |

The selector contains no IC, multiplier, candidate enumeration, threshold, zero test, or path
rank. The ratio gap is serialized only as a diagnostic.
