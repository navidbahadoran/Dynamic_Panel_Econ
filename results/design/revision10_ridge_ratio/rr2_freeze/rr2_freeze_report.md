# Revision-10 ridge-ratio theory-freeze report

## Outcome

**RR2-FREEZE.** The RR1, RR2, and RR2.5 records now form a complete theory for the Revision-10 blockwise
ridge-ratio rank selector.

## Frozen selector

One joint spectral pilot uses each reporting cap plus one and the unchanged coefficient box. With the existing
block RMS weights, define

`lambda_hat_M,j=(w_M/w_A,1)^2 sigma_j(M_hat^pil)^2/(NT)`,

`a_NT=1/log(NT)`, `lambda_hat_M,0=1`, and

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`.

Each block selects the smallest ratio minimizer over its reporting set. One final literal unpenalized
fixed-rank joint post-refit supplies all coefficients and inference.

## Pilot theorem

The pilot-only class permits ranks through `bar_r_M+1`; differences have rank at most `2bar_r_M+1`. Part (i)
of prediction identification is locally extended to this fixed class with a possibly smaller fixed `c_+>0`.
Uniform score and prediction concentration retain their Revision-9 rates because the rank enlargement changes
only fixed nuclear-norm and entropy constants.

For normalized global objective gap `o_p(zeta_NT)`, the complete basic inequality gives

`||Theta_hat^pil-Theta_0||_F=O_p(sqrt(NT zeta_NT))`

and the same upper order for the maximum block operator error.

## Rank theorem and inference

Weyl perturbation gives order-one normalized signal values and `O_p(zeta_NT)` post-truth values. The maintained
rectangular condition implies `a_NT->0` and `zeta_NT/a_NT->0`. The true ratio alone converges to zero: all
pre-truth competitors remain separated from zero and post-truth competitors converge to one. The anchor
handles rank zero and the extra pilot singular value handles truth at the reporting cap. A finite union gives

`P_0(r_hat=r_0|C)->_p1`.

On this event, the final post-refit is the supplied-true-rank estimator, so existing recovery, Riesz, target,
split-correction, spatial-variance, and feasible-normality results transfer without new target or variance
assumptions.

## Procedure and archival boundaries

The nuclear path is optional initialization only. Revision-9 IC penalties, candidate enumeration, threshold
selection, neighboring IC fits, and local completion are removed from Revision-10 statistical theory while
remaining archived for reproducibility.

An unacceptable pilot is recorded as `rank_selection_numerically_unresolved`; there is no fallback rank and
no primary selected-rank output. Observable diagnostics do not certify the theoretical global objective gap.

## Source and scope

The freeze uses the repository RR1, RR2, and RR2.5 records and the externally supplied manuscript
`E:\OneDrive\Desktop\ver7_revision9_Montecarlo_appendix_design.tex`, SHA-256
`2525d019ec5d7b28457585ff57c44670d7fafc6c199c60acc753a76e85070138`. The requested filename with `(1)` was
not present. No manuscript was substituted or edited.

This package contains theory/design documents only. It ran no Monte Carlo, performed no model fit, implemented
no selector, changed no scientific source, and tuned no constant.
