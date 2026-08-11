# RR2 formal selector, theorem, and theory-freeze report

## Decision

**RR2-PARTIAL.** Revision 10 has a complete scale-equivariant generic ridge-ratio theorem, but cannot yet
freeze the requested cap+1 implementation under the unchanged assumptions. This phase
created theory/design documents only. It did not run Monte Carlo or model fitting, alter scientific code, or
edit a manuscript.

## Sources of truth

The complete RR1 directory was used. The latest supplied manuscript is external to the repository:

`E:\OneDrive\Desktop\ver7_revision9_Montecarlo_appendix_design.tex`

SHA-256: `2525d019ec5d7b28457585ff57c44670d7fafc6c199c60acc753a76e85070138`.

The requested suffixed file `ver7_revision9_Montecarlo_appendix_design(1).tex` was absent; this fact is
recorded rather than silently substituting an older tracked manuscript.

## Exact statistical procedure

The proposed procedure would compute one joint least-squares spectral pilot subject to each fixed reporting cap plus one and the unchanged
literal coefficient box. With the existing screening weights and `s_NT=w_A,1`, form

`lambda_hat_M,j=(w_M/s_NT)^2 sigma_j(M_hat^pil)^2/(NT)`.

Set `lambda_hat_M,0=1`, `a_NT=1/log(NT)`, and

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`.

Select the smallest minimizer over `0,...,bar_r_M` independently by block, assemble the rank vector, and
compute one final literal fixed-rank unpenalized joint post-refit. All inference uses that post-refit.

## Scale resolution

Raw coefficient spectra fail finite unit equivariance. The frozen spectrum is dimensionless and uses only
already-maintained weights. Covariate scaling multiplies its screening weight and divides its coefficient by
the same magnitude; a common outcome-unit change is canceled by `w_A,1`. The scale is not calibrated to
rank or inferential outcomes. Exact finite constrained optimization comparisons use stored canonical units,
because a fixed numerical box in two different coefficient units would be two different parameter spaces.

## Theory

The theorem assumes only the existing assumptions, fixed admissible truth, and the generic operator-rate
pilot. Weyl and strong signal give normalized positive values bounded above and away from zero and post-truth
values `O_p(zeta_NT)`. The maintained rectangular condition gives

`zeta_NT/a_NT=G_NT/{(NT)^(2/(8+eta)) log(NT)}->0`.

For positive truth, pre-truth ratios stay explicitly bounded away from zero, the truth ratio tends to zero,
and post-truth ratios tend to one. For zero truth, only the anchored ratio tends to zero. The cap+1 pilot
provides a non-mechanical next singular value when truth is at the reporting cap. Each block therefore has a
unique asymptotic minimum, and a finite union proves joint vector recovery.

## Unresolved cap+1 corollary

The joint cap+1 pilot does **not** presently satisfy a derived generic rate merely from an objective gap.
Revision 9 states prediction identification only for candidates inside `R_max`. At cap truth, the reporting
class controls differences through rank `2bar_r_M`, whereas a cap+1 candidate can differ by rank
`2bar_r_M+1`. Fixedness extends covering and score bounds but not restricted curvature. A null extra-rank
direction can preserve the Revision-9 condition and still give an enlarged optimizer a strong spurious
singular value.

The `o_p(zeta_NT)` normalized objective gap is the right optimization scale **conditional on** enlarged-class
curvature, but that curvature would be a substantive assumption change. RR2 adopts no such change. The
alternatives are to strengthen identification, impose the generic rate directly, or develop another
extra-spectrum construction.

## Numerical and downstream consequences

Multistart agreement, stability, KKT/stationarity residuals, feasibility, box activity, numerical rank, and
runtime are diagnostics only; they do not certify the global infimum. Failure produces an unresolved rank and
suppresses selected-rank inference.

On exact rank recovery the final post-refit is the supplied-rank estimator. Recovery, Riesz consistency,
target expansion, two-way correction, spatial variance, and feasible normality transfer by asymptotic
equivalence. The nuclear path is optional initialization only. Revision-9 IC penalties, candidates, local
completion, and sensitivity multipliers are removed from the new statistical theory and retained only where
needed to reproduce Revision 9.

## Proposed freeze status

The normalization, ridge, anchor, ratios, generic theorem, and rectangular proof are ready. The overall
protocol is not activated, no implementation is authorized, and RR3/RR4 must wait for a paper-level resolution
of the cap+1 identification gap.
