# Authoritative Revision-10 ridge-ratio theory freeze

This protocol is activated by the combined RR1, RR2, and RR2.5 theory record. The following quantities and
rules may not change after RR5 finite-sample results are observed without a new paper-level theory decision.

1. **Reporting caps.** Each fixed block reports ranks only in `{0,...,bar_r_M}`. The fixed block count and
   caps do not grow.
2. **Pilot caps.** Every spectral-pilot cap is exactly `bar_r_M+1`. The extra rank is never reportable.
3. **Pilot fit.** Use one joint least-squares spectral pilot with all cap+1 constraints and the unchanged
   canonical-unit box `||Theta||_max<=B`.
4. **Pilot identification.** For the pilot result only, part (i) of prediction identification holds on
   `D_max^+` with a fixed `c_+>0`. Supplied-rank theory retains the original `R_max` domain.
5. **Pilot objective condition.** The normalized global enlarged-cap objective gap is `o_p(zeta_NT)`.
6. **Block weights.** Freeze the empirical RMS lag and covariate weights and `w_H=1` exactly as stated in
   `final_selector_definition.md`.
7. **Spectral normalization.** Freeze
   `lambda_hat_M,j=(w_M/w_A,1)^2 sigma_j(M_hat^pil)^2/(NT)`.
8. **Ridge.** Freeze `a_NT=1/log(NT)` with no multiplier or sensitivity version in the statistical rule.
9. **Rank-zero anchor.** Freeze `lambda_hat_M,0=1`; do not add a test or cutoff.
10. **Ratio formula.** Freeze
    `R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)` for every reportable index, using the anchor at zero.
11. **Tie rule.** Exact numerical ties select the smaller rank.
12. **Blockwise selection.** Select each block separately and assemble the full vector; do not impose a
    common rank.
13. **Cap treatment.** The true-cap ratio uses the pilot's permitted `(bar_r_M+1)`st singular value; it is not
    set to zero by reporting truncation.
14. **Rank-zero treatment.** Zero competes through the deterministic anchor only.
15. **Nuclear path.** It is not statistical selection. It may only initialize the pilot or final post-refit.
16. **Final post-refit.** Compute one literal unpenalized constrained joint fit at `r_hat`; all inference uses
    that fit and never the pilot.
17. **Numerical diagnostics.** Freeze feasibility, multistart objectives and residuals, objective stability,
    box activity, numerical rank/spectra, start provenance, runtime, and termination diagnostics. Numerical
    tolerance values must be prespecified before RR5.
18. **Failure rule.** An unacceptable pilot gives `rank_selection_numerically_unresolved`, no fallback rank,
    and no primary selected-rank point or inference output.
19. **Rectangular theory.** Retain the exact `1/N+1/T` growth condition and impose no aspect-ratio restriction.
20. **Revision-9 disposition.** IC penalties, candidates, thresholding, and local completion have no
    Revision-10 statistical role; historical code and evidence remain for reproducibility.

No pilot cap, normalization, ridge, anchor, ratio, cutoff, fallback, or diagnostic acceptance rule may be
changed in response to RR5 results without reopening the paper-level decision and documenting a new proof.
