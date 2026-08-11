# RR5f paper-level decision report

## Evidence interpretation

RR5 and RR5e do not show statistical failure of ridge-ratio rank selection. They supplied no accepted pilot spectrum and therefore evaluated no zero-rank, rank-two, exact-, under-, over-, or mixed-selection behavior. They show that the current constrained least-squares cap+1 pilot is computationally unusable under the maintained numerical definition: 0/360 pilots accepted across two fresh experiments.

The mathematical rank-consistency theorem remains valid conditional on a pilot operator rate and approximate global objective. The computation did not provide a pilot credible under those premises. Mathematical validity and computational feasibility must therefore be reported separately.

## Route conclusions

- **A:** strongest alignment of claims, theory, computation, and evidence. Supplied ranks become inputs; sensitivity is explicit; the selector theorem is removed from this submission.
- **A+:** mathematically clean only if heavily quarantined, but retains substantial notation and referee exposure for an inaccessible auxiliary procedure.
- **B:** objective relaxation affects at most five RR5e pilots; one-start acceptance at most thirty; only stationarity/KKT relaxation attacks the dominant failure, at unacceptable post-hoc risk without a new justification and independent validation.
- **C:** nuclear/full-spectrum and residual-score objects require new estimator-specific rates, normalization, dependence theory, rectangular analysis, tuning freezes, and independent validation. Prior work did not authorize either as a final selector.

The self-tuning route remains closed at ST3 because the missing observable operator-valued variance envelope under conditional spatial dependence is a major theorem, not a constant-selection issue. Revision-9 IC remains a closed NO-GO benchmark: `(0,0,0)` in 24/24, incompatible fixed multipliers across N, and 0/24 thresholded-pilot recovery.

## Contribution and practice

The paper remains coherent as a target-specific inference paper for several separately ranked coefficient matrices in a dynamic panel. Riesz representation, full-panel local-target inference, two-way correction for averages, conditional spatial variance, separate latent spaces, and rectangular asymptotics do not require automatic rank selection.

The practical answer is a prespecified compact rank grid, a transparent baseline supplied-rank specification, neighboring-rank sensitivity, and descriptive fit/singular-value diagnostics. Uncertainty is conditional on each supplied-rank fit. The paper must not claim post-selection, rank-robust, or automatic-rank inference.

## Theory, notation, and literature

Under Route A, remove Assumption 8's pilot-only extension, `R_max^+`, `D_max^+`, cap+1 notation, ridge ratios, and the joint rank-consistency theorem. Reduce Pu et al., Ahn-Horenstein, and recent eigenvalue-ratio discussion to brief context only where useful. Rank-selection literature should not dominate an operative supplied-rank paper.

## Monte Carlo and empirical consequences

Headline Monte Carlo should use supplied true ranks at balanced 100, 200, and 400; label 50 as stress; move rectangular cells to an appendix; and replace rank-recovery tables with supplied-rank sensitivity and inference tables. Empirical work should report a prespecified small grid and conditional inference across neighboring ranks. No new Monte Carlo or manuscript edit occurs in RR5f.

## Decision

**RR5f-RECOMMEND-A.** Route A minimizes unnecessary assumptions and notation, avoids post-hoc numerical tuning, preserves the feasible inference contribution, offers a credible empirical workflow, and provides the shortest scientifically defensible path to submission.
