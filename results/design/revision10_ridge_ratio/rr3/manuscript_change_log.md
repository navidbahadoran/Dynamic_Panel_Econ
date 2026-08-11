# RR3 manuscript change log

Source: `E:\OneDrive\Desktop\ver7_revision9_Montecarlo_appendix_design.tex`

SHA-256: `2525d019ec5d7b28457585ff57c44670d7fafc6c199c60acc753a76e85070138`

## Rank-selection required

- Replaced the active Revision-9 candidate/IC selector with the frozen blockwise ridge-ratio selector.
- Defined the reportable class `R_max`, pilot-only class `R_max^+`, and one joint cap+1 constrained least-squares spectral pilot.
- Added fitted-value RMS weights, normalized block spectra, `a_NT=1/log(NT)`, the rank-zero anchor and ratio, positive-rank ratios through each reporting cap, smallest-argmin tie handling, and one final literal fixed-rank post-refit.
- Added the pilot-only identification extension, the cap+1 pilot-rate corollary, the blockwise consistency theorem, and the replacement appendix proof.
- Demoted the nuclear path to optional initialization and added the frozen unresolved-pilot failure rule.
- Replaced active IC tuning, Monte Carlo, figure, and empirical diagnostic descriptions with Revision-10 designs.

## Cross-reference/wording required

- Reused `thm:rank_consistency` and `eq:rank_selector` where practical.
- Added dedicated labels for the pilot, normalized spectrum, ridge, zero ratio, positive ratios, pilot objective gap, and pilot-only identification condition.
- Repointed the modular Monte Carlo inputs to Revision-10 table shells.
- Updated the abstract, introduction, related literature, implementation, tuning, conclusion, and appendix headings only where rank-selection wording required it.

## Unrelated

None. DGP formulas, calibration, coefficient envelope, supplied-rank estimator and recovery, Riesz construction, target definitions/applicability, split correction, spatial variance, boundary/interiority rule, and rectangular growth condition are unchanged.

Revision-9 selector passages are retained only inside inactive `\iffalse ... \fi` archival blocks in the new copy; they are not part of the compiled Revision-10 manuscript.
