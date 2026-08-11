# RR3.1 consistency audit

## Frozen selector

The canonical manuscript uses one joint cap+1 constrained least-squares spectral pilot, reportable caps $\bar r_M$, pilot caps $\bar r_M+1$, the frozen fitted-value RMS weights with $w_H=1$ and reference $w_{A,1}$, normalized squared singular values, $a_{NT}=1/\log(NT)$, the rank-zero anchor, genuine cap+1 numerators, smallest-argmin tie handling, separate block ranks, and one final literal fixed-rank post-refit. There is no multiplier, zero-rank test, threshold, candidate enumeration, or IC selector.

## Identification

Assumption 7 has labeled global and local parts. The global condition covers nonzero reportable-class coefficient differences within the box. The local condition covers nonzero directions in the true tangent space. Assumption 8 is named “Pilot-only prediction identification,” refers explicitly to Assumption 7(i), and extends only its global domain to $\calD_{\max}^{+}$ with a possibly smaller fixed $c_+>0$. It does not enlarge Assumption 7(ii) or strengthen stochastic, signal, target, or growth conditions.

## Proof order

Appendix A.7 contains: the fixed-rank extension of uniform concentration to $\calD_{\max}^{+}$; the approximate-objective-gap basic inequality; Frobenius and block-operator pilot rates; scale-weight regularity; Weyl perturbation; normalized signal/noise spectral rates; rectangular ridge separation; positive-rank, zero-rank, and rank-at-cap cases; a finite union; and transfer to the supplied-rank estimator. Observable optimization diagnostics are explicitly not treated as certificates of the unknown global objective gap.

## Language and nuclear-norm disposition

No active occurrence of “nuclear” remains. The selector is described as a spectral pilot followed by a fixed-rank post-refit. The proof uses ordinary operator/Frobenius norm language. Searches found no “theory gate,” “hereditary detectability,” or “coercivity.”

## Monte Carlo and failure accounting

The main text retains rank truths $(1,1,1)$, $(2,1,1)$, and $(1,0,2)$; reporting caps $(3,3,3)$; pilot caps $(4,4,4)$; balanced and rectangular designs; and $R_{\mathrm{inference}}\le R_{\mathrm{point}}\le R_{\mathrm{attempted}}$. The exact unresolved status is `rank_selection_numerically_unresolved`, with no fallback rank. The main and appendix shells contain no numerical results.

## Unchanged scientific content

The reviewed manuscript retains the model, DGP formulas, calibration and pooled-$R^2$ treatment, $B=10$, $c_B=1$, supplied-rank estimator, targets, Riesz construction, two-way split correction, spatial variance, target applicability, stability, sequential exogeneity, spatial geometry, conditional mixing/NED, moments, signal/incoherence, loading/factor Gram condition, target regularity, and rectangular growth condition. No $N=T$, $N/T\to c$, or $N\asymp T$ restriction was introduced.
