# RR5f recommendation

## Recommended route

**Route A: supplied-rank main paper.**

This recommendation prioritizes scientific defensibility over preserving the Revision-10 selector work. The supplied-rank estimator and inference theory are intact and computationally operative. The cap+1 pilot failed its maintained numerical gate in 360/360 fresh attempts, so presenting automatic rank choice as part of the operative method would not be credible. Relaxing the gate after these outcomes creates severe post-hoc risk; replacing the pilot requires a new theory project. Retaining the theorem as auxiliary material preserves substantial assumptions and notation for a method that produces no reported estimate.

Route A yields the cleanest paper: inference conditional on separately supplied ranks, an ex-ante empirical rank grid, and transparent nearby-rank sensitivity. It removes the unnecessary pilot identification extension and keeps the focus on target-specific inference.

## Exact manuscript-change outline

No manuscript is edited in RR5f. The next separately authorized manuscript phase should make these changes:

1. **Introduction:** remove automatic rank selection from the contribution list and abstract-level claims. State clearly that ranks are specification inputs and that the empirical analysis reports prespecified sensitivity.
2. **Section 4 estimator:** present the supplied-rank joint estimator as the operative estimator. Remove the cap+1 pilot, spectral normalization, ridge-ratio definition, zero-rank anchor, and selected-rank final post-refit from the operative algorithm.
3. **Assumptions:** delete the pilot-only Assumption 8 / `D_max^+` clause; retain the original supplied-rank identification domain and all inference assumptions.
4. **Rank-selection subsection:** replace it with “rank specification and sensitivity.” Define the prespecified small rank grid, baseline rationale, neighboring specifications, and conditional interpretation.
5. **Theorem 5:** remove the joint ridge-ratio consistency theorem from the submitted paper and renumber later results. Do not replace it with an informal selection claim.
6. **Appendix A.7:** remove the cap+1 pilot and ratio-consistency proof and all cross-references. Preserve the material in repository history, not in the submitted appendix.
7. **Implementation:** describe fixed supplied-rank fitting, numerical acceptance, Riesz construction, split correction, and variance estimation. State that automatic rank selection is outside the operative implementation.
8. **Monte Carlo:** use supplied true ranks for headline inference; headline balanced cells 100, 200, and 400; retain 50 as stress; move rectangular cells to an appendix; replace rank-selection tables with supplied-rank sensitivity tables.
9. **Empirical application:** prespecify a compact rank grid, give the baseline specification rationale, report neighboring-rank estimates and fit/singular-value diagnostics, and avoid rank-robust language.
10. **References:** shorten Pu et al., Ahn-Horenstein, and related eigenvalue-ratio discussion to contextual remarks or remove citations used only by the deleted theorem. Keep focus on low-rank estimation and inference.
11. **Tables:** remove or archive ridge-ratio recovery, cap-selection, and block-ratio tables. Add supplied-rank sensitivity, numerical-validity, fit, and target-inference tables.

## Engineering preservation

Retain RR5d resumability, atomic checkpoints, fingerprints, deterministic aggregation, worker-local caching, performance instrumentation, and benchmark infrastructure. These are method-neutral simulation-engineering improvements.
