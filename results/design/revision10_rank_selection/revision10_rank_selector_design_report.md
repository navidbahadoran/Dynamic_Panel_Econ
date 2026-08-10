# Revision-10 self-tuning rank-selector design report

## Outcome

Revision 9 is frozen as a closed NO-GO benchmark. Revision 10 has an ST2 theory design, not an implemented selector. No simulation, model fit, scientific-source edit, manuscript edit, or outcome-driven tuning occurred.

## Frozen evidence

The locked IC chose zero rank 24/24 despite truth being in every candidate set. Fixed-c requirements conflict across N=50 and N=100 (CASE C). The thresholded cap pilot was exact 0/24 and over-ranked 23/24 (CASE P3). A nuclear proposal equaled truth 24/24, so screening worked but is not promoted to final selection. N=100 supplied-rank execution/inference was clean; N=50 remains a small-sample stress design with narrow KKT and split-interiority failures. Validated outer parallelism recommends `--n-jobs 12` on this machine.

## Manuscript architecture and preserved layers

The external author-supplied Revision-9 manuscript was read literally. It defines the rank vector, finite cap lattice, weighted nuclear path with gradient-based lambda_max, one rank-cap pilot, tau thresholding, R_base, neighbor-enriched R_cand, unpenalized joint fixed-rank post-refits, Q_hat, dimension d, additive IC, and IC-based local completion. The candidate-coverage proposition uses an included pilot's operator-norm accuracy and Weyl separation; the rank theorem combines coverage with IC under/overfit separation.

Revision 10 keeps screening, caps, cap pilot only as screening/coverage machinery, tau only for screening, candidates, literal box constraint, all joint post-refits, and all fixed-rank inference. It provisionally removes only the additive IC from the final statistical decision and converts local completion into statewise neighbor generation.

## Primitive gain and preferred design

For `r+=r+e_M`, define `Delta_M=RSS(r)-RSS(r+)` and

`Delta_d_M=N+T-2r_M-1`.

The preferred working statistic is `G_M=[Delta_M/Delta_d_M]/s_hat_M^2(r+)`. Start at zero. At each step evaluate all admissible additions, use matrix-specific scales for A^(ell), B^(k), and H, add the largest standardized passing coordinate with a fixed tie rule, and stop when none passes. Zero is valid; capped coordinates are skipped. This is order-defined and finite, but not yet frozen because the exact self-normalizer and boundary lack a panel theorem.

The matrix scales must differ because the block score is multiplied respectively by lagged outcomes, block covariates, or one. A global maximum envelope is possible but unnecessarily conservative and does not solve calibration.

## Bing-Wegkamp adaptation

Bing-Wegkamp contribute the principle of a conservative lower-rank start, residual/noise self-scaling, a monotone-upward update that preserves no-overfit control, and a fixed point stop. Their literal STRS uses one projected reduced-rank regression, closed-form truncated SVD, scalar rank, and iid Gaussian calibration. Those identities fail for jointly refitted dynamic-panel rank vectors under conditional temporal/spatial dependence; their formula is therefore not copied.

## Separation and assumptions

Underfit: strong order-sqrt(NT) singular values plus prediction identification should make at least one omitted block direction reduce RSS by order NT, dominating the existing stochastic rate. Overfit: existing score/concentration results limit an unnecessary optimized improvement to `O_p(zeta_NT)` over a finite cap lattice. To turn these rate statements into an observable sequential rule requires exactly two candidate new primitive conditions: uniform block self-normalizer validity and hereditary one-step detectability. All ten named maintained assumptions remain substantively unchanged.

The additive penalty sequence, multiplier, penalty-rate condition, IC gap/sensitivities, IC local completion, and penalty-dominates-zeta proof could disappear. Numerical objective accuracy would instead be measured against local gain margins. Nothing is removed until a theorem exists.

## Transfer and gate

If `P(r_hat_ST=r_0|C)->_p 1` is proved, the supplied-rank recovery, target expansion, split correction, Riesz construction, spatial variance, and selected-rank inference transfer unchanged on the recovery event. The unresolved self-normalizer and mixed-state induction make the present decision ST2. Implementation, theorem drafting, and new statistical runs remain unauthorized.
