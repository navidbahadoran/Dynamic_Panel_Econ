# Candidate self-tuning formulations

For every valid post-refit define `RSS(r)=NT Q_hat_NT(r)`. For an admissible block increment `r+=r+e_M`,

`Delta_M(r)=RSS(r)-RSS(r+)`.

Because the M contribution to dimension is `r_M(N+T-r_M)`,

`Delta_d_M(r)=(r_M+1)(N+T-r_M-1)-r_M(N+T-r_M)=N+T-2r_M-1`.

No locked outcome was used in these definitions.

Let `J=P+K+1`, and let `alpha_NT` denote a deterministic sequence tending to zero that is fixed by a future theorem, not by Revision-9 outcomes. Whenever `qhat_M(r;alpha)` appears below, it means an estimated simultaneous upper `1-alpha/J` quantile of the optimized null gain per added dimension for block M. Its existence and construction are precisely the unresolved theory question; no numerical `alpha_NT` is frozen here.

## F1: vector residual-denominator criterion

Formula: `C_lambda(r)=RSS(r)/(NT-lambda d(r))`, on positive denominators. Set `sigmahat_t^2=RSS(r^(t))/(NT-d(r^(t)))` and

`lambda_(t+1)=max_M qhat_M(r^(t);alpha_NT)/sigmahat_t^2`,

with the maximum over feasible additions. Start with `r^(0)=argmin C_(lambda_0)(r)` over the screened candidate set, where `lambda_0` uses the same display at r=0 and an extra factor two to be conservative. Update by minimizing `C_(lambda_(t+1))` over candidates coordinatewise no smaller than `r^(t)`. Stop at an unchanged vector. Rank zero is admissible and terminal if selected at the initial step; caps and positive denominators restrict the search. A fixed dimension-first lexicographic rule resolves ties, but incomparable vectors make the result sensitive to the chosen partial-order completion.

This most closely resembles STRS but lacks a canonical panel analogue of R/U, a unique minimizer, or scalar-path monotonicity. Searching a rank lattice can require the full Cartesian set. It would require uniform quantile validity plus new lattice-shape/monotonicity assumptions, likely exceeding the two-condition budget. Compatibility with cached post-refits is moderate; proof difficulty and computational cost are very high.

## F2: one-coordinate self-normalized post-refit gain (preferred working design)

Formula:

`G_M(r)={Delta_M(r)/Delta_d_M(r)}/s_hat_M^2(r+)`.

Define `s_hat_M^2(r+)=qhat_M(r;alpha_NT)`, so the normalized inclusion boundary is exactly one. Start at the zero vector. At a state r, jointly fit every admissible `r+e_M` not already cached and add a coordinate only if `G_M(r)>1`. Evaluate all coordinates simultaneously; select the largest `G_M-1`, with a deterministic manuscript-block-order tie rule. Recompute all eligible statistics after the joint refit and stop when no addition passes. Zero is terminal when no first addition passes; capped blocks are omitted. Because all neighbors are evaluated before a move and the tie rule is fixed, the algorithm is deterministic, although the proof must show that its path cannot make a false early addition.

The scale must estimate an upper envelope for the **optimized extra-rank noise gain**, not merely residual variance. A tentative construction uses the fitted residual and the normal-space block score

`S_M(r)=P_(U_M perp) [Z_M elementwise residual(r)] P_(V_M perp)`,

where Z_M is the relevant lagged outcome, covariate, or one for H. A conditional spatial-HAC/mixing-based multiplier or analytic bound would estimate the upper quantile of `||S_M||_op^2/Delta_d_M`; this is naturally matrix-specific. Ordinary variance studentization is insufficient because a fixed threshold would leave a nonvanishing false-addition probability.

This uses existing neighbor post-refits and costs at most J fits per step and at most `sum_M cap_M` steps, with cache reuse. It requires the two localized conditions in the assumption map: uniform self-normalizer validity and hereditary one-step detectability. Its proof can reuse the existing order-one underfit and `O_p(zeta_NT)` overfit separation. The unresolved ingredient is a tuning-free, uniformly conservative estimate of the data-adaptive operator-score envelope under the maintained conditional dependence. Proof difficulty is high but localized; compatibility with existing post-refits is high.

## F3: cross-fitted validation gain

On fold a choose the rank-one normal direction for M, and on the held-out fold b compute its prediction-loss improvement `Delta_(M,a->b)`. Swap the folds and set

`Delta_M^CF={Delta_(M,1->2)+Delta_(M,2->1)}/2`,

`G_M^CF=Delta_M^CF/{Delta_d_M s_hat_M,CF^2}`,

where `s_hat_M,CF^2` is the spatial-HAC standard-error scale inflated to its simultaneous `1-alpha_NT/J` upper quantile. Start at zero, evaluate all blocks, add the largest statistic exceeding one, refit/re-split at the new state, and stop when none exceeds one. Zero and caps are handled as in F2; fixed folds and block tie order make the result deterministic.

Cross-fitting reduces selection bias and gives a cleaner scalar CLT conditional on the discovery half, but an asymptotically vanishing false-addition level still requires a prespecified tail sequence and a uniform moderate-deviation result. It also adds selection splits, may weaken signal, complicates preservation of external lags, and changes the paper's one-full-panel-selection architecture. It requires at least uniform cross-fitted tail validity and split-level hereditary detectability. Proof difficulty and computation are high, and compatibility with existing full-panel candidate fits is low.

## Comparison and provisional choice

F2 best preserves the existing estimator and candidate/refit cache and respects block-specific score geometry. Matrix-specific scales are necessary: A^(ell) multiplies lagged outcomes, B^(k) multiplies its covariate, and H multiplies one; their conditional second moments and dependence envelopes are not interchangeable. A global scale could be conservatively the maximum block envelope but would lose power and still needs the same missing theorem.

F2 is a working design, not a frozen selector. No numerical boundary, quantile, bandwidth, or multiplier is selected in Phase 1.
