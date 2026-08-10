# Existing concentration inventory

All statements below are conditional on the common-shock field and uniform over its regular
realizations.

| Manuscript result | Exact scope | Strength | Order / logarithms | What it supplies | What it does not supply |
|---|---|---|---|---|---|
| Weighted concentration, Lemma `weighted_concentration` | finite deterministic scalar weight collections for conditionally centered, exponentially localized arrays | explicit high-probability Bernstein inequality, but with existential constants `c,C,K_A` | blocking radius `K_A log(NT)` and spatial-time dimension `d_s+1` | scalar finite-union control under mixing/NED | observable studentization or data-dependent optimized directions |
| Matrix consequence, equation `matrix_block_concentration` | one centered N by T array or fixed-size blocks | conditional `O_p` | `sqrt((N+T)log(NT))+M_n log(NT)^((d_s+3)/2)` | operator-score rate | stated tail quantile, sharp constant, or variance-adaptive boundary |
| Score matrix bound, equation `score_matrix_block_bound` | every fixed block M; finite block count | conditional `O_p` | same display with `M_n=b_NT^2`; bounded by `sqrt(NT zeta_NT)` | raw block score operator-norm order | observable block-specific scale |
| Uniform score bound, equation `uniform_score_bound` | continuum of all rank-capped differences `D_max` | conditional `O_p`, uniform | `sqrt(NT zeta_NT)||Delta||` | global basic inequality and finite-lattice overfit rate | a pivotal statistic or probability-indexed cutoff |
| Empirical prediction lower bound, equation `empirical_lower_bound` | all rank-capped differences | with conditional probability tending to one, uniform | `(c/2)||Delta||^2-C NT zeta_NT` | global restricted curvature outside the recovery radius | an observable `c` or a block Schur curvature estimator |
| True-tangent Gram convergence | true joint tangent space | `o_p(1)` operator convergence | displayed block rates in the proof | local fixed-rank curvature at truth | uniform curvature on estimated mixed-state normal cones |
| Rank proof overfit comparison | finite rank lattice | conditional `O_p` | normalized improvement `O_p(zeta_NT)`; RSS improvement `O_p(NT zeta_NT)` | qualitative overfit separation | observable rate-sharp envelope |
| Spatial HAC theorem | fixed finite Riesz-weighted scalar targets | consistency | bandwidth `h_N=ceil(c_h log(NT))` | scalar target variance under spatial dependence | growing N by N block-score covariance estimation |

Here

`b_NT=(NT)^(1/(8+eta)) log(NT)` and
`zeta_NT=b_NT^2 log(NT)^(d_s+2)(1/N+1/T)`.

The score result is stronger than a pointwise expectation statement: it is uniform conditional
`O_p` over the entire admissible rank-capped coefficient-difference set. Nevertheless `O_p`
means that for each probability tolerance some unknown finite constant works. It cannot be turned
into an exact data-computable inclusion boundary by setting that constant to one. The manuscript's
finite union argument removes any need for a growing model-count penalty, but it does not identify
the conditional operator-score quantile.
