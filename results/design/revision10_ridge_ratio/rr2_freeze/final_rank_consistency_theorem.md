# Final blockwise ridge-ratio rank-consistency theorem

## Theorem

Suppose `a:stab`, `a:exog`, `a:geometry`, `a:ned`, `a:moments`, `a:signal`, `a:identification`, `a:gram`, and
`a:growth` hold on their maintained domains; the block count and reporting caps are fixed;
`r_0 in R_max`; and a joint pilot satisfies

`max_M||M_hat^pil-M_0||_op=O_p(sqrt(NT zeta_NT))`

conditionally on the common-shock sigma-field. Construct the normalized spectra, deterministic ridge,
anchor, ratios, and blockwise smallest argmins exactly as in `final_selector_definition.md`. Then

`P_0(r_hat=r_0|C)->_p1`.

For the specified joint cap+1 pilot, the operator-rate premise follows from the pilot-only extension of part
(i) of prediction identification to `D_max^+` and the `o_p(zeta_NT)` normalized global objective-gap condition
in `final_cap_plus_one_corollary.md`.

## Proof

The maintained empirical block weights divided by `w_A,1` are jointly bounded above and away from zero with
probability tending to one. Weyl's inequality and the pilot rate imply, uniformly over the fixed blocks and
indices,

- `lambda_hat_M,j=Theta_p(1)`, bounded away from zero, for `j<=r_M,0`;
- `lambda_hat_M,j=O_p(zeta_NT)` for `j>r_M,0`.

The maintained growth condition gives `a_NT->0` and `zeta_NT/a_NT->0` on every admitted rectangular sequence.

If `0<r_M,0<=bar_r_M`, the ratios before truth are bounded away from zero: for suitable fixed spectral bounds
`0<underline_lambda<=overline_lambda`, the anchored competitor is at least `underline_lambda/2` and every
positive pre-truth competitor is at least `underline_lambda/(overline_lambda+1)` with probability tending to
one. At truth,

`R_M,r_M,0=O_p(a_NT)->0`.

After truth, both spectral terms are `o_p(a_NT)`, so `R_M,j->_p1`.

If `r_M,0=0`, the anchored ratio satisfies

`R_M,0={a_NT+o_p(a_NT)}/{1+a_NT}->_p0`,

whereas every `R_M,j`, `j>=1`, converges to one. No zero-rank test is needed.

If `r_M,0=bar_r_M`, the cap+1 pilot permits a genuine estimated `(bar_r_M+1)`st singular value rather than
setting it to zero by the reporting constraint. Its normalized square is `O_p(zeta_NT)`, so the true-cap ratio
still converges to zero.

Thus the true rank is each block's unique asymptotic ratio minimum. With fixed block count,

`P_0(r_hat!=r_0|C)
 <=sum_M P_0(r_hat_M!=r_M,0|C)->_p0`.

No IC, candidate-coverage, penalty-dominance, self-tuning gain, or growing-model argument is used.

## Downstream consequence

On `{r_hat=r_0}`, the final literal post-refit `Theta_hat(r_hat)` is the estimator defined under supplied true
ranks. Hence coefficient and singular-space recovery, empirical Riesz consistency, target expansion, the
two-way split correction, spatial variance consistency, and feasible conditional normality transfer by
asymptotic equivalence. There is no new target assumption or variance theory.
